#!/usr/bin/env bash

# --------------------------
# Default values
# --------------------------

THREADS=6

# --------------------------
# Helper functions
# --------------------------

usage() {
    cat <<EOF
Welcome to bulk_mrnaseq_preprocessing!

Description:

  Bulk paired-end mRNA-seq preprocessing pipeline:
    1) fastp trimming
    2) FastQC + MultiQC
    3) Kallisto pseudoalignment

Usage:

  ./bulk_mrnaseq_preprocessing.sh -i <input_dir> -o <organism> [-t <threads>]

Arguments:

  -i, --input     : input directory containing raw FASTQ files
  -o, --org       : organism (mouse, human, rat)
  -t, --threads   : number of CPU threads (default: 6)

EOF
    exit 0
}

check_fastq () {
    local dir="$1"

    if ! compgen -G "$dir"/*.fastq.gz >/dev/null 2>&1; then
        echo "No FASTQ.GZ files found in $dir"
        exit 1
    fi
}

install_package () {
    local package="$1"

    if ! command -v "$package" >/dev/null 2>&1; then
        echo "$package not found, installing..."
        sudo apt update && sudo apt install -y "$package"
        echo "$package installed successfully."
    fi
}

# --------------------------
# Parse named options
# --------------------------

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
        RAW_DIR="$2"
        shift; shift
        ;;
        -o|--org)
        ORGANISM="$2"
        shift; shift
        ;;
        -t|--threads)
        THREADS="$2"
        shift; shift
        ;;
        -h|--help)
        usage
        ;;
        *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# --------------------------
# Check mandatory args
# --------------------------

if [[ -z "$RAW_DIR" || -z "$ORGANISM" ]]; then
    echo "ERROR: --input and --org are required"
    exit 1
fi

if [ ! -d "$RAW_DIR" ]; then
    echo "ERROR: Input directory '$RAW_DIR' not found."
    exit 1
fi

case "$ORGANISM" in
    mouse|human|rat)
        echo "Proceeding"
        ;;
    *)
        echo "ERROR: Invalid organism '$ORGANISM'. Valid options: mouse, human, rat"
        exit 1
        ;;
esac

echo "Input dir: $RAW_DIR"
echo "Threads: $THREADS"
echo "Organism: $ORGANISM"



# --------------------------
# 1) fastp trimming
# --------------------------

fastp_trimming () {
    local RAW_DIR="$1"
    local THREADS="$2"
    local OUTPUT_DIR="$RAW_DIR/../1_trimmed"

    # Skip if already trimmed
    if [ -d "$OUTPUT_DIR" ] && ls "$OUTPUT_DIR"/*_R1_001.trimmed.fastq.gz >/dev/null 2>&1; then
        echo "[SKIP] fastp trimming already done."
        return
    fi

    mkdir -p "$OUTPUT_DIR"

    # Install fastp if missing
    install_package fastp

    # Check for FASTQ files
    check_fastq "$RAW_DIR"

    # Process paired-end reads
    echo "Run fastp..."
    for r1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
        # Resolve corresponding R2
        local r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

        if [[ ! -f "$r2" ]]; then
            echo "WARNING: Missing paired file for $r1"
            exit 1
        fi

        # Extract sample name
        local sample=$(basename "$r1" "_R1_001.fastq.gz")

        fastp \
            --in1 "$r1" --in2 "$r2" \
            --out1 "$OUTPUT_DIR/${sample}_R1_001.trimmed.gz" \
            --out2 "$OUTPUT_DIR/${sample}_R2_001.trimmed.gz" \
            --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
            --detect_adapter_for_pe \
            --thread "$THREADS" \
            --html "$OUTPUT_DIR/${sample}.fastp.html" \
            --json "$OUTPUT_DIR/${sample}.fastp.json"
    done
}

# --------------------------
# 2) FastQC + MultiQC
# --------------------------

fastqc_multiqc () {
    local INPUT_DIR="$1"
    local THREADS="$2"
    local OUTPUT_DIR="$INPUT_DIR/../2_fastqc"

    mkdir -p "$OUTPUT_DIR"

    # Skip step if MultiQC already exists
    if [ -d "$OUTPUT_DIR" ] && [ -f "$OUTPUT_DIR/multiqc_report.html" ]; then
        echo "[SKIP] FastQC + MultiQC already done."
        return
    fi

    # Install fastqc if missing
    install_package fastqc

    # Check for FASTQ files
    check_fastq "$INPUT_DIR"

    # Run FastQC in parallel
    echo "Running FastQC..."
    find "$INPUT_DIR" -maxdepth 1 -name "*.trimmed.gz" | parallel -j "$THREADS" fastqc {} --outdir "$OUTPUT_DIR"

    # Run MultiQC
    multiqc "$OUTPUT_DIR" --outdir "$OUTPUT_DIR"
}

# --------------------------
# 3) Kallisto pseudoalignment
# --------------------------

kallisto_pseudoalignment () {
    local INPUT_DIR="$1"
    local THREADS="$2"
    local ORGANISM="$3"
    
    local REF_DIR="$INPUT_DIR/../3_reference"
    local OUTPUT_DIR="$INPUT_DIR/../4_aligned"

    mkdir -p "$OUTPUT_DIR" "$REF_DIR"

    # Skip kallisto if quantification already exists for all samples
    if [ -f "$REF_INDEX" ] && ls "$OUTPUT_DIR"/*/abundance.tsv >/dev/null 2>&1; then
        echo "[SKIP] Kallisto pseudoalignment already done."
        return
    fi

    # Install kallisto if missing
    install_package kallisto

    # Select reference organism
    case "$ORGANISM" in
        mouse)
            REF_URL="https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
            IDX_NAME="mouse_transcriptome.idx"
            ;;
        
        human)
            REF_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
            IDX_NAME="human_transcriptome.idx"
            ;;
        
        rat)
            REF_URL="https://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
            IDX_NAME="rat_transcriptome.idx"
            ;;
    esac

    REF_FASTA="$REF_DIR/$(basename "$REF_URL")"
    REF_INDEX="$REF_DIR/$IDX_NAME"

    # Download reference file if missing
    if [[ ! -f "$REF_FASTA" ]]; then
        echo "Downloading $ORGANISM transcriptome..."
        wget -c "$REF_URL" -P "$REF_DIR"
    fi
    
    # Build index file if missing
    if [[ ! -f "$REF_INDEX" ]]; then
        echo "Building index file..."
        kallisto index -i "$REF_INDEX" "$REF_FASTA"
    fi

    # Run pseudoalignment
    for r1 in "$INPUT_DIR"/*_R1_001.trimmed.fastq.gz; do
        sample=$(basename "$r1" "_R1_001.trimmed.fastq.gz")
        r2="$INPUT_DIR/${sample}_R2_001.trimmed.fastq.gz"

        echo "Running Kallisto on $sample..."

        kallisto quant \
        -i "$REF_INDEX" \
        -o "$OUTPUT_DIR/$sample" \
        --bias \
        --threads="$THREADS" \
        --rf-stranded \
        "$r1" "$r2"
    done

    # Create report file
    # Install jq if missing
    install_package jq

    report_file="$(date +"%d%m%Y")_kallisto_pseudoalignment_report.txt"

    # Create header
    echo "Sample Name | Pseudoalignment Rate (%)" > $report_file
    echo "------------------------------------" >> $report_file

    # Loop through each output folder
    for sample in "$OUTPUT_DIR"/*; do
        sample_name=$(basename "$sample")
        pseudoaligned=$(jq .p_pseudoaligned "$sample/run_info.json")

        # Print to console
        echo "$sample_name | $pseudoaligned% pseudoaligned"

        # Append results to report file
        echo "$sample_name | $pseudoaligned%" >> "$report_file"
    done
}

# --------------------------
# 4. Run pipeline
# --------------------------

bulk_mrnaseq_preprocessing() {
    fastp_trimming "$RAW_DIR" "$THREADS"
    TRIMMED_DIR="$RAW_DIR/../trimmed"
    fastqc_multiqc "$TRIMMED_DIR" "$THREADS"
    kallisto_pseudoalignment "$TRIMMED_DIR" "$THREADS" "$ORGANISM"
    echo "### Data processing completed successfully"
}

# --------------------------
# Execute pipeline if run directly
# --------------------------
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    bulk_mrnaseq_preprocessing
fi