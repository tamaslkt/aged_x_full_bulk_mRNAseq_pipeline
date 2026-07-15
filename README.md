# Mouse Heart Bulk mRNA-seq Pipeline

This repository contains an end-to-end bulk transcriptomics workflow, featuring custom bash preprocessing scripts and structured R/Quarto analysis pipelines. It spans raw sequence quality control to advanced functional enrichment, transcription factor activity modeling, and network analysis.

## About the Data
The workflow processes Illumina-platform mRNA-seq data obtained from mouse heart tissue across three conditions:
* **Young Control** (YC, n = 4)
* **Aged Control** (AC, n = 4)
* **Aged Treated** (Treated, n = 10)

*Note: Because this data is currently under preparation for peer-reviewed publication, the name of the specific treatment is kept anonymous.*

---

## Directory Structure
```text
├── analysis/
│   ├── 1_QC.qmd                          # Data normalization, PCA, and exploratory plots
│   ├── 2_differential_gene_expression.qmd # DESeq2 modeling, custom volcanoes, & MA plots
│   ├── 3a_functional_enrichment.qmd      # ORA and GSEA (KEGG, GO, Reactome, Hallmark)
│   └── 3b_transcription_factor.qmd       # Upstream regulator & TF activity profiling
├── scripts/
│   └── bulk_mrnaseq_preprocessing.sh     # Bash pipeline (fastp, FastQC, Kallisto)
└── README.md
