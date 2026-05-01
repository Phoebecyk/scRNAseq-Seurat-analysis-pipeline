# scRNAseq-Seurat-analysis-pipeline
This repository demonstrates a reproducible workflow for analysing scRNA-seq data using R and the Seurat framework.

### Analysis Pipeline

This project follows a structured single-cell transcriptomic workflow with a sequential core pipeline followed by parallel downstream analyses.

#### Core Pipeline (sequential)

| Step | Script | Purpose |
| :--- | :--- | :--- |
| **1** | `01_preprocessing.R` | Data ingestion & gene ID mapping |
| **2** | `02_qc_filtering.R` | Quality control & cell filtering |
| **3** | `03_integration.R` | Cross-species integration (Human/Canine) |
| **4** | `04_annotation.R` | Cell-type identification & clustering |
| **5** | `05_infercnv.R` | Chromosomal CNV estimation to confirm cancer/metastatic cell populations |
| **6** | `06_DEG_analysis.R` | Differentially expressed genes analysis |

#### Downstream Analyses (parallel, run after core pipeline)

| Script | Purpose |
| :--- | :--- |
| `07_pathway_analysis.R` | Biological pathway enrichment (KEGG/GO) |
| `08_STRING_analysis.R` | Protein-protein interaction networks |
| `09_cellchat_analysis.R` | Intercellular communication analysis |
| `10_trajectory_analysis.R` | Pseudotime developmental modeling |
