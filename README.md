# scRNAseq-Seurat-analysis-pipeline
This repository demonstrates a reproducible workflow for analysing scRNA-seq data using R and the Seurat framework.
###  Analysis Pipeline

This project follows a structured single-cell transcriptomic workflow.

| Step | Script | Purpose |
| :--- | :--- | :--- |
| **1** | `01_preprocessing.R` | Data ingestion & Gene ID mapping |
| **2** | `02_qc_filtering.R` | Quality control & cell filtering |
| **3** | `03_integration.R` | Cross-species integration (Human/Canine) |
| **4** | `04_annotation.R` | Cell-type identification & clustering |
| **5** | `05_pathway_analysis.R` | Biological pathway enrichment (KEGG/GO) |
| **6** | `06_STRING_analysis.R` | Protein-protein interaction networks |
| **7** | `07_cellchat_analysis.R` | Intercellular communication analysis |
| **8** | `08_infercnv.R` | Chromosomal CNV estimation |
| **9** | `09_trajectory_analysis.R` | Pseudotime developmental modeling |
