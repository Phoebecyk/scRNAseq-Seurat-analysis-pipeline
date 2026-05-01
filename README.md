# scRNAseq Seurat Analysis Pipeline

A reproducible single-cell RNA-seq workflow for canine pheochromocytoma (PCC) using R and Seurat. Human fetal adrenal samples (GEO: GSE137804) are used as a reference in the absence of normal canine adrenal data.

---

## Requirements

- R ≥ 4.3
- Package dependencies managed via [`renv`](https://rstudio.github.io/renv/) — see `renv.lock`

---

## Setup

```r
# 1. Restore the package library from the lockfile
renv::restore()

# 2. Open config.R and verify paths and parameters before running any script
# scripts/config.R
```

---

## Data Layout

Place raw data under `data/` before running the pipeline:

```
data/
├── adrenal_scRNA_seq_data/
│   ├── GSE137804_RAW/
│   │   ├── GSM4088787_F106_gene_cell_exprs_table.xls.gz   # Human fetal (normal ref)
│   │   └── GSM4088788_F107_gene_cell_exprs_table.xls.gz
│   ├── PCCVanDijk/filtered_feature_bc_matrix/             # Canine PCC samples
│   ├── PCCDakota/filtered_feature_bc_matrix/
│   └── PCCRose/filtered_feature_bc_matrix/
```

---

## Configuration

All shared parameters live in `scripts/config.R`. Source it at the top of any script:

```r
source(here("scripts/config.R"))
```

Key parameters:

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `ORGANISM_KEGG` | `"cfa"` | KEGG organism code (Canis familiaris) |
| `MIN_FEATURES` | `500` | Minimum genes per cell (QC) |
| `MAX_FEATURES` | `6000` | Maximum genes per cell (QC) |
| `MAX_PERCENT_MT` | `10` | Maximum mitochondrial % (QC) |
| `MAX_CELLS_PER_SAMPLE` | `2000` | Downsampling ceiling per sample |
| `UMAP_DIMS` | `1:17` | PCA dimensions used for UMAP/clustering |
| `CLUSTER_RESOLUTION` | `0.3` | Seurat clustering resolution |
| `CNV_THRESHOLD` | `15` | CNV score cutoff for malignancy classification |

---

## Running the Pipeline

Scripts are run sequentially from a session rooted at the project directory. Each script reads the RDS saved by the previous step.

### Core Pipeline

```r
# Recommended: open the project in RStudio (renv activates automatically)
# Then source each script in order:

source(here("scripts/01_preprocessing.R"))  # → seurat_list (in memory)
source(here("scripts/02_qc_filtering.R"))   # → seurat_list filtered
source(here("scripts/03_integration.R"))    # → seurat_combined.rds
source(here("scripts/04_annotation.R"))     # → seurat_obj_annotated.rds
source(here("scripts/05_infercnv.R"))       # → seurat_obj_infercnv.rds (runs in cnv/)
source(here("scripts/06_DEG_analysis.R"))   # → degs_per_cluster.rds, filtered_degs_per_cluster.rds
```

| Step | Script | Input | Output |
| :--- | :--- | :--- | :--- |
| 1 | `01_preprocessing.R` | Raw count files in `data/` | `seurat_list` in memory |
| 2 | `02_qc_filtering.R` | `seurat_list` | Filtered `seurat_list`; `qc_plots_combined.png` |
| 3 | `03_integration.R` | `seurat_list` | `seurat_combined.rds` |
| 4 | `04_annotation.R` | `seurat_combined.rds` | `seurat_obj_annotated.rds` |
| 5 | `05_infercnv.R` | `seurat_obj_annotated.rds` | `seurat_obj_infercnv.rds`; CNV plots in `cnv/` |
| 6 | `06_DEG_analysis.R` | `seurat_obj_infercnv.rds` | `degs_per_cluster.rds`; CSVs and volcano plots in `deg_csv/`, `deg_vocanol/` |

> **Note — InferCNV (step 5):** The script sets its working directory to `cnv/` and expects a `filtered_gene_order.pos` gene position file there. InferCNV is computationally intensive; expect long runtimes.

### Downstream Analyses

Run independently after step 6. Each loads its inputs from saved RDS/CSV files.

```r
source(here("scripts/07_pathway_analysis.R"))   # KEGG enrichment per cell type
source(here("scripts/08_STRING_analysis.R"))     # PPI network for cell type of interest
source(here("scripts/09_cellchat_analysis.R"))   # Cell-cell communication
source(here("scripts/10_trajectory_analysis.R")) # Pseudotime / trajectory
```

| Script | Input | Output |
| :--- | :--- | :--- |
| `07_pathway_analysis.R` | `degs_per_cluster.rds`, `filtered_degs_per_cluster.rds` | KEGG dotplots and CSVs in `kegg_dotplot/`, `kegg_csv/` |
| `08_STRING_analysis.R` | `deg_csv/no_1_0_degs_<cell_type>_Normal_vs_PCC.csv` | Hub gene tables and network plots |
| `09_cellchat_analysis.R` | `seurat_obj_infercnv.rds` | CellChat objects (`.RData`), interaction plots |
| `10_trajectory_analysis.R` | `seurat_obj_annotated.rds` / steroidogenic subset | Pseudotime plots, GO enrichment in `Trajectory/` |

---

## Output Structure

```
results/          # Figures and plots
output_tables/    # CSV tables
cnv/              # InferCNV outputs
deg_csv/          # DEG tables
deg_vocanol/      # Volcano plots
kegg_csv/         # KEGG enrichment tables
kegg_dotplot/     # KEGG dotplots
Trajectory/       # Pseudotime outputs
```
