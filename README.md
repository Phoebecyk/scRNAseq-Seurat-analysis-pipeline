# scRNAseq Seurat Analysis Pipeline

A reproducible single-cell RNA-seq workflow for three canine neuroendocrine tumour types — cortisol-secreting adrenal cortical tumour (csACT), pheochromocytoma (PCC), and insulinoma — using R and Seurat v5. Human fetal adrenal samples (GEO: GSE137804) are used as a normal medullary reference for PCC in the absence of sufficient canine adrenal medullary cells.

**Translational design:** Canine samples are aligned to the CanFam3.1 reference genome, then at the end of preprocessing all canine gene names are converted to their human ortholog symbols (Ensembl BioMart, 1:1 high-confidence orthologs only). From that point onward the entire pipeline operates in human gene space — enabling direct use of human reference databases for cell type annotation, KEGG pathway enrichment (`org.Hs.eg.db`), STRING PPI (human network), and CellChat (`CellChatDB.human`). The ortholog table is cached to `results/dog_human_orthologs.csv` after the first run.

Switch between tumour types by changing one line in `scripts/config.R`; all scripts pick up the correct conditions, cell type labels, and output paths automatically.

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

Place raw data under `adrenal_scRNA_seq_data/` (or `insulinoma_scRNA_seq_data/`) before running the pipeline:

```
adrenal_scRNA_seq_data/
├── GSE137804_RAW/
│   ├── GSM4088787_F106_gene_cell_exprs_table.xls.gz   # Human fetal adrenal (PCC normal ref)
│   └── GSM4088788_F107_gene_cell_exprs_table.xls.gz
├── PCCVanDijk/filtered_feature_bc_matrix/             # Canine PCC
├── PCCDakota/filtered_feature_bc_matrix/
├── PCCRose/filtered_feature_bc_matrix/
├── NAD2016/filtered_feature_bc_matrix/                # Canine normal adrenal (csACT ref)
├── NAD2017/filtered_feature_bc_matrix/
├── csACTDeNooijer/filtered_feature_bc_matrix/         # Canine csACT
├── csACTAngel/filtered_feature_bc_matrix/
└── csACTDiesel/filtered_feature_bc_matrix/

insulinoma_scRNA_seq_data/
├── INS_P1/filtered_feature_bc_matrix/                 # Canine insulinoma — primary
├── INS_P2/filtered_feature_bc_matrix/
├── INS_P3/filtered_feature_bc_matrix/
├── INS_LM1/filtered_feature_bc_matrix/                # Liver metastasis
├── INS_LM2/filtered_feature_bc_matrix/
└── INS_LN1/filtered_feature_bc_matrix/                # Lymph node metastasis
```

---

## Configuration

All shared parameters live in `scripts/config.R`. It is sourced automatically at the top of every script.

**To switch datasets**, change the two lines at the top of `config.R`:

```r
DATASET     <- "PCC"   # "csACT" | "PCC" | "insulinoma"
IS_SNRNASEQ <- FALSE   # TRUE for Parse Bioscience FFPE nuclei data
```

Everything else — condition labels, cell type maps, marker panels, output paths — is derived from those two settings.

Key parameters:

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `DATASET` | `"PCC"` | Active tumour type; controls all dataset-specific settings |
| `IS_SNRNASEQ` | `FALSE` | If `TRUE`, uses snRNAseq QC thresholds (min features 200, MT% < 2) |
| `COND1` / `COND2` | dataset-specific | Reference / tumour condition labels (e.g. `"Normal"` / `"PCC"`) |
| `TUMOUR_CELLS` | dataset-specific | Cell types used for DEG and PPI focus analyses |
| `CELL_TYPE_MAP` | dataset-specific | Named vector mapping cluster number → cell type label |
| `ORGANISM_KEGG` | `"hsa"` | KEGG organism code (Homo sapiens) — genes are in human symbol space after ortholog conversion |
| `MIN_FEATURES` | `500` / `200` | Minimum genes per cell (scRNAseq / snRNAseq) |
| `MAX_PERCENT_MT` | `10` / `2` | Maximum mitochondrial % (scRNAseq / snRNAseq) |
| `MAX_CELLS_PER_SAMPLE` | `2000` | Downsampling ceiling per sample |
| `UMAP_DIMS` | `1:17` | PCA dimensions used for UMAP/clustering |
| `CLUSTER_RESOLUTION` | `0.3` | Seurat clustering resolution |
| `CNV_THRESHOLD` | `15` | CNV score cutoff for malignancy classification |

**After re-clustering**, update the `cell_type_map` keys for the active dataset in `DATASET_CONFIGS`. The `marker_panel` (dot plot marker genes) is also defined per dataset there and rarely needs to change.

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

| Step | Script | Input | Output (all under `results/DATASET/` or `output_tables/DATASET/`) |
| :--- | :--- | :--- | :--- |
| 1 | `01_preprocessing.R` | Raw count files | `seurat_list` in memory; `results/dog_human_orthologs.csv` (cached on first run) |
| 2 | `02_qc_filtering.R` | `seurat_list` | Filtered `seurat_list`; `qc_plots_combined.png` |
| 3 | `03_integration.R` | `seurat_list` | `seurat_combined.rds` |
| 4 | `04_annotation.R` | `seurat_combined.rds` | `seurat_obj_annotated.rds`; UMAP and dot plots |
| 5 | `05_infercnv.R` | `seurat_obj_annotated.rds` | `seurat_obj_infercnv.rds`; CNV plots in `results/DATASET/cnv/` |
| 6 | `06_DEG_analysis.R` | `seurat_obj_infercnv.rds` | `degs_per_cluster.rds`, `degs_tumour_cells.rds`; volcano plots |

> **Note — InferCNV (step 5):** Expects a `filtered_gene_order.pos` gene position file in the `cnv/` subdirectory. Computationally intensive; expect long runtimes.

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
| `07_pathway_analysis.R` | `degs_per_cluster.rds`, `filtered_degs_per_cluster.rds` | KEGG dotplots and CSVs in `results/DATASET/kegg_dotplot/`, `output_tables/DATASET/kegg_csv/` |
| `08_STRING_analysis.R` | `filtered_degs_tumour_cells.rds` | Hub gene tables and network plots per tumour cell type in `output_tables/DATASET/ppi/` |
| `09_cellchat_analysis.R` | `seurat_obj_infercnv.rds` | CellChat objects (`.RData`), interaction plots |
| `10_trajectory_analysis.R` | `seurat_obj_annotated.rds` / tumour cell subset | Pseudotime plots, GO enrichment in `Trajectory/` |

---

## Output Structure

All outputs land in dataset-specific subdirectories so results from different tumour types never overwrite each other.

```
results/
├── dog_human_orthologs.csv       # BioMart 1:1 ortholog table (cached after first run)
└── {DATASET}/                    # e.g. results/PCC/
    ├── umap_annotated.png
    ├── umap_annotated_split_condition.png
    ├── dotplot_cluster_markers.png
    ├── top_cluster_marker_heatmap.png
    ├── cnv/                      # InferCNV outputs
    ├── kegg_dotplot/
    └── kegg_dotplot_filtered/

output_tables/
└── {DATASET}/                    # e.g. output_tables/PCC/
    ├── all_markers_per_cluster.csv
    ├── top10_markers_per_cluster.csv
    ├── degs_per_cluster.rds
    ├── filtered_degs_per_cluster.rds
    ├── degs_tumour_cells.rds
    ├── filtered_degs_tumour_cells.rds
    ├── kegg_csv/
    ├── kegg_csv_filtered/
    └── ppi/                      # STRING hub gene outputs per tumour cell type
```
