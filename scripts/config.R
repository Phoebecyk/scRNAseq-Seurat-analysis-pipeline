# config.R — shared analysis parameters
# Source at the top of each script: source(here("scripts/config.R"))
#
# HOW TO USE:
#   1. Set DATASET below to the analysis you want to run: "csACT", "PCC", or "insulinoma"
#   2. Set IS_SNRNASEQ = TRUE when processing Parse Bioscience FFPE data
#   3. After each new integration run, update the cell_type_map for that dataset
#      (cluster numbers change between runs; cell type names should stay consistent)

library(here)

# ── Active dataset ─────────────────────────────────────────────────────────────
DATASET     <- "PCC"       # "csACT" | "PCC" | "insulinoma"
IS_SNRNASEQ <- FALSE       # TRUE for Parse Bioscience FFPE runs

# ── Organism ──────────────────────────────────────────────────────────────────
ORGANISM_KEGG <- "cfa"     # KEGG organism code (Canis familiaris)
ORG_DB        <- "org.Cf.eg.db"

# ── Reproducibility ───────────────────────────────────────────────────────────
SEED <- 123

# ── QC thresholds ─────────────────────────────────────────────────────────────
# scRNAseq (10X fresh/frozen)
SCRNASEQ_MIN_FEATURES   <- 500
SCRNASEQ_MAX_FEATURES   <- 6000
SCRNASEQ_MAX_PERCENT_MT <- 10

# snRNAseq (Parse Bioscience FFPE nuclei) — nuclei have near-zero MT RNA
SNRNASEQ_MIN_FEATURES   <- 200
SNRNASEQ_MAX_FEATURES   <- 6000
SNRNASEQ_MAX_PERCENT_MT <- 2

MIN_FEATURES   <- if (IS_SNRNASEQ) SNRNASEQ_MIN_FEATURES   else SCRNASEQ_MIN_FEATURES
MAX_FEATURES   <- if (IS_SNRNASEQ) SNRNASEQ_MAX_FEATURES   else SCRNASEQ_MAX_FEATURES
MAX_PERCENT_MT <- if (IS_SNRNASEQ) SNRNASEQ_MAX_PERCENT_MT else
                    SCRNASEQ_MAX_PERCENT_MT

# ── Downsampling ──────────────────────────────────────────────────────────────
MAX_CELLS_PER_SAMPLE <- 2000
MIN_CELLS_PER_SAMPLE <- 100

# ── Integration ───────────────────────────────────────────────────────────────
N_VARIABLE_FEATURES <- 3000
N_PCS               <- 20
UMAP_DIMS           <- 1:17
CLUSTER_RESOLUTION  <- 0.3

# ── DEG analysis ─────────────────────────────────────────────────────────────
DEG_LOGFC_THRESHOLD <- 0.25   # pre-filter passed to FindMarkers
DEG_MIN_PCT         <- 0.1
DEG_PADJ_CUTOFF     <- 0.05
DEG_LOG2FC_CUTOFF   <- 1      # minimum |log2FC| to call a DEG significant
DEG_STRONG_LOG2FC   <- 2      # threshold for "Strong DEG" label in volcano plots
DEG_AVG_EXPR_MIN    <- 1.5    # minimum average expression in either group
DEG_PCT_MIN         <- 0.2    # minimum fraction of cells expressing the gene

# ── Pathway analysis ──────────────────────────────────────────────────────────
KEGG_PVALUE_CUTOFF <- 0.05
KEGG_QVALUE_CUTOFF <- 0.05
TOP_N_PATHWAYS     <- 30

# ── InferCNV ──────────────────────────────────────────────────────────────────
INFERCNV_CUTOFF <- 0.1
CNV_THRESHOLD   <- 15

# ── Cell type vocabulary ───────────────────────────────────────────────────────
# Shared names used across all datasets where the same population appears.
# Tumour-type-specific names (e.g. Chromaffin, INS_high) are defined per dataset below.
SHARED_STROMAL_TYPES <- c("Endothelial", "Fibroblasts", "Pericytes", "SmoothMuscle")
SHARED_IMMUNE_TYPES  <- c("Macrophages", "T_cells", "NK_cells", "Neutrophils",
                           "Mast_cells", "pDC", "Dendritic_cells")
PROLIFERATING_TYPE   <- "Proliferating"

# ── Dataset-specific configurations ───────────────────────────────────────────
#
# cell_type_map: named character vector mapping cluster number (as character) to
#   cell type label. UPDATE THE KEYS after each new clustering run — the values
#   (cell type names) should remain stable so downstream scripts keep working.
#
# tumour_cells: the subset of cell types that are the tumour-of-interest.
#   DEG and PPI analyses are run across all cell types but results for these
#   cell types are highlighted and used for hub gene identification.
#
# cond1 / cond2: condition labels as they appear in the Seurat metadata column
#   "condition". cond1 is the reference (ident.1 in FindMarkers), so a negative
#   log2FC means upregulated in cond2 (i.e. upregulated in tumour).

DATASET_CONFIGS <- list(

  csACT = list(
    cond1 = "Normal",
    cond2 = "csACT",

    # condition assignment per sample name (used in script 01)
    sample_conditions = c(
      "NAD2016"        = "Normal",
      "NAD2017"        = "Normal",
      "csACTDeNooijer" = "csACT",
      "csACTAngel"     = "csACT",
      "csACTDiesel"    = "csACT"
    ),

    # cluster → cell type; update cluster keys after re-clustering
    cell_type_map = c(
      "1" = "Steroidogenic",
      "2" = "Dedifferentiated_Steroidogenic",
      "3" = "Endothelial",
      "4" = "Fibroblasts",
      "5" = "Macrophages",
      "6" = "Proliferating",
      "7" = "T_cells",
      "8" = "Pericytes"
    ),

    tumour_cells = c("Steroidogenic", "Dedifferentiated_Steroidogenic"),

    infercnv_ref_groups = c("NAD2016", "NAD2017"),

    # marker genes for annotation dot plot (used in script 04 only)
    marker_panel = list(
      "Steroidogenic"               = c("AKR1B1", "CYP11A1", "HSD3B2", "CYP21A2", "SF1"),
      "Dedifferentiated_Steroidogenic" = c("RIMS2", "SNAP25", "KCNIP4", "FGFR2",
                                           "RBFOX1", "MAP2", "SOX10", "PHOX2B", "NGFR"),
      "Endothelial"                 = c("PECAM1", "FLT1", "KDR", "EMCN", "EGFL7"),
      "Fibroblasts"                 = c("COL1A2", "COL3A1", "COL1A1", "COL6A1", "COL6A3"),
      "SmoothMuscle_Pericytes"      = c("ACTA2", "TAGLN", "PPP1R14A", "PLN",
                                        "RGS5", "PDGFRB", "CSPG4"),
      "Macrophages"                 = c("CD86", "PTPRC", "C1QC", "MSR1", "IL18"),
      "T_cells"                     = c("CD3E", "CD69", "ITK", "CD226"),
      "Chromaffin_ref"              = c("ENSCAFG00000024864", "CHGB", "DBH", "SYP",
                                        "PNMT", "TH")
    )
  ),

  PCC = list(
    cond1 = "Normal",
    cond2 = "PCC",

    sample_conditions = c(
      "f106"       = "Normal",
      "f107"       = "Normal",
      "PCCVanDijk" = "PCC",
      "PCCRose"    = "PCC"
      # PCCDakota excluded: <100 cells after QC
    ),

    cell_type_map = c(
      "1"  = "Neurosecretory_Chromaffin",
      "2"  = "Endothelial",
      "3"  = "Proliferating",
      "4"  = "Steroidogenic",
      "5"  = "Chromaffin",
      "6"  = "Fibroblasts",
      "7"  = "Pericytes",
      "8"  = "Macrophages",
      "9"  = "T_cells",
      "10" = "Fetal_progenitor",
      "11" = "Fetal_sympathoblasts"
    ),

    tumour_cells = c("Chromaffin", "Neurosecretory_Chromaffin"),

    infercnv_ref_groups = c("f106", "f107"),

    marker_panel = list(
      "Neurosecretory_Chromaffin" = c("SYT1", "EPHA5", "DAB1", "RORB"),
      "Chromaffin"                = c("TH", "DBH", "ENSCAFG00000024864", "PNMT", "SCG2"),
      "Fetal_progenitor"          = c("PHOX2B", "ASCL1", "ISL1", "SOX10"),
      "Fetal_sympathoblasts"      = c("MKI67", "TOP2A", "CNTN3", "CDH4"),
      "Steroidogenic"             = c("STAR", "CYP11A1", "CYP17A1", "MC2R"),
      "Endothelial"               = c("PECAM1", "FLT1", "VWF", "EMCN"),
      "Fibroblasts"               = c("COL1A2", "DCN", "MFAP5", "SERPINF1"),
      "SmoothMuscle_Pericytes"    = c("ACTA2", "TAGLN", "TPM2", "PRKG1",
                                      "RGS5", "PDGFRB", "BMPER"),
      "Macrophages"               = c("PTPRC", "CD86", "IKZF1", "C1QA"),
      "T_cells"                   = c("CD3E", "ITK", "SKAP1", "BCL11B")
    )
  ),

  insulinoma = list(
    cond1 = "Primary",
    cond2 = "Metastasis",

    # update with actual sample IDs from your data
    sample_conditions = c(
      "INS_P1"  = "Primary",
      "INS_P2"  = "Primary",
      "INS_P3"  = "Primary",
      "INS_LM1" = "Metastasis",
      "INS_LM2" = "Metastasis",
      "INS_LN1" = "Metastasis"
    ),

    cell_type_map = c(
      "1" = "INS_high",
      "2" = "INS_low",
      "3" = "T_cells",
      "4" = "Macrophages",
      "5" = "Endothelial",
      "6" = "Fibroblasts",
      "7" = "Neutrophils",
      "8" = "pDC",
      "9" = "Mast_cells"
    ),

    tumour_cells = c("INS_high", "INS_low"),

    infercnv_ref_groups = NULL,  # no matched normal for insulinoma

    marker_panel = list(
      "INS_high"    = c("INS", "ENSCAFG00000024864", "CHGB", "PDX1"),
      "INS_low"     = c("GPR39", "PRLR", "RIMS2", "ENOX1", "NPAS3", "GLIS1"),
      "T_cells"     = c("CD3E", "CD3D", "CD2", "CCL5"),
      "Macrophages" = c("TYROBP", "CD86", "AIF1", "C1QA", "C1QB"),
      "Endothelial" = c("PODXL", "PECAM1", "FLT1"),
      "Fibroblasts" = c("FAP", "COL1A2", "COL3A1"),
      "Neutrophils" = c("S100A8", "S100A12", "CSF3R"),
      "pDC"         = c("FLT3", "BCL11A", "PLD4"),
      "Mast_cells"  = c("FCER1A", "BLK", "MYB")
    )
  )
)

# ── Shared display aliases ─────────────────────────────────────────────────────
# Ensembl IDs without a human-readable symbol — used to relabel dot plot axes.
# Add entries here if new unannotated genes appear in marker panels.
ENSEMBL_ALIASES <- c("ENSCAFG00000024864" = "CHGA")

# ── Active configuration (referenced throughout all scripts) ───────────────────
CFG           <- DATASET_CONFIGS[[DATASET]]
COND1         <- CFG$cond1           # reference condition
COND2         <- CFG$cond2           # tumour condition
TUMOUR_CELLS  <- CFG$tumour_cells    # cell types for DEG/PPI focus
CELL_TYPE_MAP <- CFG$cell_type_map   # cluster → cell type name

# Output paths — all results land in dataset-specific subdirectories
RESULTS_DIR       <- here("results", DATASET)
OUTPUT_TABLES_DIR <- here("output_tables", DATASET)
