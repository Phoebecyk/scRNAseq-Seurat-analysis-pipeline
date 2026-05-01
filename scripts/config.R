# config.R — shared analysis parameters
# Source at the top of each script: source(here("scripts/config.R"))

library(here)

# ── Organism ──────────────────────────────────────────────────────────────────
ORGANISM_KEGG <- "cfa"           # KEGG organism code (Canis familiaris)
ORG_DB        <- "org.Cf.eg.db"  # Bioconductor OrgDb package

# ── Reproducibility ───────────────────────────────────────────────────────────
SEED <- 123

# ── Paths (relative to project root) ─────────────────────────────────────────
DATA_DIR          <- here("data")
RESULTS_DIR       <- here("results")
OUTPUT_TABLES_DIR <- here("output_tables")

# ── QC thresholds ─────────────────────────────────────────────────────────────
MIN_FEATURES   <- 500
MAX_FEATURES   <- 6000
MAX_PERCENT_MT <- 10

# ── Downsampling ──────────────────────────────────────────────────────────────
MAX_CELLS_PER_SAMPLE <- 2000
MIN_CELLS_PER_SAMPLE <- 100

# ── Integration ───────────────────────────────────────────────────────────────
N_VARIABLE_FEATURES <- 3000
N_PCS               <- 20
UMAP_DIMS           <- 1:17
CLUSTER_RESOLUTION  <- 0.3

# ── DEG analysis ─────────────────────────────────────────────────────────────
DEG_LOGFC_THRESHOLD <- 0.25   # Pre-filter passed to FindMarkers
DEG_MIN_PCT         <- 0.1
DEG_PADJ_CUTOFF     <- 0.05
DEG_LOG2FC_CUTOFF   <- 1      # For selecting significant DEGs
DEG_STRONG_LOG2FC   <- 2      # For "Strong DEG" label in volcano plots
DEG_AVG_EXPR_MIN    <- 1.5    # Minimum average expression in either group
DEG_PCT_MIN         <- 0.2    # Minimum fraction of cells expressing gene

# ── Pathway analysis ──────────────────────────────────────────────────────────
KEGG_PVALUE_CUTOFF <- 0.05
KEGG_QVALUE_CUTOFF <- 0.05
TOP_N_PATHWAYS     <- 30

# ── InferCNV ──────────────────────────────────────────────────────────────────
INFERCNV_CUTOFF <- 0.1
CNV_THRESHOLD   <- 15   # Score above which a steroidogenic cell is called malignant
