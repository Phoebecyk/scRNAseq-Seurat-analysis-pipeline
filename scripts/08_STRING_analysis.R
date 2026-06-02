source(here::here("scripts/config.R"))
library(readr)
library(dplyr)
library(STRINGdb)
library(igraph)

# STRING PPI + hub gene identification
# Runs over every cell type in TUMOUR_CELLS (config.R).
# For each cell type, uses the filtered DEG list produced by script 06.
# Hub genes = outliers in BOTH betweenness centrality (BC) and closeness
# centrality (CC), identified via the IQR fence method.
# A composite hub score (mean of min-max normalised BC and CC) is also
# computed and the top 5% saved separately.

# ── Load filtered DEG lists ───────────────────────────────────────────────────
if (!exists("filtered_degs_tumour")) {
  filtered_degs_tumour <- readRDS(
    file.path(RESULTS_DIR, "filtered_degs_tumour_cells.rds")
  )
}

# ── STRING database connection ────────────────────────────────────────────────
# species 9615 = Canis lupus familiaris
STRING_SPECIES        <- 9615
STRING_SCORE_THRESHOLD <- 700   # high-confidence interactions only

string_db <- STRINGdb$new(
  version         = "12",
  species         = STRING_SPECIES,
  score_threshold = STRING_SCORE_THRESHOLD
)

# ── Output directories ────────────────────────────────────────────────────────
ppi_dir <- file.path(OUTPUT_TABLES_DIR, "ppi")
dir.create(ppi_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helper: IQR outlier threshold (upper fence) ───────────────────────────────
iqr_upper <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  q[2] + 1.5 * (q[2] - q[1])
}

# ── Main loop: one STRING analysis per tumour cell type ───────────────────────
for (ct in TUMOUR_CELLS) {

  if (is.null(filtered_degs_tumour[[ct]]) ||
        nrow(filtered_degs_tumour[[ct]]) == 0) {
    cat("Skipping", ct, "— no filtered DEGs available\n")
    next
  }

  cat("\nSTRING analysis:", ct, "\n")
  deg_table  <- filtered_degs_tumour[[ct]]
  degs_list  <- deg_table$gene

  # ── Map gene symbols to STRING IDs ─────────────────────────────────────────
  mapped <- string_db$map(
    data.frame(gene = degs_list), "gene", removeUnmappedRows = TRUE
  )
  if (nrow(mapped) == 0) {
    cat("  No genes mapped to STRING for", ct, "\n")
    next
  }

  # ── Build PPI graph and compute centrality ─────────────────────────────────
  ppi_graph <- string_db$get_graph()

  bc_scores <- betweenness(ppi_graph, directed = FALSE)
  cc_scores <- closeness(ppi_graph, mode = "all")

  centrality_df <- data.frame(
    STRING_id = names(bc_scores),
    bc_score  = bc_scores,
    cc_score  = cc_scores
  ) %>%
    left_join(mapped, by = "STRING_id") %>%
    dplyr::select(gene, bc_score, cc_score) %>%
    na.omit()

  if (nrow(centrality_df) == 0) {
    cat("  No centrality scores computed for", ct, "\n")
    next
  }

  # ── IQR-based hub genes (outlier in BOTH BC and CC) ───────────────────────
  bc_thresh <- iqr_upper(centrality_df$bc_score)
  cc_thresh <- iqr_upper(centrality_df$cc_score)

  hub_genes_iqr <- centrality_df %>%
    filter(bc_score > bc_thresh & cc_score > cc_thresh)

  # ── Composite hub score: top 5% by mean of normalised BC and CC ───────────
  min_max <- function(x) (x - min(x)) / (max(x) - min(x))

  hub_genes_composite <- centrality_df %>%
    mutate(
      bc_norm   = min_max(bc_score),
      cc_norm   = min_max(cc_score),
      hub_score = (bc_norm + cc_norm) / 2
    ) %>%
    arrange(desc(hub_score)) %>%
    slice_head(n = max(1L, ceiling(0.05 * nrow(centrality_df))))

  # ── Save results ───────────────────────────────────────────────────────────
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)   # safe filename

  write.csv(
    centrality_df,
    file.path(ppi_dir, paste0("centrality_all_", ct_safe, ".csv")),
    row.names = FALSE
  )
  write.csv(
    hub_genes_iqr,
    file.path(ppi_dir, paste0("hub_genes_iqr_", ct_safe, ".csv")),
    row.names = FALSE
  )
  write.csv(
    hub_genes_composite,
    file.path(ppi_dir, paste0("hub_genes_composite_", ct_safe, ".csv")),
    row.names = FALSE
  )

  cat("  IQR hub genes (", nrow(hub_genes_iqr), "):",
      paste(hub_genes_iqr$gene, collapse = ", "), "\n")
  cat("  Top-5% composite hub genes (", nrow(hub_genes_composite), "):",
      paste(hub_genes_composite$gene, collapse = ", "), "\n")
}

cat("\nSTRING analysis complete. Results in:", ppi_dir, "\n")
