source(here::here("scripts/config.R"))

# DEG analysis: COND1 vs COND2, per cell type
# COND1 is the reference (ident.1), so:
#   negative log2FC = upregulated in COND2 (tumour)
#   positive log2FC = downregulated in COND2 (tumour)
#
# Full analysis runs across all cell types.
# Results for TUMOUR_CELLS are additionally saved with a "_tumour" suffix
# for direct use in pathway and PPI analyses (scripts 07, 08).

seurat_combined <- readRDS(file.path(RESULTS_DIR, "seurat_obj_annotated.rds"))
seurat_combined <- PrepSCTFindMarkers(seurat_combined)
DefaultAssay(seurat_combined) <- "SCT"
Idents(seurat_combined) <- "simple_cell_type"

cell_types <- unique(seurat_combined$simple_cell_type)

dir.create(file.path(OUTPUT_TABLES_DIR, "deg_csv"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS_DIR,       "deg_volcano"), showWarnings = FALSE, recursive = TRUE)

degs_per_cluster          <- list()
filtered_degs_per_cluster <- list()

for (ct in cell_types) {
  cat("Processing:", ct, "\n")

  cluster_subset <- subset(seurat_combined, idents = ct)

  if (!all(c(COND1, COND2) %in% cluster_subset$condition)) {
    cat("  Skipping —", COND1, "or", COND2, "not present\n")
    next
  }
  cluster_subset <- PrepSCTFindMarkers(cluster_subset)
  DefaultAssay(cluster_subset) <- "SCT"
  Idents(cluster_subset) <- "condition"

  de_res <- FindMarkers(
    cluster_subset,
    ident.1         = COND1,
    ident.2         = COND2,
    logfc.threshold = DEG_LOGFC_THRESHOLD,
    min.pct         = DEG_MIN_PCT,
    test.use        = "wilcox",
    recorrect_umi   = FALSE
  )

  sig_degs <- de_res %>%
    rownames_to_column("gene") %>%
    filter(abs(avg_log2FC) > DEG_LOG2FC_CUTOFF, p_val_adj < DEG_PADJ_CUTOFF)

  if (nrow(sig_degs) == 0) {
    cat("  No significant DEGs found\n")
    next
  }

  out_prefix <- file.path(
    OUTPUT_TABLES_DIR, "deg_csv",
    paste0("degs_", ct, "_", COND1, "_vs_", COND2)
  )
  write.csv(sig_degs, file = paste0(out_prefix, ".csv"), row.names = FALSE)
  degs_per_cluster[[ct]] <- sig_degs

  # ── Expression-level and cell-fraction filters ─────────────────────────────
  avg_expr <- AggregateExpression(
    cluster_subset, group.by = "condition", assays = "SCT", slot = "data"
  )
  avg_df        <- as.data.frame(avg_expr$SCT)
  avg_df$gene   <- rownames(avg_df)

  pct_exp <- sapply(c(COND1, COND2), function(grp) {
    cells <- WhichCells(cluster_subset, expression = condition == grp)
    rowMeans(GetAssayData(cluster_subset, assay = "SCT", slot = "data")[, cells] > 0)
  })
  pct_df         <- as.data.frame(pct_exp)
  colnames(pct_df) <- paste0(c(COND1, COND2), "_pct")
  pct_df$gene    <- rownames(pct_df)

  cond1_pct_col <- paste0(COND1, "_pct")
  cond2_pct_col <- paste0(COND2, "_pct")

  filtered_degs <- sig_degs %>%
    left_join(avg_df, by = "gene") %>%
    left_join(pct_df, by = "gene") %>%
    filter(
      (.data[[COND1]] > DEG_AVG_EXPR_MIN | .data[[COND2]] > DEG_AVG_EXPR_MIN),
      (.data[[cond1_pct_col]] > DEG_PCT_MIN | .data[[cond2_pct_col]] > DEG_PCT_MIN)
    ) %>%
    filter(!(pct.1 %in% c(0, 1) | pct.2 %in% c(0, 1)))

  if (nrow(filtered_degs) > 0) {
    write.csv(
      filtered_degs,
      file = paste0(out_prefix, "_filtered.csv"),
      row.names = FALSE
    )
  } else {
    cat("  No DEGs passed expression + cell-fraction filters\n")
  }

  filtered_degs_per_cluster[[ct]] <- filtered_degs

  # ── Volcano plot ───────────────────────────────────────────────────────────
  if (nrow(filtered_degs) == 0) next

  volcano_data <- filtered_degs %>%
    mutate(
      neg_log10_padj = -log10(p_val_adj),
      label          = ifelse(
        p_val_adj < DEG_PADJ_CUTOFF & abs(avg_log2FC) > DEG_STRONG_LOG2FC,
        gene, NA
      ),
      significance = case_when(
        p_val_adj < DEG_PADJ_CUTOFF & abs(avg_log2FC) > DEG_STRONG_LOG2FC ~ "Strong DEG",
        p_val_adj < DEG_PADJ_CUTOFF                                        ~ "Moderate DEG",
        TRUE                                                                ~ "Not Significant"
      )
    )

  p_volcano <- ggplot(volcano_data, aes(x = avg_log2FC, y = neg_log10_padj)) +
    geom_point(aes(color = significance), size = 2) +
    geom_text(aes(label = label), size = 3, vjust = 1.5, check_overlap = TRUE) +
    scale_color_manual(
      values = c("Strong DEG" = "red", "Moderate DEG" = "blue",
                 "Not Significant" = "gray")
    ) +
    geom_vline(xintercept = c(-DEG_STRONG_LOG2FC, DEG_STRONG_LOG2FC),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(DEG_PADJ_CUTOFF), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste0("Volcano: ", ct, " (", COND1, " vs ", COND2, ")"),
      x     = paste0("log2FC  [negative = up in ", COND2, "]"),
      y     = "-log10 Adjusted P-value",
      color = "Significance"
    )

  ggsave(
    file.path(RESULTS_DIR, "deg_volcano",
              paste0("volcano_", ct, "_", COND1, "_vs_", COND2, ".png")),
    p_volcano, width = 8, height = 6, dpi = 300
  )
}

# ── Save DEG lists for scripts 07 and 08 ─────────────────────────────────────
saveRDS(degs_per_cluster,
        file = file.path(RESULTS_DIR, "degs_per_cluster.rds"))
saveRDS(filtered_degs_per_cluster,
        file = file.path(RESULTS_DIR, "filtered_degs_per_cluster.rds"))

# ── Convenience: tumour-cell-only subsets for PPI analysis ───────────────────
# These contain only the cell types listed in TUMOUR_CELLS (config.R),
# making it easy for script 08 to focus on the relevant populations.
saveRDS(
  degs_per_cluster[intersect(names(degs_per_cluster), TUMOUR_CELLS)],
  file = file.path(RESULTS_DIR, "degs_tumour_cells.rds")
)
saveRDS(
  filtered_degs_per_cluster[intersect(names(filtered_degs_per_cluster), TUMOUR_CELLS)],
  file = file.path(RESULTS_DIR, "filtered_degs_tumour_cells.rds")
)
