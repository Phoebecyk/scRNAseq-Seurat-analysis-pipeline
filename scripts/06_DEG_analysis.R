# KEGG Pathway Analysis for DEGs (Normal vs. PCC) per simple_cell_type ####
# Prepare object and define cell types
seurat_combined <- PrepSCTFindMarkers(seurat_combined)
DefaultAssay(seurat_combined) <- "SCT"
Idents(seurat_combined) <- "simple_cell_type"
cell_types <- unique(seurat_combined$simple_cell_type)

# Create output directories
dir.create("deg_csv", showWarnings = FALSE)
dir.create("deg_vocanol", showWarnings = FALSE)

# Initialize list for DEGs
degs_per_cluster <- list()
filtered_degs_per_cluster <- list()


# Loop through each simple cell type
for (ct in cell_types) {
  cat("Processing:", ct, "\n")

  cluster_subset <- subset(seurat_combined, idents = ct)

  # Check both groups present
  if (!all(c("Normal", "PCC") %in% cluster_subset$condition)) {
    cat("  Skipping — one or both groups missing\n")
    next
  }
  cluster_subset <- PrepSCTFindMarkers(cluster_subset)
  DefaultAssay(cluster_subset) <- "SCT"
  Idents(cluster_subset) <- "condition"

  # DefaultAssay(cluster_subset) <- "RNA"
  # cluster_subset <- SCTransform(cluster_subset, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  # Idents(cluster_subset) <- "condition"

  # DEG Analysis
  de_res <- FindMarkers(cluster_subset,
                        ident.1 = "Normal",
                        ident.2 = "PCC",
                        logfc.threshold = 0.25,
                        min.pct = 0.1,
                        test.use = "wilcox"
                        ,recorrect_umi=FALSE
  )

  sig_degs <- de_res %>%
    rownames_to_column("gene") %>%
    filter(abs(avg_log2FC) > 1, p_val_adj < 0.05)

  # Skip saving if no sig DEGs
  if (nrow(sig_degs) == 0) {
    cat("  No significant DEGs found\n")
    next
  }

  # Save significant DEGs
  write.csv(sig_degs, file = paste0("deg_csv/degs_", ct, "_Normal_vs_PCC_filtered.csv"), row.names = FALSE)

  # Store in list
  degs_per_cluster[[ct]] <- sig_degs

  # ---- Filtering by expression & % expressing ---- #

  # Average expression
  avg_expr <- AggregateExpression(cluster_subset, group.by = "condition", assays = "SCT", slot = "data")
  avg_df <- as.data.frame(avg_expr$SCT)
  avg_df$gene <- rownames(avg_df)

  # % of cells expressing
  pct_exp <- sapply(c("Normal", "PCC"), function(grp) {
    cells <- WhichCells(cluster_subset, expression = condition == grp)
    rowMeans(GetAssayData(cluster_subset, assay = "SCT", slot = "data")[, cells] > 0)
  })
  pct_df <- as.data.frame(pct_exp)
  colnames(pct_df) <- paste0(colnames(pct_df), "_pct")
  pct_df$gene <- rownames(pct_df)

  # Merge & filter
  filtered_degs <- sig_degs %>%
    left_join(avg_df, by = "gene") %>%
    left_join(pct_df, by = "gene") %>%
    filter(
      (`Normal` > 1.5 | `PCC` > 1.5),
      (Normal_pct > 0.2 | PCC_pct > 0.2)
    )


  # Save filtered DEGs if any
  if (nrow(filtered_degs) > 0) {
    write.csv(filtered_degs, file = paste0("deg_csv/filtered_degs_", ct, "_Normal_vs_PCC.csv"), row.names = FALSE)
  } else {
    cat("No DEGs passed advanced filtering for", ct, "\n")
  }

  #Store in list
  filtered_degs_per_cluster[[ct]] <- filtered_degs

  # ---- Additional filter: remove genes with pct.1 or pct.2 equal to 0 or 1 ---- #
  filtered_degs <- filtered_degs %>%
    filter(!(pct.1 %in% c(0, 1) | pct.2 %in% c(0, 1)))

  # Save filtered DEGs if any
  if (nrow(filtered_degs) > 0) {
    write.csv(filtered_degs,
              file = paste0("deg_csv/no_1_0_degs_", ct, "_Normal_vs_PCC.csv"),
              row.names = FALSE)
  } else {
    cat("No DEGs passed advanced filtering (expression + %exp + pct.1/pct.2 cutoff) for", ct, "\n")
  }

  # Store in list
  filtered_degs_per_cluster[[ct]] <- filtered_degs

}

#vocanol plot for filtered Degs
for (ct in names(filtered_degs_per_cluster)) {
  volcano_data <- filtered_degs_per_cluster[[ct]] %>%
    mutate(
      neg_log10_padj = -log10(p_val_adj),
      label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 2, gene, NA),
      significance = case_when(
        p_val_adj < 0.05 & abs(avg_log2FC) > 2 ~ "Strong DEG",
        p_val_adj < 0.05 ~ "Moderate DEG",
        TRUE ~ "Not Significant"
      )
    )

  if (nrow(volcano_data) == 0) {
    cat("Skipping volcano plot for", ct, "- no filtered DEGs.\n")
    next
  }

  x_range <- range(volcano_data$avg_log2FC, na.rm = TRUE)
  y_range <- range(volcano_data$neg_log10_padj, na.rm = TRUE)

  p <- ggplot(volcano_data, aes(x = avg_log2FC, y = neg_log10_padj)) +
    geom_point(aes(color = significance), size = 2) +
    geom_text(aes(label = label), size = 3, vjust = 1.5, check_overlap = TRUE) +
    scale_color_manual(values = c("Strong DEG" = "red", "Moderate DEG" = "blue", "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    xlim(x_range[1] - 1, x_range[2] + 1) +
    ylim(0, y_range[2] + 2) +
    theme_minimal() +
    labs(
      title = paste0("Volcano: ", ct, " (Normal vs PCC)"),
      x = "log2 Fold Change", y = "-log10 Adjusted P-value", color = "Significance"
    )

  ggsave(
    filename = paste0("deg_vocanol/volcano_filtered_", ct, "_Normal_vs_PCC.png"),
    plot = p, width = 8, height = 6, dpi = 300
  )
}

# Save DEG lists for use by downstream scripts (07_pathway_analysis.R, 08_STRING_analysis.R)
saveRDS(degs_per_cluster, file = "degs_per_cluster.rds")
saveRDS(filtered_degs_per_cluster, file = "filtered_degs_per_cluster.rds")
