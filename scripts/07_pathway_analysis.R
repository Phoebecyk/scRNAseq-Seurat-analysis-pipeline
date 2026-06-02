source(here::here("scripts/config.R"))

# Load DEG lists from 06 if not already in memory
if (!exists("degs_per_cluster")) {
  degs_per_cluster <- readRDS(file.path(RESULTS_DIR, "degs_per_cluster.rds"))
}
if (!exists("filtered_degs_per_cluster")) {
  filtered_degs_per_cluster <- readRDS(
    file.path(RESULTS_DIR, "filtered_degs_per_cluster.rds")
  )
}

kegg_dir          <- file.path(OUTPUT_TABLES_DIR, "kegg_csv")
kegg_plot_dir     <- file.path(RESULTS_DIR,       "kegg_dotplot")
kegg_filt_dir     <- file.path(OUTPUT_TABLES_DIR, "kegg_csv_filtered")
kegg_filt_plot_dir <- file.path(RESULTS_DIR,      "kegg_dotplot_filtered")

dir.create(kegg_dir,           showWarnings = FALSE, recursive = TRUE)
dir.create(kegg_plot_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(kegg_filt_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(kegg_filt_plot_dir, showWarnings = FALSE, recursive = TRUE)

kegg_per_cluster <- list()

# Loop through each cell type/cluster in your differential gene expression list
for (ct in names(degs_per_cluster)) {
  deg_table <- degs_per_cluster[[ct]]
  cat("KEGG combined analysis for:", ct, "\n")

  # 1. Collect all significant DEGs
  # Since deg_table (sig_degs) is already filtered for abs(avg_log2FC) > 1 and p_val_adj < 0.05,
  # we simply pull all genes from the table.
  degs_combined <- deg_table %>% pull(gene)

  if (length(degs_combined) == 0) {
    cat("    No significant DEGs found for analysis.\n")
    next
  }

  # 2. Map gene symbols to Entrez IDs (required for enrichKEGG)
  gene_entrez <- bitr(degs_combined,
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = ORG_DB)
  if (nrow(gene_entrez) == 0) {
    cat("    No Entrez ID mapping\n")
    next
  }

  # 3. Run KEGG enrichment
  kegg_res <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                         organism     = ORGANISM_KEGG,
                         pvalueCutoff = KEGG_PVALUE_CUTOFF,
                         qvalueCutoff = KEGG_QVALUE_CUTOFF)

  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
    cat("    No significant KEGG pathways enriched\n")
    next
  }

  # 4. Process results, calculate RichFactor, and prepare for plotting
  kegg_df <- as.data.frame(kegg_res)

  # Calculate RichFactor = DEGs in pathway / all genes in pathway
  kegg_df$RichFactor <- kegg_df$Count / as.numeric(sub("/.*", "", kegg_df$BgRatio))

  # Store using the cluster name only
  kegg_per_cluster[[ct]] <- kegg_df

  write.csv(kegg_df,
            file = file.path(kegg_dir, paste0("kegg_RichFactor_", ct, ".csv")),
            row.names = FALSE)

  kegg_plot_df <- kegg_df %>%
    arrange(desc(RichFactor)) %>%
    slice_head(n = TOP_N_PATHWAYS) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  # 5. Create custom RichFactor DotPlot 📊
  rich_plot <- ggplot(kegg_plot_df, aes(x = RichFactor, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "salmon", high = "lightblue", name = "Adjusted p-value") +
    scale_size_continuous(name = "DEG count") +
    labs(
      title = paste("KEGG Pathways (RichFactor):", ct),
      x = "RichFactor (DEGs in pathway / all genes in pathway)",
      y = "Pathway"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(5, 5, 5, 5), "mm"))

  print(rich_plot)

  ggsave(
    filename = file.path(kegg_plot_dir, paste0("kegg_RichFactor_", ct, ".png")),
    plot = rich_plot, width = 8, height = 10, dpi = 300
  )
}

# KEGG (using filtered DEGs)
# Output directory (separate subfolders for clarity)
dir.create("kegg_csv_filtered", showWarnings = FALSE)
dir.create("kegg_dotplot_filtered", showWarnings = FALSE)

kegg_per_cluster_filtered <- list()

# Loop through each cell type/cluster in filtered DEGs
for (ct in names(filtered_degs_per_cluster)) {
  deg_table <- filtered_degs_per_cluster[[ct]]
  cat("KEGG filtered analysis for:", ct, "\n")

  # 1. Collect all significant DEGs (already filtered)
  degs_combined <- deg_table %>% pull(gene)

  if (length(degs_combined) == 0) {
    cat("    No filtered DEGs found for analysis.\n")
    next
  }

  # 2. Map gene symbols to Entrez IDs
  gene_entrez <- bitr(degs_combined,
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = ORG_DB)
  if (nrow(gene_entrez) == 0) {
    cat("    No Entrez ID mapping\n")
    next
  }

  # 3. Run KEGG enrichment
  kegg_res <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                         organism     = ORGANISM_KEGG,
                         pvalueCutoff = KEGG_PVALUE_CUTOFF,
                         qvalueCutoff = KEGG_QVALUE_CUTOFF)

  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
    cat("    No significant KEGG pathways enriched\n")
    next
  }

  # 4. Process results, calculate RichFactor
  kegg_df <- as.data.frame(kegg_res)
  kegg_df$RichFactor <- kegg_df$Count / as.numeric(sub("/.*", "", kegg_df$BgRatio))

  # Store results
  kegg_per_cluster_filtered[[ct]] <- kegg_df

  write.csv(kegg_df,
            file = file.path(kegg_filt_dir,
                             paste0("kegg_RichFactor_", ct, "_filtered.csv")),
            row.names = FALSE)

  kegg_plot_df <- kegg_df %>%
    arrange(desc(RichFactor)) %>%
    slice_head(n = TOP_N_PATHWAYS) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  # 5. Custom RichFactor DotPlot 📊
  rich_plot <- ggplot(kegg_plot_df, aes(x = RichFactor, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "salmon", high = "lightblue", name = "Adjusted p-value") +
    scale_size_continuous(name = "DEG count") +
    labs(
      title = paste("KEGG Pathways (RichFactor, filtered):", ct),
      x = "RichFactor (DEGs in pathway / all genes in pathway)",
      y = "Pathway"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(5, 5, 5, 5), "mm"))

  print(rich_plot)

  ggsave(
    filename = file.path(kegg_filt_plot_dir,
                         paste0("kegg_RichFactor_", ct, "_filtered.png")),
    plot = rich_plot, width = 8, height = 10, dpi = 300
  )
}
