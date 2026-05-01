# Load DEG lists from 06 if not already in memory
if (!exists("degs_per_cluster")) {
  degs_per_cluster <- readRDS("degs_per_cluster.rds")
}
if (!exists("filtered_degs_per_cluster")) {
  filtered_degs_per_cluster <- readRDS("filtered_degs_per_cluster.rds")
}

#KEGG
# Output directory
dir.create("kegg_csv", showWarnings = FALSE)
dir.create("kegg_dotplot", showWarnings = FALSE)

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
  gene_entrez <- bitr(degs_combined, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Cf.eg.db)
  if (nrow(gene_entrez) == 0) {
    cat("    No Entrez ID mapping\n")
    next
  }

  # 3. Run KEGG enrichment
  kegg_res <- enrichKEGG(gene = gene_entrez$ENTREZID,
                         organism = "cfa",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)

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

  # Save KEGG CSV (simplified filename)
  csv_file_name <- paste0("kegg_csv/kegg_RichFactor_", ct, ".csv")
  write.csv(kegg_df, file = csv_file_name, row.names = FALSE)

  # Select top 20 by RichFactor, arrange, and set factor levels for plotting order
  kegg_plot_df <- kegg_df %>%
    arrange(desc(RichFactor)) %>%
    slice_head(n = 30) %>%
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

  # Save the plot (simplified filename)
  plot_file_name <- paste0("kegg_dotplot/kegg_dotplot_RichFactor_", ct, ".png")
  ggsave(filename = plot_file_name, plot = rich_plot, width = 8, height = 10, dpi = 300)
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
  gene_entrez <- bitr(degs_combined, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Cf.eg.db)
  if (nrow(gene_entrez) == 0) {
    cat("    No Entrez ID mapping\n")
    next
  }

  # 3. Run KEGG enrichment
  kegg_res <- enrichKEGG(gene = gene_entrez$ENTREZID,
                         organism = "cfa",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)

  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
    cat("    No significant KEGG pathways enriched\n")
    next
  }

  # 4. Process results, calculate RichFactor
  kegg_df <- as.data.frame(kegg_res)
  kegg_df$RichFactor <- kegg_df$Count / as.numeric(sub("/.*", "", kegg_df$BgRatio))

  # Store results
  kegg_per_cluster_filtered[[ct]] <- kegg_df

  # Save KEGG CSV (with "_filtered" in filename)
  csv_file_name <- paste0("kegg_csv_filtered/kegg_RichFactor_", ct, "_filtered.csv")
  write.csv(kegg_df, file = csv_file_name, row.names = FALSE)

  # Select top 30 by RichFactor
  kegg_plot_df <- kegg_df %>%
    arrange(desc(RichFactor)) %>%
    slice_head(n = 30) %>%
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

  # Save the plot (with "_filtered" in filename)
  plot_file_name <- paste0("kegg_dotplot_filtered/kegg_dotplot_RichFactor_", ct, "_filtered.png")
  ggsave(filename = plot_file_name, plot = rich_plot, width = 8, height = 10, dpi = 300)
}
