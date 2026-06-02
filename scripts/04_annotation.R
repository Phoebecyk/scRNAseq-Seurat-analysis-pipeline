source(here::here("scripts/config.R"))

# ── Cluster marker genes ───────────────────────────────────────────────────────
DefaultAssay(seurat_combined) <- "SCT"
seurat_combined <- PrepSCTFindMarkers(seurat_combined)

all_markers <- FindAllMarkers(
  seurat_combined,
  only.pos          = TRUE,
  min.pct           = DEG_MIN_PCT,
  logfc.threshold   = DEG_LOGFC_THRESHOLD
)
write.csv(all_markers, file = file.path(OUTPUT_TABLES_DIR, "all_markers_per_cluster.csv"))

# Top 10 markers per cluster for heatmap
top10 <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
write.csv(top10, file = file.path(OUTPUT_TABLES_DIR, "top10_markers_per_cluster.csv"))

seurat_combined <- ScaleData(seurat_combined, features = top10$gene, assay = "SCT")
p_heatmap <- DoHeatmap(seurat_combined, features = top10$gene) + NoLegend()
ggsave(
  file.path(RESULTS_DIR, "top_cluster_marker_heatmap.png"),
  p_heatmap, width = 10, height = 15, dpi = 300
)

# ── Cell type marker dot plot ──────────────────────────────────────────────────
# marker_panel and ENSEMBL_ALIASES come from config.R — edit them there.
marker_panel <- CFG$marker_panel

p_dot <- DotPlot(
  seurat_combined,
  features   = marker_panel,
  group.by   = "seurat_clusters",
  dot.scale  = 6
) +
  scale_x_discrete(labels = ENSEMBL_ALIASES) +
  scale_color_gradientn(
    colors = c("#4169E1", "#abd9e9", "#abd9e9", "#fee090", "#fee090", "#CD2626"),
    name   = "Expression Level"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
    plot.margin  = unit(c(0, 0, 0, 2), "cm")
  ) +
  ggtitle(paste("Cluster marker expression —", DATASET))

ggsave(
  file.path(RESULTS_DIR, "dotplot_cluster_markers.png"),
  p_dot, width = 20, height = 5, dpi = 300, limitsize = FALSE
)

# ── Cell type annotation ───────────────────────────────────────────────────────
# CELL_TYPE_MAP comes from config.R: named vector of cluster → cell type label.
# To reannotate: update the map in config.R and re-run from this point.
Idents(seurat_combined) <- seurat_combined$seurat_clusters

seurat_combined <- RenameIdents(seurat_combined, CELL_TYPE_MAP)
seurat_combined$cell_type <- Idents(seurat_combined)

# simple_cell_type uses the same labels — kept as a separate column so
# downstream scripts can reference it without depending on the active ident
seurat_combined$simple_cell_type <- as.character(seurat_combined$cell_type)

# Flag tumour cells (defined in config.R) for easy subsetting downstream
seurat_combined$is_tumour_cell <- seurat_combined$simple_cell_type %in% TUMOUR_CELLS

# ── Annotated UMAP ────────────────────────────────────────────────────────────
p_umap <- DimPlot(
  seurat_combined,
  reduction  = "umap",
  label      = TRUE,
  repel      = TRUE,
  label.size = 4.5
) +
  ggtitle(paste("Cell type annotation —", DATASET))

ggsave(
  file.path(RESULTS_DIR, "umap_annotated.png"),
  p_umap, width = 12, height = 8, dpi = 300
)

# UMAP split by condition for visual QC
p_umap_split <- DimPlot(
  seurat_combined,
  reduction = "umap",
  split.by  = "condition",
  label     = TRUE,
  repel     = TRUE,
  label.size = 3
) +
  ggtitle(paste(COND1, "vs", COND2, "—", DATASET))

ggsave(
  file.path(RESULTS_DIR, "umap_annotated_split_condition.png"),
  p_umap_split, width = 14, height = 6, dpi = 300
)

saveRDS(seurat_combined, file = file.path(RESULTS_DIR, "seurat_obj_annotated.rds"))
