#Cluster marker genes###########################
# Marker identification for clusters
DefaultAssay(seurat_combined)
seurat_combined <- PrepSCTFindMarkers(seurat_combined)
DefaultAssay(seurat_combined) <- "SCT"
#Exploration / clustering / plotting UMAP → "integrated"
#DE & biological interpretation → "SCT"
all_markers <- FindAllMarkers(seurat_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(all_markers, file = "all_markers_per_cluster_de.csv")


# Top markers heatmap
top10 <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10, file = "top_10_marker_per_cluster.csv")

# Scale only the top marker genes
seurat_combined <- ScaleData(seurat_combined, features = top10$gene, assay = "SCT")
topclustermarker <-DoHeatmap(seurat_combined, features = top10$gene) + NoLegend()
ggsave("top_cluster_marker.png", topclustermarker, width = 10, height = 15, dpi = 300)

# Identifying cell types ###########
DefaultAssay(seurat_combined) <- "SCT"

cell_types <- list(
  "Endothelial Cells" = c("PECAM1", "FLT1", "KDR", "EMCN"),
  "Neurosecretory-like Cells" = c("SYT1", "EPHA5", "DAB1", "RORB"), 
  "Chromaffin Cells" = c("TH", "DBH", "ENSCAFG00000024864"),
  "Smooth Muscle Cells" = c("ACTA2", "TAGLN", "TPM2", "PRKG1"),
  "Fibroblasts" = c("COL1A2", "DCN", "MFAP5", "SERPINF1"),
  "Schwann cells" = c("MPZ","CDH19","PTPRZ1"),
  "Steroidogenic Cells" = c("STAR", "CYP11A1", "FGFR2", "RBP4"),
  "Immune Cells (APC/Lymphoid Mixed)" = c("PTPRC", "CD86", "IKZF1", "DLA-DQA1"),
  "Perivascular Cells" = c("BMPER", "EDN1", "AQP1", "GJA5"),
  "Mature Medullary Neurons" = c("CNTN3", "CDH4", "IGDCC4", "AFf106"),
  "T Cells" = c("CD3E", "ITK", "SKAP1", "BCL11B")
)

CHGA <- c(
  "ENSCAFG00000024864" = "CHGA")

# Generate a bubble plot using Seurat's built-in function
bubble_annotation <- DotPlot(seurat_combined, features = cell_types, group.by = "seurat_clusters",dot.scale = 6) +
  scale_x_discrete(labels = CHGA) +  # Apply custom labels
  theme_minimal() +
  ggtitle("Bubble Plot of Gene Expression Across Cluster")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.margin = unit(c(0, 0, 0, 2), "cm")  # Top, right, bottom, left      
  )+
  # Add custom gradient with all six colors
  scale_color_gradientn(colors = c("#4169E1", "#abd9e9", "#abd9e9", "#fee090", "#fee090", "#CD2626"),
                        name = "Expression Level")

print(bubble_annotation)
ggsave("bubble_annotation_cell_types_per_cluster.png", bubble_annotation, width = 20, height = 5, dpi = 300, limitsize = FALSE)

# Cell type annotation############
# Set the cluster identities
Idents(seurat_combined) <- seurat_combined$seurat_clusters

descriptive_names <- c(
  "1"  = "Neurosecretory Chromaffin",
  "2"  = "Endothelial",
  "3"  = "Proliferating cells",
  "4"  = "Steroidogenic",
  "5"  = "Chromaffin",
  "6"  = "Fibroblasts",
  "7"  = "Pericytes",
  "8"  = "Macrophages",
  "9"  = "T cells",
  "10" = "Fetal progenitors",
  "11" = "Fetal sympathoblasts"
)

# Simplified names for convenience
simplified_names <- c(
  "1"  = "Neurosecretory_Chromaffin",
  "2"  = "Endothelial",
  "3"  = "Proliferating",
  "4"  = "Steroidogenic",
  "5"  = "Chromaffin",
  "6"  = "Fibroblasts",
  "7"  = "Pericytes",
  "8"  = "Macrophages",
  "9"  = "T",
  "10" = "Fetal_progenitor",
  "11" = "Fetal_sympathoblasts"
)

# Rename for plotting (descriptive)
seurat_combined <- RenameIdents(seurat_combined, descriptive_names)

# Save descriptive names to metadata
seurat_combined$cell_type <- Idents(seurat_combined)

# Save simplified names to metadata
seurat_combined$simple_cell_type <- plyr::mapvalues(
  x = as.character(seurat_combined$seurat_clusters),
  from = names(simplified_names),
  to = simplified_names
)

# Plot the UMAP with descriptive names
umap_annotation <- DimPlot(seurat_combined, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4.5 ) + 
  ggtitle("UMAP with Cell Type Annotations")

ggsave(filename = "umap_annotation_update.png", plot = umap_annotation,
       width = 12, height = 8, dpi = 300, units = "in")

# Optional: Save the Seurat object with both annotations
saveRDS(seurat_combined, file = "seurat_obj_annotated.rds")
