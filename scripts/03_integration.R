# Normalise and scale and combine #####

# Normalize and find variable features for each sample
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst")
  return(obj)
})

# Scale data for each sample
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- ScaleData(obj, features = rownames(obj))
  return(obj)
})

# PCA for each sample
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  return(obj)
})

# Cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Cell cycle scoring for each sample
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  return(obj)
})

# Regress out cell cycle effects
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(obj))
  return(obj)
})

# SCTransform normalization for each sample
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- SCTransform(obj, vars.to.regress = c("S.Score", "G2M.Score"), variable.features.n = 3000)
  return(obj)
})

# Prepare for integration
insulinoma.features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = insulinoma.features)

# Find integration anchors
insulinoma.anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                             normalization.method = "SCT", 
                                             anchor.features = insulinoma.features, 
                                             dims = 1:20, 
                                             k.anchor = 5)

##Integrate datasets####
seurat_combined <- IntegrateData(anchorset = insulinoma.anchors, normalization.method = "SCT")

# Proceed with PCA, UMAP, and clustering
seurat_combined <- RunPCA(seurat_combined, npcs = 20)
elbow <- ElbowPlot(seurat_combined)
ggsave("elbow.png", elbow, width = 7, height = 7, dpi = 300)
seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:17)
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:17)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.3, algorithm = 2)

# Plot UMAP to check if balancing affected clustering
DimPlot(seurat_combined, group.by = "condition", reduction = "umap") + ggtitle("Balanced Primary vs Metastasis")


# Convert cluster labels from 0-based to 1-based indexing
seurat_combined$seurat_clusters <- as.factor(as.numeric(as.character(seurat_combined$seurat_clusters)) + 1)
Idents(seurat_combined) <- seurat_combined$seurat_clusters

# Visualize clusters and conditions
umap1 <-DimPlot(seurat_combined, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + ggtitle("Clusters")
umap2<-DimPlot(seurat_combined, reduction = "umap", group.by = "condition") + ggtitle("Conditions")
umap3<-DimPlot(seurat_combined, reduction = "umap", group.by = "sample") + ggtitle("Samples")
umap_plot <- umap1+ umap2 +umap3 + plot_layout(ncol = 3)
show(umap_plot)
ggsave("umap_unanno.png", umap1, width = 7, height = 5, dpi = 300)
ggsave("umap.png", umap_plot, width = 21, height = 7, dpi = 300)
# Save the integrated object
saveRDS(seurat_combined, file = "seurat_combined.rds")
