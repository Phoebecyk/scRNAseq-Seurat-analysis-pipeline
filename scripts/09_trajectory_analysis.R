# Load required libraries
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(viridis)
library(pheatmap)
library(RColorBrewer)

# Set seed and theme
set.seed(123)
theme_set(theme_bw(base_size = 14))

#Trajectory Analysis (Monocle3)#####
setwd("Trajectory")
seurat_combined <- readRDS("steroidogenic_subset_processed.rds")

# ---- Monocle3 Trajectory Analysis Pipeline ---- #

# 1. Convert Seurat object to Monocle3's cell_data_set object ####
# Get the count matrix
names(seurat_combined@assays$RNA@layers)

library(Matrix)

# Get only the raw count layers (those starting with "counts.")
layer_names <- grep("^counts\\.", names(seurat_combined@assays$RNA@layers), value = TRUE)

# Load each count matrix
raw_counts_list <- lapply(layer_names, function(layer) {
  LayerData(seurat_combined, assay = "RNA", layer = layer)
})
names(raw_counts_list) <- layer_names

# Find union of all genes
all_genes <- Reduce(union, lapply(raw_counts_list, rownames))

# Make sure all matrices have same rows (genes), fill missing with 0
raw_counts_list <- lapply(raw_counts_list, function(mat) {
  mat_full <- Matrix(0, nrow = length(all_genes), ncol = ncol(mat), sparse = TRUE)
  rownames(mat_full) <- all_genes
  colnames(mat_full) <- colnames(mat)
  mat_full[rownames(mat), ] <- mat
  return(mat_full)
})

# Combine all matrices by column
expression_matrix <- do.call(cbind, raw_counts_list)

# Align metadata
cell_metadata <- seurat_combined@meta.data
cell_metadata <- cell_metadata[colnames(expression_matrix), , drop = FALSE]

# Get the gene metadata
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

# Create the new cell_data_set object
cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

#Check that all cells from all samples are in expression_matrix
# Number of cells in expression matrix
n_cells_matrix <- ncol(expression_matrix)

# Number of cells in Seurat object
n_cells_seurat <- ncol(seurat_combined)

# Compare
cat("Cells in expression_matrix: ", n_cells_matrix, "\n")
cat("Cells in Seurat object:     ", n_cells_seurat, "\n")

# Check if they match
if (setequal(colnames(expression_matrix), colnames(seurat_combined))) {
  cat("All cells from the Seurat object are included in expression_matrix.\n")
} else {
  cat("⚠️ Cell mismatch detected.\n")
}


# Set the UMAP embedding from the Seurat object to the Monocle3 object
# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
list_cluster <- seurat_combined@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- seurat_combined@reductions$umap@cell.embeddings

# plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.before.trajectory

# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 2]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           cell_size = 1.5) +
  labs(title = "Pseudotime-Root: cluster2")

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Extract cluster info from Seurat object (should be a factor)
clusters_seurat <- seurat_combined@active.ident

# Add it to the cds colData
cds@colData$cluster <- clusters_seurat[colnames(cds)]

# Then create the data frame again
data.pseudo <- as.data.frame(colData(cds))

# Check it's there now
str(data.pseudo$cluster)
table(data.pseudo$cluster)

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(cluster, monocle3_pseudotime, median), fill = cluster)) +
  geom_boxplot()

# ...5. Finding genes that change as a function of pseudotime --------------------
deg <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head(30)

sig_deg <- deg %>% 
  arrange(q_value) %>% 
  filter(status == 'OK')
write.csv(sig_deg, "sig_degs_pseudotime.csv", row.names = FALSE)

DefaultAssay(seurat_combined) <- "RNA"
FeaturePlot(seurat_combined, features = c('ENSCAFG00000046992','ENSCAFG00000043730','RPL11', 'RPL18', 'RPL32','UBA52'))


# visualizing pseudotime in seurat

seurat_combined$pseudotime <- pseudotime(cds)
Idents(seurat_combined) <- seurat_combined$seurat_clusters
FeaturePlot(seurat_combined, features = "pseudotime", label = T)

# visualing psuedotime in condition
# Double-check 'condition' exists
head(seurat_combined@meta.data$condition)

# Add condition to Monocle3 object (if not already included)
cds@colData$condition <- seurat_combined@meta.data$condition[colnames(cds)]

plot_cells(
  cds,
  color_cells_by = 'condition',
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  cell_size = 1.5
) +
  labs(title = "Trajectory Colored by Condition") +
  theme(legend.position = "right")  # Ensure legend is visible



# ---- Heatmap Generation Pipeline ---- #
#Extract log-normalized expression
# Get all data.* layers
data_layer_names <- grep("^data\\.", names(seurat_combined@assays$RNA@layers), value = TRUE)

# Load each data matrix
data_list <- lapply(data_layer_names, function(layer) {
  LayerData(seurat_combined, assay = "RNA", layer = layer)
})
names(data_list) <- data_layer_names

# Union of all genes
all_genes <- Reduce(union, lapply(data_list, rownames))

# Standardize matrices to have the same genes (fill missing with 0)
library(Matrix)
data_list <- lapply(data_list, function(mat) {
  mat_full <- Matrix(0, nrow = length(all_genes), ncol = ncol(mat), sparse = TRUE)
  rownames(mat_full) <- all_genes
  colnames(mat_full) <- colnames(mat)
  mat_full[rownames(mat), ] <- mat
  return(mat_full)
})

# Combine into one big expression matrix
lognorm_expression_matrix <- do.call(cbind, data_list)

# Top genes from graph_test
top_pseudotime_genes <- deg %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  pull(gene_short_name) %>%
  head(50) %>% 
  unique()

# Intersect to ensure they exist in expression matrix
top_pseudotime_genes <- intersect(top_pseudotime_genes, rownames(lognorm_expression_matrix))

# Subset
heatmap_data <- as.matrix(lognorm_expression_matrix[top_pseudotime_genes, ])

#Order cells by pseudotime
pseudotime_values <- pseudotime(cds)
pseudotime_values <- pseudotime_values[!is.na(pseudotime_values)]

# Ensure cell order matches pseudotime
heatmap_data <- heatmap_data[, names(pseudotime_values)]
heatmap_data <- heatmap_data[, order(pseudotime_values)]

# Scale expression (row-wise z-score)
scaled_heatmap_data <- t(scale(t(heatmap_data)))

#Column annotations
column_annotations <- data.frame(
  Pseudotime = sort(pseudotime_values),
  Condition = seurat_combined$condition[names(sort(pseudotime_values))]
)
rownames(column_annotations) <- names(sort(pseudotime_values))

#Plot the heatmap
library(pheatmap)
library(RColorBrewer)

my_colour_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  color = my_colour_palette,
  main = "Top 50 Pseudotime-Associated Genes",
  cellwidth = 0.5,
  cellheight = 10,
  width = 6,
  height = 10
)

#Add GO label
# Cut the gene dendrogram into N clusters (e.g., 4–6)
gene_clusters <- cutree(hclust(dist(scaled_heatmap_data)), k = 5)  # Try k = 4 to 6

# Attach to your heatmap
row_annotation <- data.frame(GeneCluster = factor(gene_clusters))
rownames(row_annotation) <- rownames(scaled_heatmap_data)

library(clusterProfiler)
library(org.Hs.eg.db)

# Example: GO enrichment for cluster 1
genes_in_cluster1 <- rownames(row_annotation)[row_annotation$GeneCluster == 1]

# Convert to Entrez IDs (required by enrichGO)
entrez_ids <- bitr(genes_in_cluster1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Cf.eg.db)

# GO enrichment
ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Cf.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# View top pathways
head(ego)
#Repeat for each gene cluster and manually assign labels like:
#Cluster 1: "Immune activation"
#Cluster 2: "Ribosome biogenesis"
#Cluster 3: "Metabolic reprogramming"

#Add the annotations to the heatmap
#Use row_annotation (gene modules) and label each with its dominant GO/pathway manually or in a legend.
# Final annotated heatmap
pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  annotation_row = row_annotation,
  color = my_colour_palette,
  main = "Differentially Expressed Genes Across Pseudotime",
  cellwidth = 0.5,
  cellheight = 5
)


#To confirm that pseudotime captures progressive changes, check distributions:
df <- as.data.frame(colData(cds))
ggplot(df, aes(x = monocle3_pseudotime, fill = condition)) +
  geom_density(alpha = 0.5) +
  labs(title = "Pseudotime Distribution by Condition")

ggplot(df, aes(x = cluster, y = monocle3_pseudotime, fill = condition)) +
  geom_boxplot() +
  labs(title = "Pseudotime by Cluster and Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#NES as root, PCA####
# ---- Monocle3 Trajectory Analysis Pipeline (Revised) ---- #

setwd("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/steroidogenic_cells/Trajectory/NES/PCA")

library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(Matrix)
library(scater)

# 1. Convert Seurat to Monocle3 cell_data_set ---------------------

# Get raw count layers
layer_names <- grep("^counts\\.", names(seurat_combined@assays$RNA@layers), value = TRUE)
raw_counts_list <- lapply(layer_names, function(layer) {
  LayerData(seurat_combined, assay = "RNA", layer = layer)
})
names(raw_counts_list) <- layer_names

# Get union of all genes
all_genes <- Reduce(union, lapply(raw_counts_list, rownames))

# Fill missing genes with zeros
raw_counts_list <- lapply(raw_counts_list, function(mat) {
  mat_full <- Matrix(0, nrow = length(all_genes), ncol = ncol(mat), sparse = TRUE)
  rownames(mat_full) <- all_genes
  colnames(mat_full) <- colnames(mat)
  mat_full[rownames(mat), ] <- mat
  return(mat_full)
})

# Combine expression matrices
expression_matrix <- do.call(cbind, raw_counts_list)

# Align metadata
cell_metadata <- seurat_combined@meta.data
cell_metadata <- cell_metadata[colnames(expression_matrix), , drop = FALSE]

# Gene metadata
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

# Create cds
cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# 2. Preprocess and Reduce Dim (use PCA instead of UMAP) ------------------

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "PCA")

# 3. Cluster (optional; just for plotting)
cds <- cluster_cells(cds, reduction_method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 4. Learn graph (on PCA)
cds <- learn_graph(cds, use_partition = FALSE)

# 5. Use top NES-expressing cells as root -------------------------------

# Get NES expression from normalized data
normalized_counts <- logNormCounts(cds)
exprs <- assay(normalized_counts, "logcounts")  # log-normalized counts matrix
nes_expr <- exprs["NES", ]  # NES gene expression vector for all cells
top_nes_cells <- names(sort(nes_expr, decreasing = TRUE))[1:5]

# Order cells using top NES cells
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = top_nes_cells)

# 6. Plot pseudotime ----------------------------------------------------

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           cell_size = 1.5) +
  labs(title = "Pseudotime Rooted in NES-high Cells (PCA space)")

# 7. Save pseudotime + extract metadata ------------------------------

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Optionally add Seurat cluster or condition
cds$cluster <- seurat_combined$seurat_clusters[colnames(cds)]
cds$condition <- seurat_combined$condition[colnames(cds)]

# Pseudotime by condition
ggplot(data.pseudo, aes(monocle3_pseudotime, fill = condition)) +
  geom_density(alpha = 0.5) +
  labs(title = "Pseudotime Distribution by Condition")

# Pseudotime by cluster
# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Extract cluster info from Seurat object (should be a factor)
clusters_seurat <- seurat_combined@active.ident

# Add it to the cds colData
cds@colData$cluster <- clusters_seurat[colnames(cds)]

# Then create the data frame again
data.pseudo <- as.data.frame(colData(cds))

# Check it's there now
str(data.pseudo$cluster)
table(data.pseudo$cluster)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(cluster, monocle3_pseudotime, median), fill = cluster)) +
  geom_boxplot() +
  labs(title = "Pseudotime by Cluster")

# Write Seurat clusters into Monocle’s cluster slot (used by plot_cells)
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "cluster",
           label_groups_by_cluster = TRUE,  # Now works correctly
           cell_size = 1
) +
  labs(title = "Trajectory Colored by UMAP Cluster") +
  theme(legend.position = "right")  # Ensure legend is visible

ggplot(data.pseudo, aes(condition, monocle3_pseudotime, fill = condition)) +
  geom_boxplot() +
  labs(title = "Pseudotime distribution by Condition") +
  theme_minimal()

plot_cells(cds,
           genes = "NES",
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           cell_size = 1.5) +
  labs(title = "NES Expression Along Trajectory")


# 8. Find DEGs across pseudotime ---------------------------------------

deg <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

sig_deg <- deg %>%
  arrange(q_value) %>%
  filter(status == "OK")

write.csv(sig_deg, "sig_degs_pseudotime_NESroot_PCA.csv", row.names = FALSE)

deg %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head(30)

DefaultAssay(seurat_combined) <- "RNA"
FeaturePlot(seurat_combined, features = c('SLC25A6','SSR4','NHS', 'NSCAFG00000045333', 'ENSCAFG00000017688'))
#ENSCAFG00000045333 = RPL10
#ENSCAFG00000017688 = RPL36A 

# visualizing pseudotime in seurat

seurat_combined$pseudotime <- pseudotime(cds)
Idents(seurat_combined) <- seurat_combined$seurat_clusters
FeaturePlot(seurat_combined, features = "pseudotime", label = T)

# visualing psuedotime in condition
# Double-check 'condition' exists
head(seurat_combined@meta.data$condition)

# Add condition to Monocle3 object (if not already included)
cds@colData$condition <- seurat_combined@meta.data$condition[colnames(cds)]
cds@colData$condition <- seurat_combined@meta.data[rownames(colData(cds)), "condition"]


plot_cells(
  cds,
  color_cells_by = 'condition',
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  cell_size = 1
) +
  labs(title = "Trajectory Colored by Condition") +
  theme(legend.position = "right")  # Ensure legend is visible

# ---- Heatmap Generation Pipeline ---- #
#Extract log-normalized expression
# Get all data.* layers
data_layer_names <- grep("^data\\.", names(seurat_combined@assays$RNA@layers), value = TRUE)

# Load each data matrix
data_list <- lapply(data_layer_names, function(layer) {
  LayerData(seurat_combined, assay = "RNA", layer = layer)
})
names(data_list) <- data_layer_names

# Union of all genes
all_genes <- Reduce(union, lapply(data_list, rownames))

# Standardize matrices to have the same genes (fill missing with 0)
library(Matrix)
data_list <- lapply(data_list, function(mat) {
  mat_full <- Matrix(0, nrow = length(all_genes), ncol = ncol(mat), sparse = TRUE)
  rownames(mat_full) <- all_genes
  colnames(mat_full) <- colnames(mat)
  mat_full[rownames(mat), ] <- mat
  return(mat_full)
})

# Combine into one big expression matrix
lognorm_expression_matrix <- do.call(cbind, data_list)

# Top genes from graph_test
top_pseudotime_genes <- deg %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  pull(gene_short_name) %>%
  head(50) %>% 
  unique()

# Intersect to ensure they exist in expression matrix
top_pseudotime_genes <- intersect(top_pseudotime_genes, rownames(lognorm_expression_matrix))

# Subset
heatmap_data <- as.matrix(lognorm_expression_matrix[top_pseudotime_genes, ])

#Order cells by pseudotime
pseudotime_values <- pseudotime(cds)
pseudotime_values <- pseudotime_values[!is.na(pseudotime_values)]

# Ensure cell order matches pseudotime
heatmap_data <- heatmap_data[, names(pseudotime_values)]
heatmap_data <- heatmap_data[, order(pseudotime_values)]

# Scale expression (row-wise z-score)
scaled_heatmap_data <- t(scale(t(heatmap_data)))

#Column annotations
column_annotations <- data.frame(
  Pseudotime = sort(pseudotime_values),
  Condition = seurat_combined$condition[names(sort(pseudotime_values))]
)
rownames(column_annotations) <- names(sort(pseudotime_values))

#Plot the heatmap
library(pheatmap)
library(RColorBrewer)

my_colour_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  color = my_colour_palette,
  main = "Top 50 Pseudotime-Associated Genes",
  cellwidth = 0.5,
  cellheight = 10,
  width = 6,
  height = 10
)

#forcus rearrange normal first
# Get pseudotime and condition info
pseudotime_values <- pseudotime(cds)
pseudotime_values <- pseudotime_values[!is.na(pseudotime_values)]

# Get condition info aligned to pseudotime cells
condition_vec <- seurat_combined$condition[names(pseudotime_values)]

# Create a data frame for sorting
ordering_df <- data.frame(
  pseudotime = pseudotime_values,
  condition = factor(condition_vec, levels = c("Normal", "csACT"))  # Ensures Normal comes first
)

# Sort by condition first, then pseudotime
ordering_df <- ordering_df[order(ordering_df$condition, ordering_df$pseudotime), ]

# Reorder the heatmap data accordingly
heatmap_data <- heatmap_data[, rownames(ordering_df)]

# Recalculate scaled expression
scaled_heatmap_data <- t(scale(t(heatmap_data)))

# Create updated column annotations
column_annotations <- data.frame(
  Pseudotime = ordering_df$pseudotime,
  Condition = ordering_df$condition
)
rownames(column_annotations) <- rownames(ordering_df)

pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  color = my_colour_palette,
  main = "Top 50 Pseudotime-Associated Genes",
  cellwidth = 0.5,
  cellheight = 10,
  width = 6,
  height = 10
)

#see monoclone independent cluster
cds@clusters$UMAP$clusters
cds$monocle3_cluster <- cds@clusters$UMAP$clusters
table(cds$monocle3_cluster)  # Shows how many cells in each cluster

plot_cells(cds,
           color_cells_by = "monocle3_cluster",
           label_groups_by_cluster = TRUE,
           cell_size = 1) +
  labs(title = "Monocle3 UMAP Clusters")


#NES as root, UMAP####
# ---- Monocle3 Trajectory Analysis Pipeline (Revised) ---- #

setwd("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/steroidogenic_cells/Trajectory/NES/UMAP")

# 1. Convert Seurat object to Monocle3's cell_data_set object
# Get the count matrix
names(seurat_combined@assays$RNA@layers)

library(Matrix)

# Get only the raw count layers (those starting with "counts.")
layer_names <- grep("^counts\\.", names(seurat_combined@assays$RNA@layers), value = TRUE)

# Load each count matrix
raw_counts_list <- lapply(layer_names, function(layer) {
  LayerData(seurat_combined, assay = "RNA", layer = layer)
})
names(raw_counts_list) <- layer_names

# Find union of all genes
all_genes <- Reduce(union, lapply(raw_counts_list, rownames))

# Make sure all matrices have same rows (genes), fill missing with 0
raw_counts_list <- lapply(raw_counts_list, function(mat) {
  mat_full <- Matrix(0, nrow = length(all_genes), ncol = ncol(mat), sparse = TRUE)
  rownames(mat_full) <- all_genes
  colnames(mat_full) <- colnames(mat)
  mat_full[rownames(mat), ] <- mat
  return(mat_full)
})

# Combine all matrices by column
expression_matrix <- do.call(cbind, raw_counts_list)

# Align metadata
cell_metadata <- seurat_combined@meta.data
cell_metadata <- cell_metadata[colnames(expression_matrix), , drop = FALSE]

# Get the gene metadata
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

# Create the new cell_data_set object
cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

#Check that all cells from all samples are in expression_matrix
# Number of cells in expression matrix
n_cells_matrix <- ncol(expression_matrix)

# Number of cells in Seurat object
n_cells_seurat <- ncol(seurat_combined)

# Compare
cat("Cells in expression_matrix: ", n_cells_matrix, "\n")
cat("Cells in Seurat object:     ", n_cells_seurat, "\n")

# Check if they match
if (setequal(colnames(expression_matrix), colnames(seurat_combined))) {
  cat("All cells from the Seurat object are included in expression_matrix.\n")
} else {
  cat("⚠️ Cell mismatch detected.\n")
}


# Set the UMAP embedding from the Seurat object to the Monocle3 object
# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
list_cluster <- seurat_combined@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- seurat_combined@reductions$umap@cell.embeddings

# plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.before.trajectory

#3. Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

# Get NES expression from normalized data
normalized_counts <- logNormCounts(cds)
exprs <- assay(normalized_counts, "logcounts")  # log-normalized counts matrix
nes_expr <- exprs["NES", ]  # NES gene expression vector for all cells
top_nes_cells <- names(sort(nes_expr, decreasing = TRUE))[1:5]

# Order cells using top NES cells
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = top_nes_cells)

# 6. Plot pseudotime 

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           cell_size = 1.5) +
  labs(title = "Pseudotime Rooted in NES-high Cells")

# 7. Save pseudotime + extract metadata 

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Optionally add Seurat cluster or condition
cds$cluster <- seurat_combined$seurat_clusters[colnames(cds)]
cds$condition <- seurat_combined$condition[colnames(cds)]

# Pseudotime by condition
ggplot(data.pseudo, aes(monocle3_pseudotime, fill = condition)) +
  geom_density(alpha = 0.5) +
  labs(title = "Pseudotime Distribution by Condition")

# Pseudotime by cluster
# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Extract cluster info from Seurat object (should be a factor)
clusters_seurat <- seurat_combined@active.ident

# Add it to the cds colData
cds@colData$cluster <- clusters_seurat[colnames(cds)]

# Then create the data frame again
data.pseudo <- as.data.frame(colData(cds))

# Check it's there now
str(data.pseudo$cluster)
table(data.pseudo$cluster)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(cluster, monocle3_pseudotime, median), fill = cluster)) +
  geom_boxplot() +
  labs(title = "Pseudotime by Cluster")

# 8. Find DEGs across pseudotime

deg <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

sig_deg <- deg %>%
  arrange(q_value) %>%
  filter(status == "OK")

write.csv(sig_deg, "sig_degs_pseudotime_NESroot_UMAP.csv", row.names = FALSE)

deg %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head(30)

deg_filtered <- sig_deg %>%
  filter(status == 'OK') %>%
  arrange(q_value) %>%
  filter(!grepl("^RPS|^RPL|^UBA52|^FAU", rownames(.))) %>%
  head(50)
write.csv(deg_filtered, "deg_filtered.csv", row.names = FALSE)


DefaultAssay(seurat_combined) <- "RNA"
FeaturePlot(seurat_combined, features = c('ENSCAFG00000046992','ENSCAFG00000043730','RPL11', 'RPL18', 'RPL32'))

# visualizing pseudotime in seurat

seurat_combined$pseudotime <- pseudotime(cds)
Idents(seurat_combined) <- seurat_combined$seurat_clusters
FeaturePlot(seurat_combined, features = "pseudotime", label = T)

# visualing psuedotime in condition
# Double-check 'condition' exists
head(seurat_combined@meta.data$condition)

# Add condition to Monocle3 object (if not already included)
cds@colData$condition <- seurat_combined@meta.data$condition[colnames(cds)]
cds@colData$condition <- seurat_combined@meta.data[rownames(colData(cds)), "condition"]


plot_cells(
  cds,
  color_cells_by = 'condition',
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  cell_size = 1
) +
  labs(title = "Trajectory Colored by Condition") +
  theme(legend.position = "right")  # Ensure legend is visible

# ---- Heatmap Generation Pipeline ---- #
#Extract log-normalized expression
# Get all data.* layers
data_layer_names <- grep("^data\\.", names(seurat_combined@assays$RNA@layers), value = TRUE)

# Load each data matrix
data_list <- lapply(data_layer_names, function(layer) {
  LayerData(seurat_combined, assay = "RNA", layer = layer)
})
names(data_list) <- data_layer_names

# Union of all genes
all_genes <- Reduce(union, lapply(data_list, rownames))

# Standardize matrices to have the same genes (fill missing with 0)
library(Matrix)
data_list <- lapply(data_list, function(mat) {
  mat_full <- Matrix(0, nrow = length(all_genes), ncol = ncol(mat), sparse = TRUE)
  rownames(mat_full) <- all_genes
  colnames(mat_full) <- colnames(mat)
  mat_full[rownames(mat), ] <- mat
  return(mat_full)
})

# Combine into one big expression matrix
lognorm_expression_matrix <- do.call(cbind, data_list)

# Top genes from graph_test
top_pseudotime_genes <- deg_filtered %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  pull(gene_short_name) %>%
  head(50) %>% 
  unique()

# Intersect to ensure they exist in expression matrix
top_pseudotime_genes <- intersect(top_pseudotime_genes, rownames(lognorm_expression_matrix))

# Subset
heatmap_data <- as.matrix(lognorm_expression_matrix[top_pseudotime_genes, ])

#Order cells by pseudotime
pseudotime_values <- pseudotime(cds)
pseudotime_values <- pseudotime_values[!is.na(pseudotime_values)]

# Ensure cell order matches pseudotime
heatmap_data <- heatmap_data[, names(pseudotime_values)]
heatmap_data <- heatmap_data[, order(pseudotime_values)]

# Scale expression (row-wise z-score)
scaled_heatmap_data <- t(scale(t(heatmap_data)))

#Column annotations
column_annotations <- data.frame(
  Pseudotime = sort(pseudotime_values),
  Condition = seurat_combined$condition[names(sort(pseudotime_values))]
)
rownames(column_annotations) <- names(sort(pseudotime_values))

#Plot the heatmap
library(pheatmap)
library(RColorBrewer)

my_colour_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
ann_colors <- list(
  simple_cell_type = c(
    "Normal" = "salmon",
    "csACT" = "lightblue"
  )
)
pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  annotation_colors = ann_colors,
  color = my_colour_palette,
  main = "Top 50 Pseudotime-Associated Genes",
  cellwidth = 0.5,
  cellheight = 10,
  width = 6,
  height = 10
)

#forcus rearrange normal first
# Get pseudotime and condition info
pseudotime_values <- pseudotime(cds)
pseudotime_values <- pseudotime_values[!is.na(pseudotime_values)]

# Get condition info aligned to pseudotime cells
condition_vec <- seurat_combined$condition[names(pseudotime_values)]

# Create a data frame for sorting
ordering_df <- data.frame(
  pseudotime = pseudotime_values,
  condition = factor(condition_vec, levels = c("Normal", "csACT"))  # Ensures Normal comes first
)

# Sort by condition first, then pseudotime
ordering_df <- ordering_df[order(ordering_df$condition, ordering_df$pseudotime), ]

# Reorder the heatmap data accordingly
heatmap_data <- heatmap_data[, rownames(ordering_df)]

# Recalculate scaled expression
scaled_heatmap_data <- t(scale(t(heatmap_data)))

# Create updated column annotations
column_annotations <- data.frame(
  Pseudotime = ordering_df$pseudotime,
  Condition = ordering_df$condition
)
rownames(column_annotations) <- rownames(ordering_df)

pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  color = my_colour_palette,
  main = "Top 50 Pseudotime-Associated Genes",
  cellwidth = 0.5,
  cellheight = 10,
  width = 6,
  height = 10
)



#LIPA####

setwd("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/steroidogenic_cells/Trajectory/LIPA")

# 1. Convert Seurat object to Monocle3's cell_data_set object
# Get the count matrix
seurat_combined <- JoinLayers(seurat_combined, assay = "RNA")
counts_matrix <- GetAssayData(seurat_combined, assay = "RNA", layer = "counts")

# Extract cell metadata from the Seurat object.
# Monocle3 requires cell names to be in the row names of the metadata data frame.
cell_metadata <- seurat_combined@meta.data

# Extract gene metadata.
# Monocle3 requires gene short names, which we'll set as the gene symbols.
gene_metadata <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix)
)

# Create the Monocle3 cell_data_set (CDS) object.
cds <- new_cell_data_set(
  expression_data = counts_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Set the UMAP embedding from the Seurat object to the Monocle3 object
# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
list_cluster <- seurat_combined@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- seurat_combined@reductions$umap@cell.embeddings

# plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.before.trajectory

#3. Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

# Get LIPA
# (Optional but recommended) Visualise a key normal cell marker to confirm its location.
plot_cells(cds,
           genes = c("LIPA"), 
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           cell_size = 1.5)

cell_type_column   <- "condition"             # The metadata column with cell type annotations
normal_cell_type   <- "Normal" # The name of your normal/root cell population
root_gene          <- "LIPA"                  # The gene to identify the earliest cells

# --- Step 5b: Check that inputs exist ---
if(!root_gene %in% rowData(cds)$gene_short_name) {
  stop(paste("Error: Root gene '", root_gene, "' not found in the dataset.", sep=""))
}
if (!normal_cell_type %in% colData(cds)[[cell_type_column]]) {
  stop(paste("Error: The normal cell type '", normal_cell_type, "' was not found in the '", cell_type_column, "' metadata column.", sep=""))
}

# --- Step 5c: Identify the Root Cells ---
# Get the barcodes of all cells belonging to the normal cell type.
all_normal_cells <- rownames(colData(cds)[colData(cds)[[cell_type_column]] == normal_cell_type,])

# Get the expression values for the root gene across ALL cells.
lipa_expression_matrix <- cds[rowData(cds)$gene_short_name == root_gene, ]
full_expression_vec <- as.vector(exprs(lipa_expression_matrix))
names(full_expression_vec) <- colnames(exprs(lipa_expression_matrix))

# Filter the expression vector to include ONLY the normal cells.
normal_lipa_expression_vec <- full_expression_vec[names(full_expression_vec) %in% all_normal_cells]

# From this subset, find the cells with the highest expression.
# We'll select the top 1% of expressing normal cells to be the root.
quantile_threshold <- quantile(normal_lipa_expression_vec, 0.99)
root_cell_barcodes <- names(normal_lipa_expression_vec[normal_lipa_expression_vec >= quantile_threshold])

# Check if any cells were selected as roots.
if (length(root_cell_barcodes) == 0) {
  # Fallback if no cells are above the 99th percentile.
  root_cell_barcodes <- names(sort(normal_lipa_expression_vec, decreasing = TRUE)[1:10])
  warning(paste("No '", normal_cell_type, "' cells found above 99th percentile of ", root_gene, " expression. Using top 10 expressing cells as the root.", sep=""))
}

# --- Step 5d: Order the full dataset using the identified root cells ---
cds <- order_cells(cds, root_cells = root_cell_barcodes)

# 6. Plot pseudotime 

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           cell_size = 1) +
  labs(title = "Pseudotime Rooted in LIPA-high Normal Cells")

# 7. Save pseudotime + extract metadata 

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Optionally add Seurat cluster or condition
cds$cluster <- seurat_combined$seurat_clusters[colnames(cds)]
cds$condition <- seurat_combined$condition[colnames(cds)]

# Pseudotime by condition
ggplot(data.pseudo, aes(monocle3_pseudotime, fill = condition)) +
  geom_density(alpha = 0.5) +
  labs(title = "Pseudotime Distribution by Condition")

# Pseudotime by cluster
# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Extract cluster info from Seurat object (should be a factor)
clusters_seurat <- seurat_combined@active.ident

# Add it to the cds colData
cds@colData$cluster <- clusters_seurat[colnames(cds)]

# Then create the data frame again
data.pseudo <- as.data.frame(colData(cds))

# Check it's there now
str(data.pseudo$cluster)
table(data.pseudo$cluster)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(cluster, monocle3_pseudotime, median), fill = cluster)) +
  geom_boxplot() +
  labs(title = "Pseudotime by Cluster")

# 8. Find DEGs across pseudotime

deg <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

sig_deg <- deg %>%
  arrange(q_value) %>%
  filter(status == "OK")

write.csv(sig_deg, "sig_degs_pseudotime_LIPAroot_UMAP.csv", row.names = FALSE)

deg %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head(30)

deg_filtered <- sig_deg %>%
  filter(status == 'OK') %>%
  arrange(q_value) %>%
  filter(!grepl("^RPS|^RPL|^UBA52|^FAU", rownames(.))) %>%
  head(50)
write.csv(deg_filtered, "deg_filtered.csv", row.names = FALSE)


DefaultAssay(seurat_combined) <- "RNA"
FeaturePlot(seurat_combined, features = c('ENSCAFG00000046992','ENSCAFG00000043730','RPL11', 'RPL18', 'RPL32'))

# visualizing pseudotime in seurat

seurat_combined$pseudotime <- pseudotime(cds)
Idents(seurat_combined) <- seurat_combined$seurat_clusters
FeaturePlot(seurat_combined, features = "pseudotime", label = T)

# visualing psuedotime in condition
# Double-check 'condition' exists
head(seurat_combined@meta.data$condition)

# Add condition to Monocle3 object (if not already included)
cds@colData$condition <- seurat_combined@meta.data$condition[colnames(cds)]
cds@colData$condition <- seurat_combined@meta.data[rownames(colData(cds)), "condition"]


plot_cells(
  cds,
  color_cells_by = 'condition',
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  cell_size = 1
) +
  labs(title = "Trajectory Colored by Condition") +
  theme(legend.position = "right")  # Ensure legend is visible

# ---- Heatmap Generation Pipeline ---- #
#Extract log-normalized expression
lognorm_expression_matrix <- GetAssayData(seurat_combined, assay = "SCT", layer = "data")

# 2. Get top genes from your differential expression results
# Ensure 'deg_filtered' is your data frame from graph_test()
# IMPORTANT: Check that the gene names in `pull()` match the rownames of the matrix.
# For example, use pull(gene_id) if your matrix has Ensembl IDs.
top_pseudotime_genes <- deg_filtered %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  pull(gene_short_name) %>% # <-- This is the most likely source of error.
  head(50) %>% 
  unique()

# --- Sanity Check (Highly Recommended) ---
# Compare the gene name formats before intersecting
print("Gene names from DEG results:")
print(head(top_pseudotime_genes))

print("Gene names from expression matrix:")
print(head(rownames(lognorm_expression_matrix)))
# --- End Sanity Check ---

# 3. Intersect to find common genes
genes_for_heatmap <- intersect(top_pseudotime_genes, rownames(lognorm_expression_matrix))

# Check if you have any genes left
if (length(genes_for_heatmap) == 0) {
  stop("No overlapping genes found between DEG results and expression matrix. Check your gene identifiers (e.g., gene symbols vs. Ensembl IDs).")
}

# 4. Subset the matrix and create the heatmap
heatmap_data <- as.matrix(lognorm_expression_matrix[genes_for_heatmap, ])

#Order cells by pseudotime
pseudotime_values <- pseudotime(cds)
pseudotime_values <- pseudotime_values[!is.na(pseudotime_values)]

# Ensure cell order matches pseudotime
heatmap_data <- heatmap_data[, names(pseudotime_values)]
heatmap_data <- heatmap_data[, order(pseudotime_values)]

# Scale expression (row-wise z-score)
scaled_heatmap_data <- t(scale(t(heatmap_data)))

#Column annotations
column_annotations <- data.frame(
  Pseudotime = sort(pseudotime_values),
  Condition = seurat_combined$condition[names(sort(pseudotime_values))]
)
rownames(column_annotations) <- names(sort(pseudotime_values))

#Plot the heatmap
library(pheatmap)
library(RColorBrewer)

my_colour_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
ann_colors <- list(
  simple_cell_type = c(
    "Normal" = "salmon",
    "csACT" = "lightblue"
  )
)
pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  annotation_colors = ann_colors,
  color = my_colour_palette,
  main = "Top 50 Pseudotime-Associated Genes",
  cellwidth = 0.5,
  cellheight = 10,
  width = 6,
  height = 10
)

#forcus rearrange normal first
# Get pseudotime and condition info
pseudotime_values <- pseudotime(cds)
pseudotime_values <- pseudotime_values[!is.na(pseudotime_values)]

# Get condition info aligned to pseudotime cells
condition_vec <- seurat_combined$condition[names(pseudotime_values)]

# Create a data frame for sorting
ordering_df <- data.frame(
  pseudotime = pseudotime_values,
  condition = factor(condition_vec, levels = c("Normal", "csACT"))  # Ensures Normal comes first
)

# Sort by condition first, then pseudotime
ordering_df <- ordering_df[order(ordering_df$condition, ordering_df$pseudotime), ]

# Reorder the heatmap data accordingly
heatmap_data <- heatmap_data[, rownames(ordering_df)]

# Recalculate scaled expression
scaled_heatmap_data <- t(scale(t(heatmap_data)))

# Create updated column annotations
column_annotations <- data.frame(
  Pseudotime = ordering_df$pseudotime,
  Condition = ordering_df$condition
)
rownames(column_annotations) <- rownames(ordering_df)

pheatmap(
  scaled_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = column_annotations,
  color = my_colour_palette,
  main = "Top 50 Pseudotime-Associated Genes",
  cellwidth = 0.5,
  cellheight = 10,
  width = 6,
  height = 10
)

#Slingshot####
#Trajectory Analysis (Monocle3)#####
setwd("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/steroidogenic_cells/Trajectory/slingshot")
seurat_combined <- readRDS("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/steroidogenic_cells/steroidogenic_subset_processed.rds")
# Required packages
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tradeSeq)   # optional: for testing genes along lineages
library(mgcv)       # for GAM if you prefer

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_combined, assay = "SCT")

# Use PCA embedding (from Seurat) OR diffusion map (destiny)
# Option A: PCA
pca_emb <- Embeddings(seurat_combined, "pca")[, 1:20]   # adjust PCs if you prefer
reducedDim(sce, "PCA") <- pca_emb

# Option B (recommended to try both): diffusion map via destiny
if (!requireNamespace("destiny", quietly = TRUE)) {
  BiocManager::install("destiny")
}
library(destiny)
dm <- DiffusionMap(pca_emb)
dmap_emb <- as.matrix(dm@eigenvectors[, 1:10])   # top diffusion components
reducedDim(sce, "DM") <- dmap_emb

# CLUSTER labels for slingshot: use your Seurat clusters (numeric)
# make sure cluster labels are integers starting at 1
cluster_labels <- as.numeric(as.character(seurat_combined$seurat_clusters))
table(cluster_labels)
colData(sce)$cluster <- cluster_labels

# fit Slingshot on PCA (no start specified first)
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'PCA', 
                 start.clus = NULL, stretch = 0)  # stretch=0 avoids line extension beyond data

# If you want to fit on diffusion map instead:
sce_dm <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'DM', start.clus = NULL, stretch = 0)

# Inspect lineages found
lineages_pca <- slingLineages(sce)
lineages_dm  <- slingLineages(sce_dm)
print(lineages_pca)
print(lineages_dm)

# Extract pseudotimes (one pseudotime per lineage)
pt_pca <- slingPseudotime(sce)
pt_dm  <- slingPseudotime(sce_dm)

# Check dimensions first
dim(pt_pca)   # rows = cells, cols = lineages

# View first few rows of all available lineages
head(pt_pca)

# Same for diffusion map run
dim(pt_dm)
head(pt_dm)

# Add primary pseudotime (choose lineage 1 or combine) back to Seurat metadata.
# Common approach: if only one lineage, take that. If multiple, you can take the min or use lineage-specific analyses.
# Example: take first lineage pseudotime:
seurat_combined$slingshot_pt_lineage1 <- pt_pca[,1]
seurat_combined$slingshot_pt_lineage2 <- if (ncol(pt_pca) >= 2) pt_pca[,2] else NA

# Visualize pseudotime on UMAP (if UMAP exists in seurat_sub)
p1 <- DimPlot(seurat_combined, reduction = "umap", group.by = "condition")
p2 <- FeaturePlot(seurat_combined, features = "slingshot_pt_lineage1", reduction = "umap") + ggtitle("Slingshot pseudotime (lineage1)")
p1; p2

# Plot slingshot curves on PCA space
plot(reducedDim(sce, "PCA")[,1], reducedDim(sce, "PCA")[,2],
     col = as.factor(seurat_combined$condition), pch = 16, cex = 0.6,
     xlab = "PC1", ylab = "PC2", main = "Slingshot on PCA")
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')  # draws trajectories

# If using DM:
plot(reducedDim(sce_dm, "DM")[,1], reducedDim(sce_dm, "DM")[,2],
     col = as.factor(seurat_combined$condition), pch = 16, cex = 0.6,
     xlab = "DC1", ylab = "DC2", main = "Slingshot on Diffusion Map")
lines(SlingshotDataSet(sce_dm), lwd = 2, col = 'black')

# Compare pseudotime distributions between Normal and csACT
library(ggpubr)
df_pt <- data.frame(cell = colnames(sce),
                    pt1 = pt_pca[,1],
                    condition = seurat_combined$condition)
ggplot(df_pt, aes(x = condition, y = pt1)) +
  geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.4) +
  labs(title = "Pseudotime (lineage1) by condition")

# Correlate with Monocle3 pseudotime if you have it:
if ("monocle3_pseudotime" %in% colnames(seurat_combined@meta.data)) {
  df_cor <- data.frame(monocle_pt = seurat_combined$monocle3_pseudotime,
                       sling_pt = seurat_combined$slingshot_pt_lineage1)
  cor_val <- cor(df_cor$monocle_pt, df_cor$sling_pt, use = "complete.obs", method = "spearman")
  cat("Spearman correlation between Monocle3 PT and Slingshot PT (lineage1):", cor_val, "\n")
}

# Test genes along pseudotime: simple GAM approach for a gene list
genes_to_test <- c("ZNF157")  # update with genes of interest
res_gam <- data.frame(gene = genes_to_test, p.value = NA)
for (i in seq_along(genes_to_test)) {
  g <- genes_to_test[i]
  expr <- as.numeric(GetAssayData(seurat_combined, assay = "SCT", slot = "data")[g, ])
  df <- data.frame(expr = expr, pt = seurat_combined$slingshot_pt_lineage1)
  # remove NAs
  df <- df[!is.na(df$pt), ]
  if (nrow(df) > 20) {
    fit <- mgcv::gam(expr ~ s(pt, bs = "cs"), data = df)
    res_gam$p.value[i] <- summary(fit)$s.table[1, 4]
  } else {
    res_gam$p.value[i] <- NA
  }
}
res_gam

library(ggplot2)

# Extract metadata into a data frame
meta <- seurat_combined@meta.data

# Make sure the pseudotime and cell type columns exist
head(meta[, c("slingshot_pt_lineage1", "condition")])

# Plot pseudotime distribution by cell type
ggplot(meta, aes(x = cell_type, y = slingshot_pt_lineage1, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Distribution of Slingshot pseudotime by Condition",
    x = "Cell Type",
    y = "Pseudotime (lineage 1)"
  ) +
  scale_fill_manual(values = c("Normal" = "skyblue", "csACT" = "salmon"))
