# Load required libraries
library(Seurat)
library(infercnv)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(ggpubr)

set.seed(123)
theme_set(theme_bw(base_size = 14))

# Set working directory and load your Seurat object
setwd("cnv")
seurat_combined <- readRDS("/steroidogenic_cells.rds")


# if plot not show, close all graphics devices
while (dev.cur() > 1) dev.off()

# 1. Extract Raw Counts Matrix
seurat_combined<- JoinLayers(seurat_combined, assay = "RNA")
counts_matrix <- GetAssayData(seurat_combined, assay = "RNA", slot = "counts")

# Verify dimensions
dim(counts_matrix)
colnames(counts_matrix) # Should match all cell barcodes

# 2. Create Annotations File
annotations <- data.frame(
  CellID = colnames(seurat_combined),
  Sample = seurat_combined$sample
)

write.table(annotations, file = "infercnv_annotations_steroidogenic.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 4. Create InferCNV Object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = "infercnv_annotations_steroidogenic.txt",
  gene_order_file = "filtered_gene_order.pos",
  ref_group_names = c("NAD2016","NAD2017")
)

# 5. Run InferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # 10x Genomics data
  out_dir = "infercnv_output",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  noise_filter = 0.1,
  HMM = FALSE
)

#How InferCNV Clustered the Cells
# Load cell cluster assignments
# Reload with no col names, and check structure
cell_clusters <- read.table(
  "infercnv_output/infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt",
  header = FALSE,
  skip = 1,
  stringsAsFactors = FALSE
)
colnames(cell_clusters) <- c("CellID", "Cluster", "Color1", "Number", "Color2")
cell_clusters_sub <- cell_clusters[, c("CellID", "Cluster")]

annotations <- read.table("infercnv_annotations_steroidogenic.txt", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
colnames(annotations) <- c("CellID", "Condition")

cluster_annot <- merge(cell_clusters_sub, annotations, by = "CellID")

table(cluster_annot$Cluster, cluster_annot$Condition)

#What CNVs Were Detected in Each Cluster

#add back to seurat####
# Ensure the rownames match Seurat cell names
rownames(cluster_annot) <- cluster_annot$CellID

# Add inferCNV cluster to Seurat metadata
#seurat_combined$infercnv_cluster <- cluster_annot[colnames(seurat_combined), "Cluster"]
seurat_combined$infercnv_subclone <- cluster_annot$Cluster[match(colnames(seurat_combined), cluster_annot$CellID)]
DimPlot(seurat_combined, group.by = "infercnv_subclone", label = TRUE) 
DimPlot(seurat_combined, group.by = "simple_cell_type", label = TRUE) 

#Define "Malignant Cells"
seurat_combined$malignancy <- ifelse(is.na(seurat_combined$infercnv_subclone), "Normal-like", "Malignant")
DimPlot(seurat_combined, group.by = "malignancy", label = TRUE)

#To confirm subclones only exist in tumors:
table(seurat_combined$infercnv_subclone, seurat_combined$condition)

#Subset tumor cells for further analysis
tumor_cells <- colnames(seurat_combined)[seurat_combined$malignancy == "Malignant"]
tumor_seurat <- subset(seurat_combined, cells = tumor_cells)


FeaturePlot(seurat_combined, features = c("TOP2A", "MKI67", "CDK1", "SF1", "STAR"), split.by = "malignancy")

table(seurat_combined$infercnv_subclone)

# Marker genes for cortical layers
zG_markers <- c("DACH1", "VSNL1")
zF_markers <- c("CCN3", "NCAM1")
zR_markers <- c("CYB5A", "SULT2A1")
all_markers <- c(zG_markers, zF_markers, zR_markers)

# Average expression per cluster
avg_exp <- AverageExpression(
  seurat_combined,
  features = all_markers,
  group.by = "infercnv_cluster"
)$RNA

print(avg_exp)

#Differential Expression (DE) Between Tumor Clones####
Idents(tumor_seurat) <- tumor_seurat$infercnv_subclone
de_s1_vs_s2 <- FindMarkers(tumor_seurat, ident.1 = "csACT_s1", ident.2 = "csACT_s2", logfc.threshold = 0.25, min.pct = 0.1)
head(de_s1_vs_s2)
clone_ids <- unique(tumor_seurat$infercnv_subclone)
for (i in 1:(length(clone_ids)-1)) {
  for (j in (i+1):length(clone_ids)) {
    id1 <- clone_ids[i]
    id2 <- clone_ids[j]
    result <- FindMarkers(tumor_seurat, ident.1 = id1, ident.2 = id2)
    write.csv(result, paste0("DE_", id1, "_vs_", id2, ".csv"))
  }
}

#Annotate Inferred Subclones####
clone_markers <- FindAllMarkers(tumor_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(clone_markers, file = "clone_markers.csv")

top10 <- clone_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

topclustermarker <-DoHeatmap(tumor_seurat, features = top10$gene) + NoLegend()
ggsave("top_cluster_marker.png", topclustermarker, width = 10, height = 15, dpi = 300)

# Identifying cell types
DefaultAssay(tumor_seurat) <- "SCT"


cell_types <- list(
  "Progenitor"= c("NES","NR0B1","NR5A1"),
  "Proliferative" = c("MKI67", "TOP2A"),
  "Steroidogenic" = c("STAR", "ENSCAFG00000001285", "CYP17A1"),
  "Stress-responding" = c("JUN", "FOS","DUSP1")
)

simple_labels <- c("NR0B1" = "DAX1(NR0B1)",
                   "NR5A1" = "SF-1(NR5A1)",
                   "ENSCAFG00000001285" = "CYP11B2")

# Generate a bubble plot using Seurat's built-in function
bubble_annotation <- DotPlot(tumor_seurat, features = cell_types, group.by = "infercnv_subclone",dot.scale = 6) +
  scale_x_discrete(labels = simple_labels) +  # Apply custom labels
  theme_minimal() +
  ggtitle("Bubble Plot of Gene Expression Across Cluster")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0, 0, 0, 2), "cm"))+
  # Add custom gradient with all six colors
  scale_color_gradientn(colors = c("#4169E1", "#abd9e9", "#abd9e9", "#fee090", "#fee090", "#CD2626"),
                        name = "Expression Level")

ggsave("bubble_annotation.png", bubble_annotation, width = 8, height = 3, dpi = 300)

FeaturePlot(tumor_seurat, features = unlist(cell_types))

# Extract the data used in the plot
expr_table <- bubble_annotation$data
head(expr_table)
expr_table$gene <- as.character(expr_table$features.plot)
expr_table$gene_label <- simple_labels[expr_table$gene]
expr_table$gene_label[is.na(expr_table$gene_label)] <- expr_table$gene[is.na(expr_table$gene_label)]
print(expr_table)
# Save to CSV
write.csv(expr_table, "bubble_annotation_table.csv", row.names = FALSE)


#Monocle3 Pseudotime Across Tumor Clones####
library(monocle3)

# Use tumor_seurat from above
cds <- as.cell_data_set(tumor_seurat)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# Set root cells (optional but improves pseudotime)
# Use cells from less proliferative / more steroidogenic clone (e.g., s1) as root
root_cells <- colnames(tumor_seurat)[tumor_seurat$infercnv_subclone == "csACT_s1"]
cds <- order_cells(cds, root_cells = root_cells)
plot_cells(cds, color_cells_by = "pseudotime",cell_size = 2)
plot_cells(cds,
           color_cells_by = "infercnv_subclone",
           label_groups_by_cluster = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           cell_size = 2)


#Gene Expression Trends Along Pseudotime
deg_pseudo <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
pseudo_genes <- rownames(subset(deg_pseudo, q_value < 0.05))

deg_pseudo %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head(30)

plot_genes_in_pseudotime(cds[pseudo_genes[1:6], ])

#filter away ribsome gene
#Often light up in pseudotime because they scale with translational demand, cell state, stress, or growth
#May not reflect specific tumor biology, but rather global shifts (e.g., quiescent to active)
deg_filtered <- deg_pseudo %>%
  filter(status == 'OK') %>%
  arrange(q_value) %>%
  filter(!grepl("^RPS|^RPL|^UBA52|^FAU", rownames(.))) %>%
  head(50)

print(deg_filtered )

rowData(cds)$gene_short_name <- rownames(cds)
top_genes <- rownames(deg_filtered)[1:10]
plot_genes_in_pseudotime(cds[top_genes, ])

plot_genes_in_pseudotime(cds[c("CD79A", "DYRK1B", "PLCB1", "CACNA1C"), ])


cds <- preprocess_cds(cds, num_dim = 50)  # or whatever number of dimensions is appropriate
cds <- reduce_dimension(cds, reduction_method = "UMAP")
gene_module_df <- find_gene_modules(cds, resolution = 0.001)
plot_cells(cds, genes = gene_module_df, show_trajectory_graph = FALSE)

head(gene_module_df)
rowData(cds)$module <- gene_module_df$module[match(rownames(cds), gene_module_df$id)]

# Extract gene names for module 1
module1_genes <- gene_module_df %>%
  filter(module == 1) %>%
  pull(id)

# Now plot those genes
plot_cells(cds,
           genes = module1_genes,
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE)

#endo_as_ref####
# Set working directory and load your Seurat object
setwd("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/steroidogenic_cells/cnv")
seurat_combined <- readRDS("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/seurat_obj_annotated_updated.rds")

# 1. Ensure the correct identities are active
Idents(seurat_combined) <- "simple_cell_type"

# 2. Get the cell barcodes for each group
ref_cells <- WhichCells(seurat_combined, idents = "Endothelial", expression = condition == "Normal")
normal_steroid_cells <- WhichCells(seurat_combined, idents = c("Dedifferentiated_Steroidogenic", "Steroidogenic"), expression = condition == "Normal")
tumour_steroid_cells <- WhichCells(seurat_combined, idents = c("Dedifferentiated_Steroidogenic", "Steroidogenic"), expression = condition == "csACT") # Assuming 'csACT' is your tumour condition

# 3. Combine all cells for the analysis
analysis_cells <- c(ref_cells, normal_steroid_cells, tumour_steroid_cells)
infercnv_subset_new <- subset(seurat_combined, cells = analysis_cells)

# 4. Create a new annotation column with THREE distinct groups
# This is the corrected and most important step
infercnv_subset_new$infercnv_group <- "Tumour_Steroidogenic" # Default label
infercnv_subset_new$infercnv_group[colnames(infercnv_subset_new) %in% ref_cells] <- "Endothelial_Ref"
infercnv_subset_new$infercnv_group[colnames(infercnv_subset_new) %in% normal_steroid_cells] <- "Normal_Steroidogenic"


# Collapse the layers for the RNA assay in your new subset
infercnv_subset_new <- JoinLayers(infercnv_subset_new, assay = "RNA")
counts_matrix <- GetAssayData(infercnv_subset_new, assay = "RNA", layer = "counts")

# 5. Create the new annotations file
annotations <- data.frame(
  CellID = colnames(infercnv_subset_new),
  Group = infercnv_subset_new$infercnv_group
)

write.table(annotations,
            file = "infercnv_annotations_3groups.txt", # Use a new name
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. Create and run infercnv (using the new annotations file)
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix, 
  annotations_file = "infercnv_annotations_3groups.txt",
  gene_order_file = "filtered_gene_order.pos",
  ref_group_names = c("Endothelial_Ref") # Reference group name must match
)


# Run InferCNV
infercnv_obj_final <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "infercnv_output_endo_as_ref", # Use a new output directory
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE
)

#infercnv_obj_final <- readRDS("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/steroidogenic_cells/cnv/infercnv_output_endo_as_ref/run.final.infercnv_obj")
# The final matrix is in the 'expr.data' slot of the object
cnv_matrix <- infercnv_obj_final@expr.data

# Check the dimensions to make sure it looks correct (genes x cells)
dim(cnv_matrix)

# Calculate the CNV score for each cell
cnv_scores <- apply(cnv_matrix, 2, function(x) sum((x - 1)^2))


# Step 4: Create a data frame with cell barcodes and their scores
scores_df <- data.frame(
  cell_id = names(cnv_scores),
  cnv_score = cnv_scores
)

# Step 5: Get the group annotations for each cell from your Seurat object
# Assumes 'infercnv_subset_new' is the Seurat object you ran the analysis on
annotations_df <- data.frame(
  cell_id = colnames(infercnv_subset_new),
  group = infercnv_subset_new$infercnv_group # The 3-level annotation
)

# Step 6: Clean up barcodes and merge scores with annotations
scores_df$cell_id <- gsub("\\.", "-", scores_df$cell_id)
results_df <- merge(scores_df, annotations_df, by = "cell_id")

# Step 7: Create the violin plot
# Ensure the groups are ordered logically for the plot
results_df$group <- factor(results_df$group, levels = c("Endothelial_Ref", "Normal_Steroidogenic", "Tumour_Steroidogenic"))

cnv_plot <- ggplot(results_df, aes(x = group, y = cnv_score, fill = group)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_classic() +
  labs(
    title = "Comparison of CNV Scores",
    x = "Cell Group",
    y = "CNV Score (Sum of Squared Deviations)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  # Optional: Add statistical comparisons (e.g., Wilcoxon rank-sum test)
  stat_compare_means(comparisons = list(c("Endothelial_Ref", "Normal_Steroidogenic"),
                                        c("Normal_Steroidogenic", "Tumour_Steroidogenic"),
                                        c("Endothelial_Ref", "Tumour_Steroidogenic")),
                     method = "wilcox.test", label = "p.signif")

# Print the plot
print(cnv_plot)

# Save the plot
ggsave("cnv_score_violin_plot.png", plot = cnv_plot, width = 6, height = 6, dpi = 300)

# Plot the histogram to help you decide
hist(seurat_combined$cnv_score, breaks = 100)

# Load dplyr for the case_when function
library(dplyr)

# Set your threshold
cnv_threshold <- 15

# Create the new classification column directly in the object's metadata
seurat_combined$final_class <- case_when(
  seurat_combined$simple_cell_type %in% c("Steroidogenic", "Dedifferentiated_Steroidogenic") & seurat_combined$condition == "csACT" & seurat_combined$cnv_score > cnv_threshold ~ "Malignant",
  seurat_combined$simple_cell_type %in% c("Steroidogenic", "Dedifferentiated_Steroidogenic") & seurat_combined$condition == "Normal" & seurat_combined$cnv_score > cnv_threshold ~ "Pre-Malignant",
  seurat_combined$simple_cell_type %in% c("Steroidogenic", "Dedifferentiated_Steroidogenic") & seurat_combined$cnv_score <= cnv_threshold ~ "Diploid Steroidogenic",
  !seurat_combined$simple_cell_type %in% c("Steroidogenic", "Dedifferentiated_Steroidogenic") ~ "Diploid Stroma",
  TRUE ~ "Unclassified" # A catch-all
)

# Check the new classifications
table(seurat_combined$final_class)

# Set the active identity to your new classification
Idents(seurat_combined) <- "final_class"

# Plot the UMAP
DimPlot(seurat_combined, label = TRUE, repel = TRUE, pt.size = 0.5) +
  ggtitle("Final Cell Classification") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot the UMAP, split by the 'condition' column
DimPlot(seurat_combined, split.by = "condition", pt.size = 0.5)

# Plot the CNV score directly onto the UMAP
FeaturePlot(seurat_combined, features = "cnv_score", pt.size = 0.5) +
  scale_colour_viridis_c() +
  ggtitle("CNV Score Across Cells")

# Split by condition and colour by your final cell classification
DimPlot(seurat_combined, split.by = "condition", group.by = "final_class", pt.size = 0.5) +
  ggtitle("Cell Classification by Condition") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))








#human gene pos####
# 1. Extract Raw Counts Matrix
seurat_combined<- JoinLayers(seurat_combined, assay = "RNA")
counts_matrix <- GetAssayData(seurat_combined, assay = "RNA", slot = "counts")

# Verify dimensions
dim(counts_matrix)
colnames(counts_matrix) # Should match all cell barcodes

# 2. Create Annotations File
annotations <- data.frame(
  CellID = colnames(seurat_combined),
  Sample = seurat_combined$sample
)

write.table(annotations, file = "infercnv_annotations_steroidogenic.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 4. Create InferCNV Object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = "infercnv_annotations_steroidogenic.txt",
  gene_order_file = "infercnv_gene_order_human.txt",
  ref_group_names = c("NAD2016","NAD2017")
)

# 5. Run InferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # 10x Genomics data
  out_dir = "infercnv_output_humanref",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  noise_filter = 0.1,
  HMM = FALSE
)
