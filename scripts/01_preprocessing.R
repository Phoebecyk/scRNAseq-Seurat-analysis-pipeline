# Load required libraries ####
library(Seurat)
library(sctransform)
library(tidyverse)
library(sleepwalk)
library(SCINA)
library(harmony)  # For integration across samples
library(patchwork)
library(ggplot2)
library(cowplot)
library(clusterProfiler)
library(org.Cf.eg.db)  # Canis familiaris (dog) organism database
library(ggplot2)
library(enrichplot)
library(dplyr)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(DOSE)
library(data.table)
library(org.Hs.eg.db)

# Set seed and theme
set.seed(123)
theme_set(theme_bw(base_size = 14))

# Setup relative paths
if (!require("here")) install.packages("here")
library(here)

# Use here() to automatically find the project root
base_path <- here()

#dev.off()  # Turn off current graphics device (repeat if multiple open)

# load human fetal sample
create_seurat_from_geo_xls <- function(file_path, project_name) {
  # Print a message to track progress
  cat("Processing:", file_path, "\n")
  
  # Use data.table's fread for fast and efficient file reading
  counts_matrix <- fread(file_path)
  gene_names <- counts_matrix[[1]]
  
  # Check for duplicated gene names
  if (any(duplicated(gene_names))) {
    cat("Warning: Duplicated gene names found in", project_name, ". Making them unique.\n")
    gene_names <- make.unique(gene_names)
  }
  
  # Convert the remaining columns to a numeric matrix
  # The `-1` removes the first column (gene names)
  counts_matrix <- as.matrix(counts_matrix[, -1])
  
  # Assign the gene names to the rows of the matrix
  rownames(counts_matrix) <- gene_names
  
  # Create the Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = project_name)
  
  cat("Successfully created Seurat object for", project_name, "with", ncol(seurat_obj), "cells.\n\n")
  
  return(seurat_obj)
}


# Specify the paths to your downloaded files ---
# Make sure these files are in your R working directory, or provide the full path.
file_f106 <- "adrenal_scRNA_seq_data/GSE137804_RAW/GSM4088787_F106_gene_cell_exprs_table.xls.gz"
file_f107 <- "adrenal_scRNA_seq_data/GSE137804_RAW/GSM4088788_F107_gene_cell_exprs_table.xls.gz"


# Create individual Seurat objects for each sample ---
# Call the function for each file.
f106_seurat <- create_seurat_from_geo_xls(file_path = file_f106, project_name = "f106")
f107_seurat <- create_seurat_from_geo_xls(file_path = file_f107, project_name = "f107")


convert_ensembl_to_symbol <- function(seurat_obj) {
  # 1. Extract Ensembl IDs
  ensembl_ids <- gsub("\\.\\d+$", "", rownames(seurat_obj))
  
  # 2. Map Ensembl IDs to Gene Symbols
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # 3. Handle Duplicates and Missing Values
  print(paste("Number of unmapped genes:", sum(is.na(gene_symbols))))
  unique_gene_symbols <- make.unique(ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols))
  
  # 4. Update the Seurat Object
  # UPDATED: Use 'layer' instead of the deprecated 'slot'
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts") 
  rownames(counts) <- unique_gene_symbols
  
  new_seurat_obj <- CreateSeuratObject(counts = counts, meta.data = seurat_obj@meta.data)
  
  # CORRECTED: Use standard R subsetting for features (genes)
  # The format is [rows_to_keep, columns_to_keep]
  genes_to_keep <- !grepl("^ENSG", rownames(new_seurat_obj))
  new_seurat_obj <- new_seurat_obj[genes_to_keep, ]
  
  return(new_seurat_obj)
}

# Convert Gene IDs ---
f106_seurat <- convert_ensembl_to_symbol(f106_seurat)
f107_seurat <- convert_ensembl_to_symbol(f107_seurat)

#load canine
canine_sample_paths <- list(
  "PCCVanDijk"     = "adrenal_scRNA_seq_data/PCCVanDijk/filtered_feature_bc_matrix/",
  "PCCDakota"      = "adrenal_scRNA_seq_data/PCCDakota/filtered_feature_bc_matrix/",
  "PCCRose"        = "adrenal_scRNA_seq_data/PCCRose/filtered_feature_bc_matrix/"
)

canine_seurat_list <- lapply(names(canine_sample_paths), function(sample_name) {
  counts <- Read10X(data.dir = canine_sample_paths[[sample_name]])
  CreateSeuratObject(counts = counts, project = sample_name)
})
names(canine_seurat_list) <- names(canine_sample_paths)

# --- 2. Combine all samples into a single list ---
seurat_list <- c(canine_seurat_list, list(f106 = f106_seurat, f107 = f107_seurat))

# --- Add Metadata (sample and condition) to all objects ---
seurat_list <- lapply(names(seurat_list), function(sample_name) {
  obj <- seurat_list[[sample_name]]
  obj$sample <- sample_name
  
  # Assign condition based on sample name
  obj$condition <- case_when(
    sample_name %in% c("f106", "f107") ~ "Normal",
    sample_name %in% c("PCCVanDijk", "PCCDakota", "PCCRose") ~ "PCC",
    TRUE ~ "Unknown"
  )
  return(obj)
})
names(seurat_list) <- names(c(canine_sample_paths, list(f106 = "f106", f107 = "f107")))