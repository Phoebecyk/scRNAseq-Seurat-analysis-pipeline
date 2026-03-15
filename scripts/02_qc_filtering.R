# Define the list of human mitochondrial gene symbols (kept outside the loop for clarity, though it works inside)
human_mito_genes <- c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8")
human_mito_pattern <- paste0("^", human_mito_genes, "$", collapse = "|")

seurat_list <- lapply(seurat_list, function(obj) {
  
  # 1. Determine the correct mitochondrial gene pattern based on the condition
  mito_pattern <- if (obj$condition[1] == "Normal") {
    human_mito_pattern
  } else if (obj$condition[1] == "Canine_PCC") {
    "^MT-" # Canine pattern
  } else {
    "^MT-" # Default fallback
  }
  
  # 2. Find mitochondrial genes using the correct pattern for this object
  mito_genes_in_obj <- grep(
    pattern = mito_pattern, 
    x = rownames(obj), 
    value = TRUE, 
    ignore.case = TRUE
  )
  
  # 3. Robustly calculate percent.mt
  
  # A. Check for samples with no cells (0-cell objects will break PercentageFeatureSet)
  if (ncol(obj) == 0) {
    warning(paste0("Skipping percent.mt calculation for ", obj$sample[1], ": 0 cells."))
    obj[["percent.mt"]] <- 0 # Add a placeholder
    return(obj)
  }
  
  # B. Calculate percent.mt, handling errors (e.g., from empty feature list)
  obj[["percent.mt"]] <- tryCatch({
    # Suppress warnings about layers, which are common but not critical here
    suppressWarnings(
      PercentageFeatureSet(obj, features = mito_genes_in_obj)
    )
  }, error = function(e) {
    # If PercentageFeatureSet throws an error (e.g., if mito_genes_in_obj is empty 
    # in some Seurat versions/setups), set to zero.
    warning(paste0("PercentageFeatureSet error for ", obj$sample[1], ". Setting percent.mt to 0."))
    return(rep(0, ncol(obj)))
  })
  
  # C. CRITICAL STEP: Replace any remaining NaN (Not a Number) or NA with 0. 
  # This specifically addresses the 'missing value where TRUE/FALSE needed' error source.
  obj$percent.mt[is.nan(obj$percent.mt) | is.na(obj$percent.mt)] <- 0
  
  # 4. Also ensure nCount_RNA and nFeature_RNA are clean (if they contained NaN/NA)
  # This addresses the secondary error observed
  obj$nCount_RNA[is.nan(obj$nCount_RNA) | is.na(obj$nCount_RNA)] <- 0
  obj$nFeature_RNA[is.nan(obj$nFeature_RNA) | is.na(obj$nFeature_RNA)] <- 0
  
  return(obj)
})

# --- 3. Generate and collect QC plots (REVISION for Plotting Robustness) ---
plot_list <- lapply(names(seurat_list), function(sample) {
  obj <- seurat_list[[sample]]
  
  # ROBUSTNESS CHECK: Skip plotting for samples with 0 cells or constant metrics
  # A constant metric (e.g., all percent.mt = 0) can cause the VlnPlot/FeatureScatter error
  if (ncol(obj) < 2) {
    message(paste0("Skipping VlnPlot for ", sample, ": Too few cells or data is constant."))
    return(NULL)
  }
  
  # Check if any feature to be plotted is constant (i.e., min == max)
  # The 'all(data[, feature] == data[, feature][1])' check fails if feature is constant and plotting
  # a violin plot is meaningless in this case anyway.
  qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  if(any(sapply(obj@meta.data[, qc_features], function(x) min(x) == max(x)))){
    message(paste0("Skipping VlnPlot for ", sample, ": One or more QC features are constant."))
    return(NULL)
  }
  
  # Violin plot
  p1 <- VlnPlot(obj, features = qc_features) +
    ggtitle(paste0(sample, " - Violin Plot")) +
    theme(plot.title = element_text(size = 10))
  
  # ggplot QC scatter plot (nCount_RNA vs nFeature_RNA, colored by percent.mt)
  # Note: Using Seurat::as.data.frame to avoid 'as_tibble' dependency
  #qc.metrics <- Seurat::as.data.frame(obj[[]], row.names = "Cell.Barcode")
  # *** REVISED LINE: Use as_tibble (requires dplyr/tidyverse to be loaded) ***
  qc.metrics <- as_tibble(obj[[]], rownames = "Cell.Barcode") 
  
  # Use base R syntax for ggplot (avoids tidyverse dependency/loading issues)
  p3 <- ggplot2::ggplot(qc.metrics, ggplot2::aes(nCount_RNA, nFeature_RNA, colour = percent.mt)) + 
    ggplot2::geom_point(size = 0.5) + 
    ggplot2::scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
    ggplot2::ggtitle(paste0(sample, " - QC Metrics")) +
    ggplot2::geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
    ggplot2::geom_hline(yintercept = 6000, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
  
  # Combine the two plots for this sample into one row (requires 'patchwork')
  combined_plot <- p1 + p3 + patchwork::plot_layout(ncol = 4) # Changed ncol to 2
  
  return(combined_plot)
})

# Remove NULL elements (skipped plots) after the loop
plot_list <- plot_list[!sapply(plot_list, is.null)]

# Check QC metrics across all objects for NaN/NA issues
invisible(lapply(seurat_list, function(obj) {
  sample_name <- obj$sample[1]
  
  # Check 1: Ensure the object has cells
  if (ncol(obj) == 0) {
    cat(paste0("CRITICAL ERROR: Sample ", sample_name, " has 0 cells.\n"))
    return(NULL) # Skip to next object
  }
  
  # Check 2: Check 'percent.mt'
  mt_vector <- obj$percent.mt
  if (is.null(mt_vector)) {
    cat(paste0("ERROR: Sample ", sample_name, " is missing 'percent.mt' metadata.\n"))
  } else if (any(is.na(mt_vector)) || any(is.nan(mt_vector))) {
    # This should ideally be fixed by the prior code, but checking anyway
    cat(paste0("WARNING: Sample ", sample_name, " still contains NA/NaN in 'percent.mt'.\n"))
  } else if (all(mt_vector == mt_vector[1])) {
    # This condition directly mirrors the problematic internal VlnPlot check
    # It means all cells have the *exact* same percent.mt value (e.g., all are 0)
    cat(paste0("INFO: Sample ", sample_name, " has identical 'percent.mt' across all cells (Value: ", mt_vector[1], ").\n"))
  }
  
  # Check 3: Check 'nCount_RNA'
  ncount_vector <- obj$nCount_RNA
  if (all(ncount_vector == ncount_vector[1])) {
    cat(paste0("WARNING: Sample ", sample_name, " has identical 'nCount_RNA' across all cells. This indicates a potential issue with the raw data import.\n"))
  }
  
  # Check 4: Check 'nFeature_RNA'
  nfeature_vector <- obj$nFeature_RNA
  if (all(nfeature_vector == nfeature_vector[1])) {
    cat(paste0("WARNING: Sample ", sample_name, " has identical 'nFeature_RNA' across all cells. This indicates a potential issue with the raw data import.\n"))
  }
}))

# Combine all sample plots into one large figure (stack vertically)
all_plots <- wrap_plots(plot_list, ncol = 1)

#Save the combined plot to a file
ggsave(filename = "qc_plots_combined.png", 
       plot = all_plots, 
       width = 14,  # Width in inches (adjust as needed)
       height = 5 * length(seurat_list),  # Height scales with number of samples (4 inches per sample)
       dpi = 300,  # High resolution
       units = "in",
       limitsize = FALSE)

# After reviewing the saved plot, adjust thresholds here and filter
seurat_list <- lapply(seurat_list, function(obj) {
  # Initial suggested thresholds (adjust based on plot review)
  cells_to_filter <- rownames(subset(obj, subset = nFeature_RNA > 500 & 
                                       nFeature_RNA < 6000 & 
                                       percent.mt < 10)@meta.data)
  obj$keep <- rownames(obj@meta.data) %in% cells_to_filter
  obj <- subset(obj, subset = keep)
  return(obj)
})
rm(all_plots,plot_list)
##Check cell counts and variable features for each sample

cat("Sample Summary After Filtering:\n")
sample_summary <- lapply(names(seurat_list), function(sample) {
  n_cells <- ncol(seurat_list[[sample]])
  n_features <- nrow(seurat_list[[sample]])
  cat(sprintf("Sample %s: %d cells, %d features\n", sample, n_cells, n_features))
  return(data.frame(sample = sample, n_cells = n_cells, n_features = n_features))
})
sample_summary <- do.call(rbind, sample_summary)

# Remove samples with low cell counts if any (e.g., <100 cells) ####
min_cells <- 100
seurat_list <- seurat_list[sapply(seurat_list, ncol) >= min_cells]

# Define the maximum number of cells for downsampling
max_cells_per_sample <- 2000

cat("Downsampling all samples to a maximum of", max_cells_per_sample, "cells...\n")

# Use lapply to iterate through the list of Seurat objects
# This version iterates directly over the objects, which is slightly cleaner
seurat_list <- lapply(seurat_list, function(obj) {
  
  # Check if the current object has more cells than the maximum allowed
  if (ncol(obj) > max_cells_per_sample) {
    # Set a seed for reproducible random sampling
    set.seed(42)
    
    # Randomly sample the cell barcodes to keep
    cells_to_keep <- sample(colnames(obj), size = max_cells_per_sample)
    
    # Subset the object to keep only the sampled cells
    obj <- subset(obj, cells = cells_to_keep)
    
    # Print a confirmation message using the project name from the object
    cat(sprintf(" - %s downsampled to %d cells\n", obj@project.name, ncol(obj)))
    
  } else {
    # If the object is already at or below the threshold, do nothing
    cat(sprintf(" - %s kept with %d cells\n", obj@project.name, ncol(obj)))
  }
  
  # Return the (potentially modified) object
  return(obj)
})
names(seurat_list) <- names(seurat_list)  # Preserve sample names