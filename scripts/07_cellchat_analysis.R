# CellChat v2 Analysis: Comparing Normal and Tumour (csACT) Conditions
# This script follows the protocol by Jin et al., Nature Protocols (2025)

library(sctransform)
library(CellChat)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(circlize) # This package provides the colorRamp2 function


# Set seed and theme
set.seed(123)
theme_set(theme_bw(base_size = 14))
setwd("")
#seurat_combined <- readRDS("seurat_obj_annotated_updated.rds")
#load everythings when new session start
cellchat_normal <- readRDS("cellchat_Normal.rds")
cellchat_normal <- netAnalysis_computeCentrality(cellchat_normal, slot.name = "netP")
cellchat_csACT <- readRDS("cellchat_csACT.rds")
cellchat_csACT <- netAnalysis_computeCentrality(cellchat_csACT, slot.name = "netP")
object.list <- list(Normal = cellchat_normal, csACT = cellchat_csACT)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# --- 2. Data Preparation ---####
seurat_combined<- JoinLayers(seurat_combined, assay = "RNA")
# Split the Seurat object by condition to create two separate objects.
seurat_normal <- subset(seurat_combined, subset = condition == "Normal")
seurat_csACT <- subset(seurat_combined, subset = condition == "csACT")

# --- 3. Create and Process CellChat Objects for Each Condition ---####
# This follows the initial steps of Procedure 1 for each dataset.

# -- For Normal Samples --
# Create CellChat object from Seurat object.
# group.by should be your cell type annotation column.

cellchat_normal <- createCellChat(object = seurat_normal, group.by = "simple_cell_type", assay = "RNA")

# Set the ligand-receptor interaction database.
# Use CellChatDB.human or CellChatDB.mouse depending on your species.
cellchat_normal@DB <- CellChatDB.human

# Pre-processing: Subset data and identify over-expressed genes/interactions.
cellchat_normal <- subsetData(cellchat_normal)
# Use multiple cores if available to speed up the process.
future::plan("multisession", workers = 4)
cellchat_normal <- identifyOverExpressedGenes(cellchat_normal)
cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)

# Infer cell-cell communication network.
cellchat_normal <- computeCommunProb(cellchat_normal, type = "triMean")
cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)
cellchat_normal <- computeCommunProbPathway(cellchat_normal)
cellchat_normal <- aggregateNet(cellchat_normal)
cellchat_normal <- netAnalysis_computeCentrality(cellchat_normal, slot.name = "netP")
# Save the processed Normal CellChat object
saveRDS(cellchat_normal, file = "cellchat_Normal.rds")

# -- For csACT (Tumour) Samples --
# Repeat the same process for the tumour condition.
cellchat_csACT <- createCellChat(object = seurat_csACT, group.by = "simple_cell_type", assay = "RNA")
cellchat_csACT@DB <- CellChatDB.human
cellchat_csACT <- subsetData(cellchat_csACT)
future::plan("multisession", workers = 4)
cellchat_csACT <- identifyOverExpressedGenes(cellchat_csACT)
cellchat_csACT <- identifyOverExpressedInteractions(cellchat_csACT)
cellchat_csACT <- computeCommunProb(cellchat_csACT, type = "triMean")
cellchat_csACT <- filterCommunication(cellchat_csACT, min.cells = 10)
cellchat_csACT <- computeCommunProbPathway(cellchat_csACT)
cellchat_csACT <- aggregateNet(cellchat_csACT)
# Save the processed csACT CellChat object
saveRDS(cellchat_csACT, file = "cellchat_csACT.rds")

# --- 4. Merge CellChat Objects for Comparative Analysis --- ####
# This is the core of Procedure 2.
object.list <- list(Normal = cellchat_normal, csACT = cellchat_csACT)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# Save the list of objects and the merged object
save(object.list, file = "cellchat_object.list_Normal_vs_csACT.RData")
save(cellchat_merged, file = "cellchat_merged_Normal_vs_csACT.RData")


# --- 5. Visualise and Compare Signalling Changes ---####

# -- Compare total number of interactions and interaction strength --####
# This gives a high-level overview of whether communication is generally
# increased or decreased in the tumour condition.
gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
# Plot side-by-side
p1 <- gg1 + gg2
print(p1)
ggsave(filename = "compareInteractions.png",  
       plot = p1,
       width = 6, height = 6,  dpi = 300)
# The plots will show the total number and strength of interactions for
# Normal (dataset 1) vs. csACT (dataset 2).

# -- Differential interaction visualisation --####
# Circle plot showing differential number of interactions.
# Red edges: increased signalling in csACT. Blue edges: decreased.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")

# Heatmap showing differential interactions between cell types.
# Red: increased in csACT. Blue: decreased in csACT.
gg1 <- netVisual_heatmap(cellchat_merged)
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
p2 <- gg1 + gg2
print(p2)

#Circle Plot showing interaction count in Normal and csACT
# Calculate the maximum interaction count across both datasets to ensure a consistent scale
weight.max <- getMaxWeight(object.list, slot.name = "net", attribute = "count")
# Set up a 1x2 plotting grid
par(mfrow = c(1,2), xpd=TRUE)
# Loop through each object in the list (Normal and csACT) to create a circle plot
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max, # Use the calculated max weight
    edge.width.max = 12,
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# 2D Scatter Plot: Major Sources and Targets
#(A) Identify cell populations with significant changes in sending or receiving signals
# First, compute the network centrality scores for each object in the list
# # This step is required before creating the scatter plot.
# for (i in 1:length(object.list)) {
#   object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
# }
# Calculate the total number of links for each cell type across both datasets
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
# Determine the min and max link counts to standardise dot sizes across plots
weight.MinMax <- c(min(num.link), max(num.link))
# diff scale
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
}
patchwork::wrap_plots(plots = gg)
# same scale
gg <- list()
for (i in 1:length(object.list)) {
  # Create the base plot first
  plot <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
  
  # Now, add the axis limits to the created ggplot object
  gg[[i]] <- plot + coord_cartesian(xlim = c(0, 20), ylim = c(0, 13))
}
patchwork::wrap_plots(plots = gg)

#(B) Identify the signaling changes of specific cell populations
# You can find your cell type names with: levels(cellchat_merged@idents$Normal)
# First, get a list of all unique cell type names from your dataset
cell_types <- levels(cellchat_merged@idents$Normal)

# Create an empty list to store the ggplot objects
plot_list <- list()

# Loop through each cell type
for (cell in cell_types) {
  # Generate the scatter plot for the current cell type
  plot_list[[cell]] <- netAnalysis_signalingChanges_scatter(
    cellchat_merged, 
    idents.use = cell
  )
}
# Arrange all the generated plots into a single figure using patchwork
print(patchwork::wrap_plots(plots = plot_list, ncol = 3))

# Plot for the first cell type of interest
gg1 <- netAnalysis_signalingChanges_scatter(
  cellchat_merged, # Use the merged object for this function
  idents.use = "Fibroblast", signaling.exclude = c("COLLAGEN","LAMININ")
)
print (gg1)
# Plot for the second cell type of interest
gg2 <- netAnalysis_signalingChanges_scatter(
  cellchat_merged,
  idents.use = "Endothelial", signaling.exclude = c("COLLAGEN","LAMININ")
)
print(gg2)

# Plot side-by-side
patchwork::wrap_plots(plots = list(gg1, gg2))


# -- Identify conserved and condition-specific signalling pathways --####
# This helps to pinpoint which pathways are most affected.
# RankNet visualises pathways that are turned on/off or have increased/decreased strength.
rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)


# -- Compare individual pathways or ligand-receptor pairs --
#Circle plot
# Example: Compare the CXCL pathway between conditions.
pathways.show <- c("CXCL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control for edge weights across plots
par(mfrow = c(1,2), xpd=TRUE)
for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

pathways.show <- c("COLLAGEN")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control for edge weights across plots
par(mfrow = c(1,2), xpd=TRUE)
for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

pathways.show <- c("CD6")
netVisual_aggregate(
  object.list[["csACT"]],
  signaling = pathways.show,
  layout = "circle",
  edge.width.max = 10,
  signaling.name = paste(pathways.show, "csACT")
)


# --- Create and Save Comparative Circle Plots for Key Pathways ---

# (Optional but recommended) Create a directory to save the plots
if (!dir.exists("circle_plots_keypathway")) {
  dir.create("circle_plots_keypathway")
}

# 1. Define the list of all pathways you want to visualise
pathways_to_plot <- c(
  # ECM & Adhesion Pathways
  "LAMININ", "COLLAGEN", "PTPRM", 
  # Angiogenesis Pathways
  "VEGF", "SEMA3", "ANGPT", 
  # Core Cancer Pathways
  "NOTCH", "FGF", "TGFb" 
)

# 2. Loop through the list to generate and save a plot for each pathway
for (pathway in pathways_to_plot) {
  
  # Open a PNG graphics device to save the file
  png(filename = paste0("circle_plots_keypathway/", pathway, "_comparison.png"), 
      width = 10, height = 5, units = "in", res = 300)
  
  # Calculate the max weight for consistent edge scaling across the pair of plots
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathway)
  
  # Set up a 1x2 plotting grid
  par(mfrow = c(1,2), xpd=TRUE)
  
  # Loop to draw the 'Normal' and 'csACT' plots side-by-side
  for (i in seq_along(object.list)) {
    netVisual_aggregate(
      object.list[[i]],
      signaling = pathway,
      layout = "circle",
      edge.weight.max = weight.max[1],
      edge.width.max = 10,
      signaling.name = paste(pathway, names(object.list)[i])
    )
  }
  
  # Close the graphics device to finalise saving the file
  dev.off()
}

# --- Create and Save DIFFERENCE Circle Plots for Key Pathways ---

# The directory should already be created by the previous step, but we check again
if (!dir.exists("circle_plots_keypathway_diff_circlize")) {
  dir.create("circle_plots_keypathway_diff_circlize")
}

# Load the 'circlize' package if it is not already loaded
# library(circlize) # Already in your initial script

# Define the list of all pathways you want to visualise
pathways_to_plot <- c(
  "LAMININ", "COLLAGEN", "PTPRM",
  "VEGF", "SEMA3", "ANGPT",
  "NOTCH", "FGF", "TGFb"
)

# Create a directory to save the differential plots
# Define a consistent colour scale for the differential plots
diff_col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Loop through the list to calculate the difference and plot
# Loop through the list to calculate the difference and plot
for (pathway in pathways_to_plot) {
  
  cat(paste("\nProcessing pathway:", pathway, "...\n"))
  
  # --- A. Matrix Extraction and Calculation (unchanged, as it works) ---
  mat1_lr <- subsetCommunication(object.list[["Normal"]], signaling = pathway, slot.name = "net")
  mat1 <- aggregate(mat1_lr$prob, by = list(source = mat1_lr$source, target = mat1_lr$target), FUN = sum)
  mat1 <- tapply(mat1$x, list(mat1$source, mat1$target), sum)
  mat1[is.na(mat1)] <- 0
  
  mat2_lr <- subsetCommunication(object.list[["csACT"]], signaling = pathway, slot.name = "net")
  mat2 <- aggregate(mat2_lr$prob, by = list(source = mat2_lr$source, target = mat2_lr$target), FUN = sum)
  mat2 <- tapply(mat2$x, list(mat2$source, mat2$target), sum)
  mat2[is.na(mat2)] <- 0
  
  all_cell_types <- levels(object.list[["Normal"]]@idents)
  mat1_full <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types), dimnames = list(all_cell_types, all_cell_types))
  mat2_full <- mat1_full
  mat1_full[rownames(mat1), colnames(mat1)] <- mat1
  mat2_full[rownames(mat2), colnames(mat2)] <- mat2
  
  mat_diff <- mat2_full - mat1_full
  
  idx <- which(rowSums(abs(mat_diff)) > 0 | colSums(abs(mat_diff)) > 0)
  
  if (length(idx) == 0) {
    cat("No *differential* interactions found (after L-R aggregation). Skipping plot.\n")
    next
  }
  
  mat_diff_filtered <- mat_diff[idx, idx, drop = FALSE]
  
  # --- B. Visualise the Difference using circlize::chordDiagram ---
  
  png(filename = paste0("circle_plots_keypathway_diff_circlize/", pathway, "_DIFFERENCE_Strength_Chord.png"),
      width = 10, height = 10, units = "in", res = 300)
  
  # INCREASE MARGINS: Provide more space at the top and bottom for the title and legend
  par(mfrow = c(1,1), mar = c(10, 5, 10, 5), xpd=FALSE) # Increased bottom (4) and top (4) margins
  
  # Create the chord diagram
  chordDiagram(
    mat_diff_filtered,
    directional = 1,
    transparency = 0.2,
    col = diff_col_fun,
    grid.col = "grey50",
    # Use no annotation tracks to keep the drawing area clean
    annotationTrack = "grid" # Only keep the sector grid line
  )
  
  # Manually plot the cell type names as a separate track (on the edge of the circle)
  # INCREASE LABEL OFFSET: Increased cm_h(3) to cm_h(6) to move text further out
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    cell.sector = get.cell.meta.data("sector.index")
    xcenter = get.cell.meta.data("xcenter")
    # Move the text outwards
    circos.text(xcenter, CELL_META$ylim[1] + cm_h(15), cell.sector,
                # CRITICAL CHANGE 1: Set to FALSE for horizontal text
                niceFacing = FALSE, 
                # CRITICAL CHANGE 2: Adjust justification for placement
                adj = c(1, 0), 
                cex = 0.8)
  }, bg.border = NA)
  
  # Move the TITLE to the top margin (using the increased margin space)
  mtext(paste("Differential Signalling Strength:", pathway, " (csACT - Normal)"), 
        side = 3, line = 2, outer = FALSE, font = 2) # line=2 uses space in the margin
  
  # Add the Legend (will be placed in the bottom margin space)
  lgd_links = Legend(
    at = c(-1, 0, 1),
    col_fun = diff_col_fun,
    title = "Interaction Strength Change",
    direction = "horizontal"
  )
  
  # Draw the legend (placed slightly higher in the bottom margin)
  draw(lgd_links, x = unit(0.5, "npc"), y = unit(0.2, "cm"), just = "center") # y=unit(0.2, "cm") places it near the bottom edge
  
  circos.clear()
  
  dev.off()
}


#heatmap
pathways.show <- c("COLLAGEN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

pathways.show <- c("VEGF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

pathways.show <- c("TGFb") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Ensure you are in the correct project directory first using setwd()

# Create a "heatmaps" sub-directory if it doesn't already exist
if (!dir.exists("heatmaps_keypathway")) {
  dir.create("heatmaps_keypathway")
}

# Define the list of pathways to visualise
pathways_to_plot <- c(
  "LAMININ", "COLLAGEN", "PTPRM", 
  "VEGF", "SEMA3", "ANGPT", 
  "NOTCH", "FGF", "TGFb" 
)

# Loop through the list to generate and save a plot for each pathway
for (pathway in pathways_to_plot) {
  
  # Generate the heatmap objects
  ht1 <- netVisual_heatmap(object.list[[1]], signaling = pathway, color.heatmap = "Reds", title.name = paste(pathway, names(object.list)[1]))
  ht2 <- netVisual_heatmap(object.list[[2]], signaling = pathway, color.heatmap = "Reds", title.name = paste(pathway, names(object.list)[2]))
  
  # Combine the heatmaps
  ht_list <- ht1 + ht2
  
  # Open a PNG device to save the file
  png(filename = paste0("heatmaps_keypathway/", pathway, "_comparison.png"), 
      width = 8, height = 6, units = "in", res = 300)
  
  # Draw the plot to the file
  ComplexHeatmap::draw(ht_list, ht_gap = unit(0.5, "cm"))
  
  # Close the device to finalise saving
  dev.off()
}


#Heatmap: Comparing Outgoing/Incoming Signalling Patterns
# Combine all identified signaling pathways from both datasets
pathway.union <- union(object.list[["Normal"]]@netP$pathways, object.list[["csACT"]]@netP$pathways)

# Generate the initial heatmap objects to extract their data 
# We won't plot these directly; we just need them as a source for the data matrices and annotations.
ht1_source <- netAnalysis_signalingRole_heatmap(
  object.list[["Normal"]],
  pattern = "outgoing",
  signaling = pathway.union
)

ht2_source <- netAnalysis_signalingRole_heatmap(
  object.list[["csACT"]],
  pattern = "outgoing",
  signaling = pathway.union
)

# Extract data and create a shared color scale
# Extract the raw data matrices from the source objects
mat1 <- ht1_source@matrix
mat2 <- ht2_source@matrix

# Find the maximum value across BOTH matrices to create a global scale
global_max <- max(c(mat1, mat2), na.rm = TRUE)

# Create a color mapping function that ranges from 0 to our global maximum
# This ensures the color scale is identical for both plots.
col_fun <- colorRamp2(c(0, global_max), c("white", "darkgreen")) # From light grey to dark red

# Create NEW Heatmap objects with the shared color scale
# We now build the heatmaps from their components, giving us full control.
ht1_final <- Heatmap(mat1,
                     name = "Relative strength", # This will be the legend title
                     col = col_fun,
                     cluster_rows = FALSE, # Keep CellChat's original ordering
                     cluster_columns = FALSE,
                     column_title = names(object.list)[1], # "Normal"
                     top_annotation = ht1_source@top_annotation,
                     right_annotation = ht1_source@right_annotation)

ht2_final <- Heatmap(mat2,
                     col = col_fun,
                     show_heatmap_legend = FALSE, # Hide the legend for the second plot
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     column_title = names(object.list)[2], # "csACT"
                     top_annotation = ht2_source@top_annotation,
                     right_annotation = ht2_source@right_annotation)

# Draw the final, synchronized heatmaps 
draw(ht1_final + ht2_final, ht_gap = unit(0.5, "cm"))

# -- Identify up-regulated and down-regulated signalling ligand-receptor pairs --####
# Use differential expression analysis to find significant changes.
# This identifies pairs where the ligand in the sender and/or the receptor in the receiver are significantly changed.
# --- Identify Up- and Down-regulated Signaling via a Single DE Analysis ---

# 1. Define the positive dataset for comparison (the tumour condition)
pos.dataset <- "csACT"
# Define a name for storing the differential expression results
features.name <- "csACT.DE"

# 2. Perform differential expression analysis ONCE
# This compares every cell type in 'csACT' vs 'Normal' and stores all results.
cellchat_merged <- identifyOverExpressedGenes(
  cellchat_merged, 
  group.dataset = "datasets", 
  pos.dataset = pos.dataset, 
  features.name = features.name, 
  only.pos = FALSE, 
  thresh.pc = 0.1, 
  thresh.fc = 0.05, # As recommended in the v2 tutorial for the faster test
  thresh.p = 0.05
)

# 3. Map all DE results (both up and down) onto the network
# The key here is `variable.all = TRUE`
net <- netMappingDEG(
  cellchat_merged, 
  features.name = features.name, 
  variable.all = TRUE
)

# 4. Subset the network to get UP-regulated interactions in csACT
# This selects interactions in the 'csACT' dataset where the ligand has a positive logFC.
net.up <- subsetCommunication(
  cellchat_merged, 
  net = net, 
  datasets = "csACT",
  ligand.logFC = 0.05, # Positive logFC threshold
  receptor.logFC = NULL # Can also filter by receptor if needed
)

# 5. Subset the network to get DOWN-regulated interactions in csACT
# This selects interactions in the 'Normal' dataset where the ligand has a negative logFC 
# (meaning it was more highly expressed in 'Normal' than 'csACT').
net.down <- subsetCommunication(
  cellchat_merged, 
  net = net, 
  datasets = "Normal",
  ligand.logFC = -0.05, # Negative logFC threshold
  receptor.logFC = NULL
)

# --- (Optional but Recommended) Extract the specific genes involved ---
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merged)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merged)

# --- Save the results to CSV files ---
write.csv(net.up, file = "net_up_csACT.csv", row.names = FALSE)
write.csv(net.down, file = "net_down_csACT.csv", row.names = FALSE)

# Visualise up-regulated signalling in csACT

netVisual_bubble(
  cellchat_merged,
  pairLR.use = net.up[, "interaction_name", drop = F],
  sources.use = NULL, # Check all sources
  targets.use = NULL, # Check all targets
  comparison = c(1, 2),
  angle.x = 45,
  remove.isolate = TRUE,
  title.name = "Up-regulated signaling in csACT"
)



# Visualise the down-regulated signalling in csACT.
# We use `max.dataset = 1` to highlight interactions where the first dataset ("Normal") has the higher strength.
netVisual_bubble(
  cellchat_merged,
  pairLR.use = net.down[, "interaction_name", drop = F],
  sources.use = NULL, # Check all sources
  targets.use = NULL, # Check all targets
  comparison = c(1, 2),
  max.dataset = 1, # This highlights interactions stronger in the first dataset
  angle.x = 45,
  remove.isolate = TRUE,
  title.name = "Down-regulated signaling in csACT"
)

# -- Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs-- #####
# For a single CellChat object
levels(cellchat_normal@idents)

#(A) Bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.up, sources.use = "Fibroblast", targets.use = "SmoothMuscle", comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.down, sources.use = "Fibroblast", targets.use = "SmoothMuscle", comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

netVisual_bubble(cellchat_merged, sources.use = 4, targets.use = c(5:6),  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat_merged, sources.use = 4, targets.use = c(5:6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in csACT", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = 4, targets.use = c(5:6),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in csACT", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

#(B) Chord diagram
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = 3, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = 3, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway

#prioritise your net.up and net.down results.
#1. Quantitative Ranking and Filtering
# Sort net.up by p-value, then by probability (strength)
net.up.sorted <- net.up[order(net.up$pval, -net.up$prob), ]

# Sort net.down similarly
net.down.sorted <- net.down[order(net.down$pval, -net.down$prob), ]

# View the top 10 most important up-regulated interactions
head(net.up.sorted, 10)

#2. Pathway-Level Aggregation
# Count occurrences of each pathway in your sorted lists
up_pathway_counts <- as.data.frame(table(net.up.sorted$pathway_name))
colnames(up_pathway_counts) <- c("Pathway", "NumberOfInteractions")

# Sort to see which pathways have the most up-regulated L-R pairs
up_pathway_counts <- up_pathway_counts[order(-up_pathway_counts$NumberOfInteractions), ]

# View the top pathways
head(up_pathway_counts,20)

#3. Biological Context: Sender-Receiver Analysis 
# Signals sent by tumour cells (source).
# Signals received by immune cells (target).
# Interactions between fibroblasts and endothelial cells.
# Example: Find up-regulated signals sent FROM Fibroblasts TO TCells
fibroblast_to_tcell_up <- subset(net.up.sorted, source == "Fibroblast" & target == "TCells")

# View the results
print(fibroblast_to_tcell_up)

#4. Visualisation of Top Hits
# Visualise only the top 20 up-regulated interactions
top_20_up <- head(net.up.sorted, 20)

netVisual_bubble(
  cellchat_merged,
  pairLR.use = top_20_up[, "interaction_name", drop = F],
  sources.use = NULL,
  targets.use = NULL,
  comparison = c(1, 2),
  angle.x = 45,
  remove.isolate = TRUE,
  title.name = "Top 20 Up-regulated Signals in csACT"
)

## --- Filter out major ECM pathways to see other signals ---

# Create a vector of the dominant pathways you want to temporarily exclude
pathways_to_exclude <- c("LAMININ", "COLLAGEN", "FN1", "THBS", "TENASCIN")

# Filter your pathway counts list
up_pathway_counts_filtered <- subset(up_pathway_counts, !(Pathway %in% pathways_to_exclude))

# View the next most abundant up-regulated pathways
head(up_pathway_counts_filtered, 20)


# --- Filter your main interaction list ---
net.up.sorted.filtered <- subset(net.up.sorted, !(pathway_name %in% pathways_to_exclude))

# View the top non-ECM interactions
head(net.up.sorted.filtered, 20)

# Bar Chart of Up-regulated Interaction Counts
# Select the top 20 pathways to display
top_pathways_df <- head(up_pathway_counts, 20)

# Create the plot
ggplot(data = top_pathways_df, aes(x = reorder(Pathway, NumberOfInteractions), y = NumberOfInteractions)) +
  geom_bar(stat = "identity", fill = "#b2182b") +
  coord_flip() + # Flips the axes to make labels readable
  labs(
    title = "Top Up-regulated Signalling Pathways in csACT",
    x = "Signalling Pathway",
    y = "Number of Significantly Up-regulated Interactions"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10)
  )

#--  extract all the ligand-receptor pairs for each of your selected pathways --####
# Access the main interaction database from CellChat
interaction_db <- CellChatDB.human$interaction

# Your list of key pathways
pathways_to_plot <- c(
  "LAMININ", "COLLAGEN", "PTPRM", 
  "VEGF", "SEMA3", "ANGPT", 
  "NOTCH", "FGF", "TGFb"
)

# Create an empty list to store the results
genes_by_pathway <- list()

# Loop through each pathway and extract its associated L-R pairs
for (pathway in pathways_to_plot) {
  genes_by_pathway[[pathway]] <- subset(interaction_db, pathway_name == pathway)
}

# Print the list to view all the genes, organised by pathway
print(genes_by_pathway)

# -- Check with DEG list -- ####
# 1. Define your list of genes and receptors of interest
gene_list <- c(
  "ACVR1_ACVR2A", "ACVR1_BMPR2", "ACVR1_TGFbR", "ACVR1B_TGFbR2", "ADM", "ADGRE5", 
  "ADGRL2", "ADGRL3", "APP", "AXL", "BMPR1A_ACVR2A", "BMPR1A_BMPR2", 
  "BMPR1B_ACVR2A", "BMPR1B_BMPR2", "CADM1", "CCL7", "CDH3", "CD44", "CD47", 
  "CD74", "CLDN1", "CLSTN1", "CLSTN2", "COL1A1", "COL3A1", "COL4A1", "COL4A2", 
  "COL6A1", "COL6A3", "DAG1", "EPHA2", "EPHA3", "FGFR1", "ICAM2", "IGF1R", 
  "IL34", "ITGA1_ITGB1", "ITGA5_ITGB1", "ITGA9_ITGB1", "ITGA11_ITGB1", 
  "ITGAV_ITGB1", "ITGAV_ITGB3", "ITGAV_ITGB5", "ITGAV_ITGB8", "JAG1", "JAM3", 
  "LAMA2", "LAMA4", "LAMB1", "LAMC1", "LAMC3", "MERTK", "MMP2", 
  "MMP14_ITGAV_ITGB3", "MPZL1", "NECTIN2", "NOTCH2", "NR1H4", "NRG1", "NRG3", 
  "NRP1_NRP2", "NRP1_PLXNA1", "NRP2_PLXNA1", "NRXN3", "NTRK3", "PCDHGC3", 
  "PDGFRA", "PDGFRB", "PGF", "PLAU", "PLXNA1", "PLXND1", "PTPRM", "PTPRS", 
  "ROBO1", "RORA", "SDC2", "SDC4", "SEMA5A", "SEMA6D", "TENM4", 
  "Testosterone-HSD17B12", "TGFB1", "TGFbR1_R2", "TGM2", "THBS1", "THBS2", 
  "UNC5B", "VCAM1", "VEGFA", "VEGFB", "VSIR", "VTN"
)

# 2. Read your CSV file into a data frame
# Make sure the file is in your current working directory, or provide the full path.
degs_fibroblast <- read.csv("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/deg_csv/degs_Fibroblast_normal_vs_csACT_filtered.csv")

# 3. Check which genes from your list are in the "gene" column of the data frame
genes_found <- gene_list[gene_list %in% degs_fibroblast$gene]
genes_not_found <- gene_list[!(gene_list %in% degs_fibroblast$gene)]

# 4. Display the results
cat("--- Check Complete ---\n\n")
cat(paste("Found", length(genes_found), "out of", length(gene_list), "genes in the file.\n\n"))

cat("Genes FOUND:\n")
print(genes_found)
# > print(genes_found)
# [1] "CD44"   "COL1A1" "COL4A1" "COL4A2" "COL6A1" "COL6A3" "EPHA3"  "LAMB1"  "LAMC1"  "MMP2"   "NR1H4"  "SDC2"   "SEMA5A" "TGFbR1_R2"  "TGM2"   "THBS1" 
# [16] "UNC5B"  "VEGFA" 

cat("\nGenes NOT FOUND:\n")
print(genes_not_found)

# Plot Fibroblast-Specific Pathways
# 1. Define the list of your differentially expressed genes of interest
genes_found <- c("CD44", "COL1A1", "COL4A1", "COL4A2", "COL6A1", "COL6A3", 
                 "EPHA3", "LAMB1", "LAMC1", "MMP2", "NR1H4", "SDC2", 
                 "SEMA5A","TGFbR1_R2", "TGM2", "THBS1", "UNC5B", "VEGFA")

# 2. Identify the unique pathways associated with these genes
# We search for these genes in both the ligand and receptor columns of the database.
interaction_db <- CellChatDB.human$interaction
pathways_of_interest <- unique(
  interaction_db[interaction_db$ligand %in% genes_found | 
                   grepl(paste(genes_found, collapse="|"), interaction_db$receptor), 'pathway_name']
)

# 3. Filter the up-regulated interactions (net.up) for these pathways AND involving Fibroblasts
net.up.fibroblast <- subset(net.up, 
                            (source == "Fibroblast" | target == "Fibroblast") & 
                              pathway_name %in% pathways_of_interest)

# 4. Filter the down-regulated interactions (net.down) similarly
net.down.fibroblast <- subset(net.down, 
                              (source == "Fibroblast" | target == "Fibroblast") & 
                                pathway_name %in% pathways_of_interest)

# 5. Create the bubble plots for visualisation

# Plot for UP-regulated interactions involving Fibroblasts
if (nrow(net.up.fibroblast) > 0) {
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.up.fibroblast[, "interaction_name", drop = F],
      sources.use = "Fibroblast",
      targets.use = NULL, # Set to NULL to see all targets
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Up-regulated Fibroblast SIGNALLING (Sender)"
    )
  )
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.up.fibroblast[, "interaction_name", drop = F],
      sources.use = NULL, # Set to NULL to see all sources
      targets.use = "Fibroblast",
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Up-regulated Fibroblast SIGNALLING (Receiver)"
    )
  )
} else {
  print("No up-regulated interactions found for the specified genes and Fibroblasts.")
}


# Plot for DOWN-regulated interactions involving Fibroblasts
if (nrow(net.down.fibroblast) > 0) {
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.down.fibroblast[, "interaction_name", drop = F],
      sources.use = "Fibroblast",
      targets.use = NULL,
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Down-regulated Fibroblast SIGNALLING (Sender)"
    )
  )
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.down.fibroblast[, "interaction_name", drop = F],
      sources.use = NULL,
      targets.use = "Fibroblast",
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Down-regulated Fibroblast SIGNALLING (Receiver)"
    )
  )
} else {
  print("No down-regulated interactions found for the specified genes and Fibroblasts.")
}

# --- Automated Gene Overlap Check ---

# 1. Define the cell type you want to analyse
target_cell_type <- "Endothelial"

# 2. Load the master list of all up-regulated interactions
net_up <- read.csv("net_up_csACT.csv")

# 3. Filter for interactions where the target cell type is the sender OR receiver
interactions_filtered <- subset(net_up, source == target_cell_type | target == target_cell_type)

# 4. Extract all unique gene names from the filtered interactions
# This handles both single ligands and complex multi-subunit receptors
ligands <- interactions_filtered$ligand
receptors_split <- unlist(strsplit(interactions_filtered$receptor, "_"))
endothelial_gene_list <- unique(c(ligands, receptors_split))

# 5. Load the corresponding DEG file for the target cell type
# The file path is constructed dynamically using the 'target_cell_type' variable
deg_file_path <- paste0("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/deg_csv/degs_", target_cell_type, "_normal_vs_csACT_filtered.csv")
degs_endothelial <- read.csv(deg_file_path)

# 6. Check which genes from the CellChat list are in the DEG list
genes_found <- endothelial_gene_list[endothelial_gene_list %in% degs_endothelial$gene]
genes_not_found <- endothelial_gene_list[!(endothelial_gene_list %in% degs_endothelial$gene)]

# 7. Display the results for the target cell type
cat(paste("--- Overlap Check for:", target_cell_type, "---\n\n"))
cat(paste("Found", length(genes_found), "out of", length(endothelial_gene_list), "genes from CellChat interactions in the DEG file.\n\n"))

cat("Genes FOUND:\n")
print(genes_found)

cat("\nGenes NOT FOUND:\n")
print(genes_not_found)

# Plot Endothelial-Specific Pathways

# 1. Define the list of your differentially expressed genes of interest for Endothelial cells
genes_found_endothelial <- c("PDGFD", "VEGFC", "PGF", "CCL7", "ANGPT2", "SEMA3A", "ADM", "SLIT2", 
                             "LAMA4", "LAMC1", "TNC", "VWF", "NECTIN2", "GJA1", "FLT4", "ITGB1", 
                             "KIT", "ROBO1", "ITGA8", "ITGA9", "CD74", "CD34")

# 2. Identify the unique pathways associated with these genes
interaction_db <- CellChatDB.human$interaction
pathways_of_interest_endo <- unique(
  interaction_db[interaction_db$ligand %in% genes_found_endothelial | 
                   grepl(paste(genes_found_endothelial, collapse="|"), interaction_db$receptor), 'pathway_name']
)

# 3. Filter the up-regulated interactions (net.up) for these pathways AND involving Endothelial cells
net.up.endothelial <- subset(net.up, 
                             (source == "Endothelial" | target == "Endothelial") & 
                               pathway_name %in% pathways_of_interest_endo)

# 4. Filter the down-regulated interactions (net.down) similarly
net.down.endothelial <- subset(net.down, 
                               (source == "Endothelial" | target == "Endothelial") & 
                                 pathway_name %in% pathways_of_interest_endo)

# 5. Create the bubble plots for visualisation

# Plot for UP-regulated interactions involving Endothelial cells
if (nrow(net.up.endothelial) > 0) {
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.up.endothelial[, "interaction_name", drop = F],
      sources.use = "Endothelial",
      targets.use = NULL, # Set to NULL to see all targets
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Up-regulated Endothelial SIGNALLING (Sender)"
    )
  )
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.up.endothelial[, "interaction_name", drop = F],
      sources.use = NULL, # Set to NULL to see all sources
      targets.use = "Endothelial",
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Up-regulated Endothelial SIGNALLING (Receiver)"
    )
  )
} else {
  print("No up-regulated interactions found for the specified genes and Endothelial cells.")
}


# Plot for DOWN-regulated interactions involving Endothelial cells
if (nrow(net.down.endothelial) > 0) {
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.down.endothelial[, "interaction_name", drop = F],
      sources.use = "Endothelial",
      targets.use = NULL,
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Down-regulated Endothelial SIGNALLING (Sender)"
    )
  )
  print(
    netVisual_bubble(
      cellchat_merged,
      pairLR.use = net.down.endothelial[, "interaction_name", drop = F],
      sources.use = NULL,
      targets.use = "Endothelial",
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "Down-regulated Endothelial SIGNALLING (Receiver)"
    )
  )
} else {
  print("No down-regulated interactions found for the specified genes and Endothelial cells.")
}

# --- Automated CellChat Analysis for Multiple Cell Types ---

# 1. Define the list of all cell types you want to analyse
cell_types_to_analyse <- c(
  "Dedifferentiated_Steroidogenic",
  "Endothelial",
  "Fibroblast",
  "Macrophage",
  "NeuralCrest",
  "SmoothMuscle",
  "Steroidogenic"
)

# Create a directory to save the plots if it doesn't exist
if (!dir.exists("cell_type_bubble_plots")) {
  dir.create("cell_type_bubble_plots")
}

# Load the master up- and down-regulated interaction lists (if not already in memory)
net_up <- read.csv("net_up_csACT.csv")
net_down <- read.csv("net_down_csACT.csv")

# Access the main interaction database
interaction_db <- CellChatDB.human$interaction

# 2. Start the main loop to process each cell type
for (cell_type in cell_types_to_analyse) {
  
  cat(paste("\n\n--- Processing Cell Type:", cell_type, "---\n"))
  
  # --- Part A: Find Overlapping DEGs ---
  
  # Filter for interactions where the cell type is the sender OR receiver
  interactions_filtered <- subset(net_up, source == cell_type | target == cell_type)
  
  # Extract all unique gene names from the filtered interactions
  ligands <- interactions_filtered$ligand
  receptors_split <- unlist(strsplit(interactions_filtered$receptor, "_"))
  cellchat_gene_list <- unique(c(ligands, receptors_split))
  
  # Construct the DEG file path dynamically
  deg_file_path <- paste0("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/N_csACT/deg_csv/degs_", cell_type, "_normal_vs_csACT_filtered.csv")
  
  # Load the corresponding DEG file
  degs_df <- read.csv(deg_file_path)
  
  # Find the genes present in both the CellChat interactions and the DEG list
  genes_found <- cellchat_gene_list[cellchat_gene_list %in% degs_df$gene]
  
  cat(paste("Found", length(genes_found), "overlapping DEG(s) for", cell_type, "\n"))
  
  # Print the actual overlapping gene names
  cat("Overlapping DEGs:\n")
  print(genes_found)
  
  # --- Part B: Plotting based on Overlapping Genes ---
  
  if (length(genes_found) == 0) {
    cat(paste("No overlapping DEGs found for", cell_type, ". Skipping plots.\n"))
    next # Skip to the next cell type in the loop
  }
  
  # Identify the unique pathways associated with the found genes
  pathways_of_interest <- unique(
    interaction_db[interaction_db$ligand %in% genes_found |
                     grepl(paste(genes_found, collapse="|"), interaction_db$receptor), 'pathway_name']
  )
  
  # Filter up-regulated interactions for these pathways AND the current cell type
  net.up.filtered <- subset(net.up,
                            (source == cell_type | target == cell_type) &
                              pathway_name %in% pathways_of_interest)
  
  # Filter down-regulated interactions similarly
  net.down.filtered <- subset(net.down,
                              (source == cell_type | target == cell_type) &
                                pathway_name %in% pathways_of_interest)
  
  # --- Visualisation and Saving ---
  
  # Plot and save for UP-regulated interactions
  if (nrow(net.up.filtered) > 0) {
    p1 <- netVisual_bubble(cellchat_merged, pairLR.use = net.up.filtered[, "interaction_name", drop = F], sources.use = cell_type, targets.use = NULL, comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE, title.name = paste("Up-regulated", cell_type, "SIGNALLING (Sender)"))
    p2 <- netVisual_bubble(cellchat_merged, pairLR.use = net.up.filtered[, "interaction_name", drop = F], sources.use = NULL, targets.use = cell_type, comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE, title.name = paste("Up-regulated", cell_type, "SIGNALLING (Receiver)"))
    
    # Save sender plot
    png(filename = paste0("cell_type_bubble_plots/Up-regulated_", cell_type, "_Sender.png"), width = 10, height = 20, units = "in", res = 300)
    print(p1)
    dev.off()
    
    # Save receiver plot
    png(filename = paste0("cell_type_bubble_plots/Up-regulated_", cell_type, "_Receiver.png"), width = 10, height = 20, units = "in", res = 300)
    print(p2)
    dev.off()
  }
  
  # Plot and save for DOWN-regulated interactions
  if (nrow(net.down.filtered) > 0) {
    p3 <- netVisual_bubble(cellchat_merged, pairLR.use = net.down.filtered[, "interaction_name", drop = F], sources.use = cell_type, targets.use = NULL, comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE, title.name = paste("Down-regulated", cell_type, "SIGNALLING (Sender)"))
    p4 <- netVisual_bubble(cellchat_merged, pairLR.use = net.down.filtered[, "interaction_name", drop = F], sources.use = NULL, targets.use = cell_type, comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE, title.name = paste("Down-regulated", cell_type, "SIGNALLING (Receiver)"))
    
    # Save sender plot
    png(filename = paste0("cell_type_bubble_plots/Down-regulated_", cell_type, "_Sender.png"), width = 10, height = 8, units = "in", res = 300)
    print(p3)
    dev.off()
    
    # Save receiver plot
    png(filename = paste0("cell_type_bubble_plots/Down-regulated_", cell_type, "_Receiver.png"), width = 10, height = 8, units = "in", res = 300)
    print(p4)
    dev.off()
  }
}

print("--- Automated analysis complete. All plots saved to 'cell_type_bubble_plots' folder. ---")
#TGF found in dediff and fibro DEGs as TFGBR3 and TGFBR2
