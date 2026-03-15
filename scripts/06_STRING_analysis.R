#STRING (trying to use R to find bc and cc score: not the same)
# 1. Load Libraries
library(readr)
library(dplyr)
library(STRINGdb)
library(igraph)

# 2. Load Your DEG List
# Ensure the file path is correct
setwd("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/humanfetalN_PCC/STRING/neurochromaffin")
degs_df <- read_csv("/Users/phoebechan/Documents/adrenal_scRNA_seq_data/analysis/humanfetalN_PCC/deg_csv/no_1_0_degs_Neurosecretory_Chromaffin_Normal_vs_PCC.csv")
# Assuming your gene symbols are in a column named "gene"
degs_list <- degs_df$gene

#3. Connect to STRING Database and Map Genes
# We'll use human (species ID 9606) and a medium confidence score (400)
#string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400)
string_db <- STRINGdb$new(version="12", species=9615, score_threshold=700)
mapped_genes <- string_db$map(data.frame(gene = degs_list), "gene", removeUnmappedRows = TRUE)

# 4. Get the PPI Network as an igraph Object
ppi_graph <- string_db$get_graph()

# 5. Calculate Betweenness and Closeness Centrality
# Note: igraph uses STRING IDs as vertex names
bc_scores <- betweenness(ppi_graph, directed = FALSE)
cc_scores <- closeness(ppi_graph, mode = "all") # 'mode = "all"' for undirected graph

# 6. Create a Data Frame with Scores and Gene Names
# Match the scores (indexed by STRING ID) back to the mapped gene symbols
centrality_df <- data.frame(
  STRING_id = names(bc_scores),
  bc_score = bc_scores,
  cc_score = cc_scores
) %>%
  left_join(mapped_genes, by = c("STRING_id" = "STRING_id")) %>%
  dplyr::select(gene, bc_score, cc_score) %>% # Explicitly call dplyr's select function
  na.omit() # Remove any rows with missing data

# 7. Define Your Outlier Function and Identify Hub Genes
# This is the same function you provided
find_outliers <- function(scores) {
  Q1 <- quantile(scores, 0.25, na.rm = TRUE)
  Q3 <- quantile(scores, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  threshold <- Q3 + (1.5 * IQR)
  return(threshold)
}

Calculate thresholds
bc_threshold <- find_outliers(centrality_df$bc_score)
cc_threshold <- find_outliers(centrality_df$cc_score)

# Filter for hub genes
hub_genes <- centrality_df %>%
  filter(bc_score > bc_threshold & cc_score > cc_threshold)

# 8. Print and Save Results
print("Identified Hub Genes (Outliers in both BC and CC):")
print(hub_genes)

write.csv(hub_genes, "identified_hub_genes_Neurosecretory_Chromaffin.csv", row.names = FALSE)

# Save the full list of centrality scores to a CSV file
write.csv(centrality_df, "all_centrality_scores_Neurosecretory_Chromaffin.csv", row.names = FALSE)

Rank-Based Composite Score
1. Normalise scores using the 'centrality_df' object
# Min-Max normalisation scales the data to a [0, 1] range
centrality_df$bc_norm <- (centrality_df$bc_score - min(centrality_df$bc_score)) / (max(centrality_df$bc_score) - min(centrality_df$bc_score))
centrality_df$cc_norm <- (centrality_df$cc_score - min(centrality_df$cc_score)) / (max(centrality_df$cc_score) - min(centrality_df$cc_score))

# 2. Calculate a composite hub score (the average of the normalised scores)
centrality_df$hub_score <- (centrality_df$bc_norm + centrality_df$cc_norm) / 2

# 3. Identify hub genes by selecting the top 5% based on the composite score
top_n <- ceiling(0.05 * nrow(centrality_df))

hub_genes_composite <- centrality_df %>%
  arrange(desc(hub_score)) %>%
  head(top_n)

# 4. Print the results
print(paste("Identified Hub Genes (Top", top_n, "by Composite Score):"))
print(hub_genes_composite)

# 5. Save the identified hub genes to a CSV file
write.csv(hub_genes_composite, "ranked_hub_genes.csv", row.names = FALSE)

#Trying to mimic STRING/cytoscope way in calculating cc and bc score: not the same
# 1. Convert STRING confidence scores to distance weights
# This step is the same as before. Lower distance = higher confidence.
E(ppi_graph)$weight <- 1 - (E(ppi_graph)$combined_score / 1000)

# 2. Calculate both centrality scores using the edge weights
# The 'weights' argument tells igraph to find the "cheapest" paths, not just the shortest.
bc_scores_weighted <- betweenness(ppi_graph, directed = FALSE, weights = E(ppi_graph)$weight)
cc_scores_weighted <- closeness(ppi_graph, mode = "all", weights = E(ppi_graph)$weight)

# 3. Create a new, clean data frame with the weighted scores
centrality_df_weighted <- data.frame(
  STRING_id = names(bc_scores_weighted),
  bc_score = bc_scores_weighted,
  cc_score = cc_scores_weighted
) %>%
  left_join(mapped_genes, by = c("STRING_id" = "STRING_id")) %>%
  dplyr::select(gene, bc_score, cc_score) %>% # Use dplyr::select to avoid conflicts
  na.omit()

# 4. Display the top rows of the new data frame
print("Centrality scores calculated with weighted edges:")
head(centrality_df_weighted)

# 5. (Optional) Save the new weighted scores to a file
write.csv(centrality_df_weighted, "weighted_centrality_scores.csv", row.names = FALSE)