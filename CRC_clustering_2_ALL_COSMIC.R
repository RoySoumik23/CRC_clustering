#ALL_COSMIC

source("/home/shivalik/Downloads/R_codes/CRC_clustering_1.R")

# Define the file path for storing outputs
filepath <- "/home/shivalik2/Soumik/my_codes/Avik_bhaiya_project/all_COSMIC_variant/"

#======================================================================
# Log Normalisation DBS2-COSMIC Signatures
#======================================================================
# Log normalization + z-score (row-wise)
log_data <- log1p(df_COSMIC)  # log(1 + x) to handle zeros
df_scaled <- t(scale(t(log_data)))  # z-score per patient (row-wise)
df_scaled <- as.data.frame(df_scaled)

# Save the scaled dataset
write.csv(t(df_scaled), paste0(filepath, "df_ALL_COSMIC.csv"), row.names = T)

#======================================================================
# K-Means Clustering and Logging
#======================================================================

# Load the K-Means clustering script
source("/home/shivalik/Downloads/R_codes/CRC_clustering_3_Kmeans.R")

# Process K-Means clustering results
all_COSMIC_cluster <- kmean_clusters_df
all_COSMIC_cluster$PATIENT_IDs <- rownames(all_COSMIC_cluster)
colnames(all_COSMIC_cluster)[1] <- "Cluster"
rownames(all_COSMIC_cluster) <- 1:nrow(all_COSMIC_cluster)
all_COSMIC_cluster <- all_COSMIC_cluster[, c(2, 1)]  # Reorder columns for clarity

# Save the clustering results
write_csv(all_COSMIC_cluster, paste0(filepath, "kmean_clusters_all_COSMIC.csv"))

#======================================================================
# Kaplan-Meier (KM) Survival Plot
#======================================================================

# Generate the KM survival analysis plot
source("/home/shivalik/Downloads/R_codes/CRC_clustering_4_KMplot.R")
message("KM PLOT: PLOTTED AND SAVED")

#======================================================================
# Sankey Diagram (Cluster Comparison)
#======================================================================

# Load the Sankey diagram script to visualize clustering relationships
source("/home/shivalik/Downloads/R_codes/CRC_clustering_5_shanky.R")
message("Shanky diagram: PLOTTED AND SAVED")

#======================================================================
# Complex Heatmap (Cluster Comparison)
#======================================================================

df_htmp <- data.frame(fread(paste0(filepath, "df_ALL_COSMIC.csv"), sep=',', header=T), check.names = F, row.names=1)
meta_data_htmp <- data.frame(fread(paste0(filepath, "merged_cluster.csv"), sep = ",", header = TRUE))
source("/home/shivalik/Downloads/R_codes/ComplexHeatmap_CRC.R")
message("ComplexHeatmap: PLOTTED AND SAVED")
