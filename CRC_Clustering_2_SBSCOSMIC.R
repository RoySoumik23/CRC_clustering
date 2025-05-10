#======================================================================
# Prognostic genome and transcriptome signatures in colorectal cancers
#======================================================================
# PART 2_SBSCOSMIC #
#======================================================================

# Load the necessary clustering results from Part 1
source("/home/shivalik/Downloads/R_codes/CRC_clustering_1.R")

# Define the file path for storing outputs
filepath <- "/home/shivalik2/Soumik/my_codes/Avik_bhaiya_project/sbs2_COSMIC_variant/"

#======================================================================
# Log Normalisation SBS2-COSMIC Signatures
#======================================================================
# Log normalization + z-score (row-wise)
log_data <- log1p(df_sbs2)  # log(1 + x) to handle zeros
df_scaled <- t(scale(t(log_data)))  # z-score per patient (row-wise)
df_scaled <- as.data.frame(df_scaled)

# Save the normalized dataset
write.csv(t(df_scaled), paste0(filepath, "df_scaled_sbs2_zscore.csv"), row.names = T)

#======================================================================
# K-Means Clustering and Logging
#======================================================================

# Load K-Means clustering script
source("/home/shivalik/Downloads/R_codes/CRC_clustering_3_Kmeans.R")

# Process K-Means clustering results
sbs2_cluster <- kmean_clusters_df
sbs2_cluster$PATIENT_IDs <- rownames(sbs2_cluster)
colnames(sbs2_cluster)[1] <- "Cluster"
rownames(sbs2_cluster) <- 1:nrow(sbs2_cluster)
sbs2_cluster <- sbs2_cluster[, c(2, 1)]  # Reorder columns for clarity

# Save the clustering results
write_csv(sbs2_cluster, file.path(filepath, "kmean_clusters_sbs2.csv"))

#======================================================================
# Kaplan-Meier (KM) Survival Plot
#======================================================================

# Load KM plot script and generate the survival analysis plot
source("/home/shivalik/Downloads/R_codes/CRC_clustering_4_KMplot.R")
message("KM PLOT: PLOTTED AND SAVED")

#======================================================================
# Sankey Diagram (Cluster Comparison)
#======================================================================

# Load script for visualizing the relationship between clustering results
source("/home/shivalik/Downloads/R_codes/CRC_clustering_5_shanky.R")
message("Shanky diagram: PLOTTED AND SAVED")

#======================================================================
# Complex Heatmap (Cluster Comparison)
#======================================================================

df_htmp <- as.data.frame(fread(paste0(filepath, "df_scaled_sbs2_zscore.csv"), sep = ",", header = TRUE))
meta_data_htmp <- as.data.frame(fread(paste0(filepath, "merged_cluster.csv"), sep = ",", header = TRUE))
source("/home/shivalik/Downloads/R_codes/ComplexHeatmap_CRC.R")
message("ComplexHeatmap: PLOTTED AND SAVED")
