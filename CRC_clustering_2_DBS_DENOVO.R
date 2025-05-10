#======================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
#======================================================================
# PART 2_DBS78
#======================================================================

# Load the necessary clustering results from Part 1
source("/home/shivalik/Downloads/R_codes/CRC_clustering_1.R")

# Define the file path for storing outputs
filepath <- "/home/shivalik2/Soumik/my_codes/Avik_bhaiya_project/dbs1_denovo_variant/"

#======================================================================
# Log Normalisation DBS1-COSMIC Signatures
#======================================================================
# Log normalization + z-score (row-wise)
log_data <- log1p(df_dbs1)  # log(1 + x) to handle zeros
df_scaled <- t(scale(t(log_data)))  # z-score per patient (row-wise)
df_scaled <- as.data.frame(df_scaled)

# Save the processed data
write.csv(t(df_scaled), paste0(filepath, "df_dbs_denovo.csv"), row.names = T)

#======================================================================
# K-Means Clustering and Logging
#======================================================================

# Load the K-Means clustering script
source("/home/shivalik/Downloads/R_codes/CRC_clustering_3_Kmeans.R")

# Process K-Means clustering results
dbs1_cluster <- kmean_clusters_df
dbs1_cluster$PATIENT_IDs <- rownames(dbs1_cluster)
colnames(dbs1_cluster)[1] <- "Cluster"
rownames(dbs1_cluster) <- 1:nrow(dbs1_cluster)
dbs1_cluster <- dbs1_cluster[, c(2, 1)]  # Reorder columns for clarity

# Save the clustering results
write_csv(dbs1_cluster, paste0(filepath, "kmean_clusters_dbs1.csv"))

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

df_htmp <- as.data.frame(fread(paste0(filepath, "df_dbs_denovo.csv"), sep = ",", header = TRUE))
meta_data_htmp <- as.data.frame(fread(paste0(filepath, "merged_cluster.csv"), sep = ",", header = TRUE))
source("/home/shivalik/Downloads/R_codes/ComplexHeatmap_CRC.R")
message("ComplexHeatmap: PLOTTED AND SAVED")
