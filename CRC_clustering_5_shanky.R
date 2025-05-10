#=========================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
#=========================================================================
# PART 5: Sankey Diagram (Alluvial Plots)
#=========================================================================

#====================
# Data Preparation
#====================

# Extract relevant columns from the main dataset for previous clustering classifications
df_old_clustering <- df_whole_data[, c(2, 43:45)]
colnames(df_old_clustering)[1] <- "Sample_IDs"  # Rename the first column to "Sample_IDs"

# Prepare K-means clustering results
df_kmean <- kmean_clusters_df
df_kmean$Sample_IDs <- rownames(kmean_clusters_df)  # Assign row names as sample IDs
colnames(df_kmean)[1] <- "Cluster_kmean"  # Rename the first column for clarity

# Merge K-means clustering results with previous classifications
cluster_dfs_merged <- merge(df_kmean, df_old_clustering, by = "Sample_IDs")

# Export merged clustering dataset to CSV
write.csv(cluster_dfs_merged, paste0(filepath, "merged_cluster.csv"), row.names = FALSE)

# Standardize column names for consistency
colnames(cluster_dfs_merged)[2] <- "K_Mean Tumour"

# Clean classification labels by removing prefixes (e.g., "CMS", "iCMS", "CRPS")
cluster_dfs_merged$`CMS Tumour`  <- gsub("^CMS", "", cluster_dfs_merged$`CMS Tumour`)
cluster_dfs_merged$`iCMS Tumour` <- gsub("^iCMS", "", cluster_dfs_merged$`iCMS Tumour`)
cluster_dfs_merged$`CRPS Tumour` <- gsub("^CRPS", "", cluster_dfs_merged$`CRPS Tumour`)

# Replace "Undefined" values with 0 for uniformity
cluster_dfs_merged[cluster_dfs_merged == "Undefined"] <- "0"

#===============================
# Alluvial Plot: K-Means, CMS, iCMS, and CRPS
#===============================

# Group data and calculate frequency counts for visualization
df_shakny_4 <- cluster_dfs_merged %>% 
  group_by(`K_Mean Tumour`, `CMS Tumour`, `iCMS Tumour`, `CRPS Tumour`) %>% 
  summarize(Freq = n())

# Generate and save the alluvial plot for all four classifications
plot_path <- file.path(filepath, "alluvial_plot_all_four.jpeg")

jpeg(filename = paste0(filepath, "alluvial_plot_all_four.jpeg"), res = 300, height = 5, width = 10, units = "in")

p <- ggplot(data = df_shakny_4,
       aes(axis1 = `K_Mean Tumour`, axis2 = `CMS Tumour`, axis4 = `iCMS Tumour`, axis3 = `CRPS Tumour`, y = Freq)) +
  geom_alluvium(aes(fill = as.factor(`K_Mean Tumour`)), curve_type = "sigmoid") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("K_Mean Tumour", "CMS Tumour", "CRPS Tumour", "iCMS Tumour"), expand = c(0.15, 0.05)) +
  theme_minimal()

print(p)
if (dev.cur() != 1) dev.off()

#===============================
# Alluvial Plot: CMS vs. CRPS
#===============================

df_shakny_CMS_CRPS <- cluster_dfs_merged %>% 
  group_by(`CMS Tumour`, `CRPS Tumour`) %>% 
  summarize(Freq = n())

jpeg(filename = paste0(filepath, "alluvial_plot_cms_crps.jpeg"), res = 300, height = 5, width = 5, units = "in")

p <- ggplot(data = df_shakny_CMS_CRPS,
       aes(axis2 = `CMS Tumour`, axis1 = `CRPS Tumour`, y = Freq)) +
  geom_alluvium(aes(fill = as.factor(`CRPS Tumour`)), curve_type = "sigmoid") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CRPS Tumour", "CMS Tumour"), expand = c(0.15, 0.05)) +
  theme_minimal()

print(p)
if (dev.cur() != 1) dev.off()

#===============================
# Alluvial Plot: K-Means vs. CRPS
#===============================

df_shakny_kmean_crps <- cluster_dfs_merged %>% 
  group_by(`K_Mean Tumour`, `CRPS Tumour`) %>% 
  summarize(Freq = n())

jpeg(filename = paste0(filepath, "alluvial_plot_kmean_crps.jpeg"), res = 300, height = 5, width = 5, units = "in")

p <- ggplot(data = df_shakny_kmean_crps,
       aes(axis2 = `CRPS Tumour`, axis1 = `K_Mean Tumour`, y = Freq)) +
  geom_alluvium(aes(fill = as.factor(`K_Mean Tumour`)), curve_type = "sigmoid") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("K_Mean Tumour", "CRPS Tumour"), expand = c(0.15, 0.05)) +
  theme_minimal()

print(p)
if (dev.cur() != 1) dev.off()

#===============================
# Alluvial Plot: K-Means vs. CMS
#===============================

# Group data and count occurrences
df_shakny_kmean_cms <- cluster_dfs_merged %>% 
  group_by(`K_Mean Tumour`, `CMS Tumour`) %>% 
  summarize(Freq = n())

# Save alluvial plot as JPEG
jpeg(filename = paste0(filepath, "alluvial_plot_kmean_cms.jpeg"), 
     res = 300, height = 5, width = 5, units = "in")

# Generate alluvial plot while suppressing warnings
p <- ggplot(data = df_shakny_kmean_cms,
            aes(axis2 = `CMS Tumour`, axis1 = `K_Mean Tumour`, y = Freq)) +
  geom_alluvium(aes(fill = as.factor(`K_Mean Tumour`)), curve_type = "sigmoid") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("K_Mean Tumour", "CMS Tumour"),
                   expand = c(0.15, 0.05)) +
  theme_minimal()

print(p)
if (dev.cur() != 1) dev.off()
