#======================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
#======================================================================
# PART 3: K-Means Clustering & Heatmap Visualization
#======================================================================

# Set seed for reproducibility
set.seed(43)

# Perform K-means clustering with 4 clusters and 25 random starts for stability
kmeans_res <- kmeans(df_scaled, centers = num_clusters, nstart = 25)

# Convert clustering results into a dataframe
kmean_clusters_df <- as.data.frame(kmeans_res$cluster)

# Replace "-" with "." in row names (necessary processing to maintain the flow of the code)
rownames(kmean_clusters_df) <- gsub("\\.", "-", rownames(kmean_clusters_df))

# Prepare cluster annotations for visualization
annotations <- data.frame(Cluster = as.factor(kmean_clusters_df$`kmeans_res$cluster`))
rownames(annotations) <- rownames(df_scaled)

# Define a palette of colors
color_palette <- c(
  "#FF0000",  # Pure Red
  "#0000FF",  # Pure Blue
  "#008000",  # Pure Green
  "#FFA500",  # Bright Orange
  "#800080",  # Purple
  "#00FFFF",  # Cyan
  "#000000",  # Black
  "#FFFF00",  # Yellow
  "#A52A2A",  # Brown
  "#FFC0CB"   # Hot Pink
)

# Dynamically select only the number of colors needed
selected_colors <- setNames(color_palette[1:num_clusters], as.character(1:num_clusters))

# Use selected colors in annotation
annotation_colors <- list(Cluster = color_palette)


# Order the dataset so patients from the same cluster are grouped together
ordered_indices <- order(annotations$Cluster)
df_scaled <- df_scaled[ordered_indices, ]
annotations <- annotations[ordered_indices, , drop = FALSE]

# Save heatmap as a high-resolution PNG
png(filename = file.path(filepath, "heatmap_plot_clustered.png"), width = 20, height = 7, units = "in", res = 300)

# Define a custom color gradient for the heatmap
color_fun <- colorRamp2(
  breaks = c(min(df_scaled), 0, max(df_scaled)),  # Mapping data range to colors
  colors = c("#E817D4", "white", "#17E82B")       # Pink (low) → White (mid) → Green (high)
)

# Generate and save the heatmap
pheatmap(
  df_scaled,
  annotation_row = annotations,                                               # Add cluster annotations
  cluster_rows = F,                                                           # Keep rows in predefined cluster order
  cluster_cols = T,                                                           # Allow hierarchical clustering of columns
  annotation_colors = annotation_colors,                                      # Apply cluster-specific colors
  color = color_fun(seq(min(df_scaled), max(df_scaled), length.out = 100)),   # Apply custom color gradient
  show_rownames = F,                                                          # Hide row names for clarity
  show_colnames = T                                                           # Display column names
)

if (dev.cur() != 1) dev.off()
