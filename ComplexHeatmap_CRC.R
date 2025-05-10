#======================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
#======================================================================
# PART 6: Complex HeatMap generation
#======================================================================

# Filter to only samples present in df_htmp
meta_data_htmp <- meta_data_htmp[meta_data_htmp$Sample_IDs %in% colnames(df_htmp), ]
# Set rownames so we can align by sample ID
rownames(meta_data_htmp) <- meta_data_htmp$Sample_IDs
# Reorder rows of metadata to match column order in df_htmp
meta_data_htmp <- meta_data_htmp[colnames(df_htmp), ]
# Remove rows with NA
meta_data_htmp <- na.omit(meta_data_htmp)
# Keep only matching columns in df_htmp
df_htmp <- df_htmp[, rownames(meta_data_htmp)]
# Final check
stopifnot(all(colnames(df_htmp) == rownames(meta_data_htmp)))


df_htmp_numeric <- df_htmp[, sapply(df_htmp, is.numeric)]  # Keep only numeric columns

# Convert Cluster_kmean to factor and order data accordingly
meta_data_htmp$Cluster_kmean <- as.factor(meta_data_htmp$Cluster_kmean)
ord <- order(meta_data_htmp$Cluster_kmean)
df_htmp <- df_htmp[, ord]
meta_data_htmp <- meta_data_htmp[ord, ]
stopifnot(all(colnames(df_htmp) == rownames(meta_data_htmp)))  # Ensure alignment

df_htmp <- as.matrix(df_htmp)  # Convert data frame to matrix for heatmap

# Define heatmap annotation
# Only Clusters are used in the annotation
# Other annotations are commented out but can be added as needed

# Define color scale for heatmap values (ensuring range is correctly mapped)
min_val <- min(df_htmp, na.rm = TRUE)
max_val <- max(df_htmp, na.rm = TRUE)
mid_val <- 0  # Centering at zero
buffer <- (max_val - min_val) * 0.05  # Add a small buffer for better scaling

expr_col_fun <- colorRamp2(
  c(min_val - buffer, mid_val, max_val + buffer),  # Extend range slightly
  c("blue", "white", "red")  # Blue for low, white for mid, red for high values
)

# Define distinct colors for Clusters
cluster_levels <- levels(meta_data_htmp$Cluster_kmean)  # Get the unique cluster levels
cluster_colors <- setNames(c("#AF2FD0", "#2FA0D0", "#50D02F", "#D05F2F", "#D0AF2F", "#2FD0AF")[1:length(cluster_levels)], 
                           cluster_levels)  # Ensure correct mapping

# Define heatmap annotation with corrected color mapping
top_ha <- HeatmapAnnotation(
  Clusters = meta_data_htmp$Cluster_kmean,
  col = list(Clusters = cluster_colors),  # Assign correctly named color mapping
  annotation_name_side = "left",
  gap = unit(1, "mm")
)

# Create heatmap
ht <- Heatmap(
  df_htmp,
  name = "Expression",
  col = expr_col_fun,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split = meta_data_htmp$Cluster_kmean,
  top_annotation = top_ha,
  heatmap_legend_param = list(
    title = "Z-Score",
    at = c(round(min_val, 1), mid_val, round(max_val, 1)),  # Ensure zero is included
    legend_direction = "vertical"
  ),
  row_title = "Signature",
  column_title = "Samples"
)

# Save heatmap as JPEG
jpeg(file.path(filepath, "ComplexHeatmap.jpeg"), height = 13, width = 9, units = "in", res = 300)
draw(ht,
     merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     gap = unit(2, "mm"))
if (dev.cur() != 1) dev.off()
