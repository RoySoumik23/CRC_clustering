#===============================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
# Final Part: Data Preprocessing & Preparation for final visualisation
#===============================================================================

#===============================================================================
# Load required libraries and set up the environment
#===============================================================================
source("https://raw.githubusercontent.com/RoySoumik23/R_scripts/refs/heads/main/libraries_and_setup.R")

#===============================================================================
# Load Data
#===============================================================================
df_immune <- read_excel("/home/shivalik2/Soumik/my_codes/Avik_bhaiya_project/PROGNOSTIC_Supplementary/Supplementary_Table_29.xlsx", skip = 3)

#===============================================================================
# Log Normalisation ALL-COSMIC Signatures
#===============================================================================
# Separate Patient ID column
patient_ids <- df_immune[, 1]
data_only <- df_immune[, -1]

# Log normalization and row-wise z-score
log_data <- log1p(data_only)  # log(1 + x)
normalized_data <- t(scale(t(log_data)))  # z-score row-wise

# Combine back with Patient ID
df_immune <- data.frame(Patient_ID = patient_ids, normalized_data, check.names = FALSE)

# df_whole_data ----
df_whole_data_new <- df_whole_data[, c(2,16,18,22,25,27,28,30,35,36,37,38,39,40)]
colnames(df_whole_data_new)[1] <- "Patient_ID"

#===============================================================================
# Working with signatures
#===============================================================================
df_signatures <- df_COSMIC
df_signatures <- data.frame(Patient_ID = rownames(df_signatures), df_signatures, row.names = NULL)
df_signatures$Total_mutation <- rowSums(df_signatures[, 2:ncol(df_signatures)], na.rm = TRUE)
df_signatures[, 2:(ncol(df_signatures)-1)] <- df_signatures[, 2:(ncol(df_signatures)-1)] / df_signatures$Total_mutation
df_signatures <- df_signatures[, !(names(df_signatures) %in% "Total_mutation")]

#===============================================================================
# Cluster
#===============================================================================
df_cluster <- kmean_clusters_df
df_cluster <- data.frame(Patient_ID = rownames(df_cluster), df_cluster, row.names = NULL)
colnames(df_cluster)[2] <- "Cluster"

#===============================================================================
# Merging df_whole_data_new & df_signatures
#===============================================================================
df_signatures[, 2:ncol(df_signatures)] <- lapply(df_signatures[, 2:ncol(df_signatures)], function(x) as.numeric(as.character(x)))

df_joined_sig_whole_data <- df_whole_data_new %>%
  left_join(df_cluster, by = "Patient_ID") %>%
  left_join(df_signatures, by = "Patient_ID")

df_joined_sig_whole_data[is.na(df_joined_sig_whole_data)] <- 0
write.csv(df_joined_sig_whole_data, file = paste0(filepath, "Signatures_Whole_data.csv"), quote = TRUE, row.names = FALSE)

#===============================================================================
# Fixing df_immune and Merging df_whole_data_new & df_immune
#===============================================================================
df_immune$Patient_ID <- gsub("\\.", "-", df_immune$Patient_ID)
df_immune[, 2:ncol(df_immune)] <- lapply(df_immune[, 2:ncol(df_immune)], function(x) as.numeric(as.character(x)))

df_joined_immune_whole_data <- df_whole_data_new %>%
  left_join(df_cluster, by = "Patient_ID") %>%
  left_join(df_immune, by = "Patient_ID")

df_joined_immune_whole_data[is.na(df_joined_immune_whole_data)] <- 0
write.csv(df_joined_immune_whole_data, file = paste0(filepath, "Immune_Whole_data.csv"), quote = TRUE, row.names = FALSE)



#===============================================================================
#===============================================================================
# VISUALISATION                                                                #
#===============================================================================
#===============================================================================


# Standardize column names and set Cluster as a categorical factor
names(df_joined_sig_whole_data) <- make.names(names(df_joined_sig_whole_data))
df_joined_sig_whole_data$Cluster <- factor(df_joined_sig_whole_data$Cluster)

#===============================================================================
# SHANKY PLOT
#===============================================================================

# 1. Cluster vs MSI Status----
jpeg(filename = paste0(filepath, "msi_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")

p <- ggplot(df_joined_sig_whole_data, 
            aes(axis1 = `MSI.Status`, axis2 = as.factor(`Cluster`))) +
  geom_alluvium(aes(fill = `MSI.Status`), width = 1/12, curve_type = "sigmoid") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("MSI.Status", "Cluster"), expand = c(.05, .05)) +
  labs(title = "Flow from MSI Status to Cluster", y = "Number of Patients") +
  theme_minimal()

print(p)
if (dev.cur() != 1) dev.off()
 
# 1.1 Bar plot of MSI Status vs clusters----
jpeg(filename = paste0(filepath, "MSI_Status_by_Cluster.jpeg"), 
     res = 300, height = 5, width = 5, units = "in")

p <- df_joined_sig_whole_data %>%
  group_by(Cluster, `MSI.Status`) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = Cluster, y = freq, fill = `MSI.Status`)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "MSI Status by Cluster", y = "Proportion") +
  scale_y_continuous(labels = scales::percent)

# Print plot
print(p)
if (dev.cur() != 1) dev.off()

# 2. Cluster vs Hypermutation status----
jpeg(filename = paste0(filepath, "Hypermutation_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")

p <- ggplot(df_joined_sig_whole_data, 
            aes(axis1 = `Hypermutation.status`, axis2 = as.factor(`Cluster`))) +
  geom_alluvium(aes(fill = `Hypermutation.status`), width = 1/15, curve_type = "sigmoid") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  scale_x_discrete(limits = c("Hypermutation.status", "Cluster"), expand = c(.05, .05)) +
  labs(title = "Flow from Hypermutation status to Cluster", y = "Number of Patients") +
  theme_minimal()

print(p)
if (dev.cur() != 1) dev.off()

#===============================================================================
# Violin Plot
#===============================================================================
# Convert Cluster to factor (if not already)
df_joined_sig_whole_data$Cluster <- as.factor(df_joined_sig_whole_data$Cluster)

# Calculate max y after conversion to MB
max_y <- max(df_joined_sig_whole_data$`Total.Mutation.Count`, na.rm = TRUE) / 1e6

# Define pairwise comparisons (adjust as needed)
comparisons <- list(
  c("1", "2"),
  c("1", "3"),
  c("1", "4"),
  c("2", "3"),
  c("2", "4"),
  c("3", "4")
)

jpeg(filename = paste0(filepath, "MutCount_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")

# Plot with y-axis in MB
p <- ggplot(df_joined_sig_whole_data, 
            aes(x = Cluster, y = `Total.Mutation.Count` / 1e6, fill = Cluster)) +  # scale here
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_minimal() +
  labs(title = "Mutation Burden Across Clusters",
       x = "Cluster",
       y = "Total Mutation Count (MB)") +  # updated label
  theme(legend.position = "none") +
  
  stat_compare_means(method = "kruskal.test", label.y = max_y * 1.05) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif")

# Print plot
print(p)
if (dev.cur() != 1) dev.off()

#===============================================================================
# Barplots for categorical phenotypes by Cluster
#===============================================================================

# Sex by Cluster----
pval <- chisq.test(table(df_joined_sig_whole_data$Sex, df_joined_sig_whole_data$Cluster))$p.value
p <- ggplot(df_joined_sig_whole_data, aes(x = Sex, fill = Cluster)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Sex", y = "Count", fill = "Cluster",
       subtitle = paste0("Chi-sq p = ", signif(pval, 3))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

jpeg(filename = paste0(filepath, "Sex_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Age.group by Cluster----
pval <- chisq.test(table(df_joined_sig_whole_data$Age.group, df_joined_sig_whole_data$Cluster))$p.value

p <- ggplot(df_joined_sig_whole_data, aes(x = Age.group, fill = Cluster)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Age Group", y = "Count", fill = "Cluster", subtitle = paste0("Chi-sq p = ", signif(pval, 3))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
jpeg(filename = paste0(filepath, "AgeGroup_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Bar Plot (Simple)
# Create a bar plot to show the distribution of Age.group across the clusters
p<- ggplot(df_joined_sig_whole_data, aes(x = Cluster, fill = Age.group)) +
  geom_bar(position = "fill") +  # position = "fill" makes the bars proportionate
  labs(title = "Age Group Distribution Across Clusters",
       x = "Cluster", y = "Proportion") +
  scale_fill_manual(values = c("≤ 65" = "skyblue", "66-79" = "orange", "≥ 80" = "purple")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate axis labels for better readability
jpeg(filename = paste0(filepath, "Age_percentage_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Anatomic.Organ.Subdivision by Cluster----
pval <- chisq.test(table(df_joined_sig_whole_data$Anatomic.Organ.Subdivision, df_joined_sig_whole_data$Cluster))$p.value

p <- ggplot(df_joined_sig_whole_data, aes(x = Anatomic.Organ.Subdivision, fill = Cluster)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Anatomic Organ Subdivision", y = "Count", fill = "Cluster", subtitle = paste0("Chi-sq p = ", signif(pval, 3))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
jpeg(filename = paste0(filepath, "OrganSubdivision_vs_cluster.jpeg"), res = 300, height = 5, width = 7, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Tumour.Grade by Cluster----
pval <- chisq.test(table(df_joined_sig_whole_data$Tumour.Grade, df_joined_sig_whole_data$Cluster))$p.value

p <- ggplot(df_joined_sig_whole_data, aes(x = Tumour.Grade, fill = Cluster)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Tumour Grade", y = "Count", fill = "Cluster", subtitle = paste0("Chi-sq p = ", signif(pval, 3))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
jpeg(filename = paste0(filepath, "TumourGrade_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Tumour.Stage by Cluster----
pval <- chisq.test(table(df_joined_sig_whole_data$Tumour.Stage, df_joined_sig_whole_data$Cluster))$p.value

p <- ggplot(df_joined_sig_whole_data, aes(x = Tumour.Stage, fill = Cluster)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Tumour Stage", y = "Count", fill = "Cluster", subtitle = paste0("Chi-sq p = ", signif(pval, 3))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
jpeg(filename = paste0(filepath, "TumourStage_vs_cluster.jpeg"), res = 300, height = 5, width = 6, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Hypermutation.status by Cluster----
pval <- chisq.test(table(df_joined_sig_whole_data$Hypermutation.status, df_joined_sig_whole_data$Cluster))$p.value

p <- ggplot(df_joined_sig_whole_data, aes(x = Hypermutation.status, fill = Cluster)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Hypermutation Status", y = "Count", fill = "Cluster", subtitle = paste0("Chi-sq p = ", signif(pval, 3))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
jpeg(filename = paste0(filepath, "HypermutationStatus_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# MSI.Status by Cluster----
pval <- chisq.test(table(df_joined_sig_whole_data$MSI.Status, df_joined_sig_whole_data$Cluster))$p.value

p <- ggplot(df_joined_sig_whole_data, aes(x = MSI.Status, fill = Cluster)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "MSI Status", y = "Count", fill = "Cluster", subtitle = paste0("Chi-sq p = ", signif(pval, 3))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
jpeg(filename = paste0(filepath, "MSIStatus_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

#===============================================================================
# Boxplots/Violin plots for continuous phenotypes by Cluster
#===============================================================================

# Total Mutation Count (Boxplot)----
pval <- kruskal.test(Total.Mutation.Count ~ Cluster, data = df_joined_sig_whole_data)$p.value
pval_label <- paste0("p = ", signif(pval, 3))

p <- ggplot(df_joined_sig_whole_data, aes(x = Cluster, y = Total.Mutation.Count, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cluster", y = "Total Mutation Count") +
  annotate("text", x = 1, y = max(df_joined_sig_whole_data$Total.Mutation.Count, na.rm = TRUE) * 1.05,
           label = pval_label, hjust = 0, size = 4) +
  theme_minimal() +
  theme(legend.position = "none")

jpeg(filename = paste0(filepath, "MutCount_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Mutated Driver Genes (Boxplot)----
pval <- kruskal.test(Mutated.Driver.Genes ~ Cluster, data = df_joined_sig_whole_data)$p.value
pval_label <- paste0("p = ", signif(pval, 3))

p <- ggplot(df_joined_sig_whole_data, aes(x = Cluster, y = Mutated.Driver.Genes, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cluster", y = "Mutated Driver Genes") +
  annotate("text", 
           x = 1, 
           y = max(df_joined_sig_whole_data$Mutated.Driver.Genes, na.rm = TRUE) * 1.05,
           label = pval_label, 
           hjust = 0, size = 4) +
  theme_minimal() +
  theme(legend.position = "none")

jpeg(filename = paste0(filepath, "MutatedDriverGenes_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Copy Number Segments (Violin)----
pval <- kruskal.test(Copy.Number.Segments ~ Cluster, data = df_joined_sig_whole_data)$p.value
pval_label <- paste0("p = ", signif(pval, 3))

p <- ggplot(df_joined_sig_whole_data, aes(x = Cluster, y = Copy.Number.Segments, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cluster", y = "Copy Number Segments") +
  annotate("text", 
           x = 1, 
           y = max(df_joined_sig_whole_data$Copy.Number.Segments, na.rm = TRUE) * 1.05, 
           label = pval_label, 
           hjust = 0, size = 4) +
  theme_minimal() +
  theme(legend.position = "none")

jpeg(filename = paste0(filepath, "CopyNumber_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

# Structural Variants (Boxplot)----
pval <- kruskal.test(Structural.Variants ~ Cluster, data = df_joined_sig_whole_data)$p.value
pval_text <- paste0("Kruskal-Wallis p = ", signif(pval, 3))

p <- ggplot(df_joined_sig_whole_data, aes(x = Cluster, y = Structural.Variants, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cluster", y = "Structural Variants") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = Inf, y = Inf, label = pval_text, hjust = 1.1, vjust = 1.5, size = 4)

jpeg(filename = paste0(filepath, "StructuralVariants_vs_cluster.jpeg"), res = 300, height = 5, width = 5, units = "in")
print(p)
if (dev.cur() != 1) dev.off()

#===============================================================================
# Ridge Plot
#===============================================================================
pval <- kruskal.test(Total.Mutation.Count ~ Cluster, data = df_joined_sig_whole_data)$p.value
pval_text <- paste0("Kruskal-Wallis p = ", signif(pval, 3))

p <- ggplot(df_joined_sig_whole_data, aes(x = Total.Mutation.Count, y = Cluster, fill = Cluster)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Total Mutation Count", y = "Cluster", subtitle = pval_text) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10, hjust = 0)
  ) +
  # Add as plot annotation in case subtitle doesn't render in saved image
  annotate("text", x = Inf, y = -Inf, label = pval_text, hjust = 1.1, vjust = -1.5, size = 4)

jpeg(filename = paste0(filepath, "Ridge_MutCount_vs_cluster.jpeg"), res = 300, height = 5, width = 6, units = "in")
print(p)
if (dev.cur() != 1) dev.off()


#===============================================================================
# t-sne plot
#===============================================================================

df_scaled_sigs <- df_COSMIC
df_clusters_kmeans <- kmean_clusters_df

df_scaled_sigs <- rownames_to_column(df_scaled_sigs, var = "Patient_ID") #tibble
df_clusters_kmeans <- rownames_to_column(df_clusters_kmeans, var = "Patient_ID") #tibble
df_t_sne <- left_join(df_scaled_sigs, df_clusters_kmeans, by = "Patient_ID")
colnames(df_t_sne)[colnames(df_t_sne) == "kmeans_res$cluster"] <- "Cluster"

# Plot
# Remove non-numeric columns (like IDs)
numeric_data <- df_t_sne[, -1]  # Exclude first column (ID)

# Run t-SNE
set.seed(42)  # for reproducibility
tsne_result <- Rtsne(as.matrix(numeric_data), dims = 2, perplexity = 2, verbose = TRUE)

# Create a data frame for ggplot
tsne_df <- data.frame(
  X = tsne_result$Y[,1],
  Y = tsne_result$Y[,2],
  Label = df_t_sne$Patient_ID,
  cluster = factor(df_t_sne$Cluster)
)

# Define custom colors for each cluster
custom_colors <- c(
  "1" = "#AF2FD0",
  "2" = "#2FA0D0",
  "3" = "#50D02F",
  "4" = "#D05F2F"
)

# Plot
p <- ggplot(tsne_df, aes(x = X, y = Y, color = cluster)) +  # map color to cluster variable
  geom_point(size = 0.7) +
  ggtitle("t-SNE Plot") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(color = "Cluster") +
  scale_color_manual(values = custom_colors)

jpeg(filename = paste0(filepath, "t-sne plot.jpeg"), res = 300, height = 5, width = 6, units = "in")

print(p)
if (dev.cur() != 1) dev.off()



#===============================================================================
# Immune based clustering
#===============================================================================
df_immune_clustered <- left_join(df_immune, df_cluster, by = "Patient_ID")
df_immune_clustered <- df_immune_clustered[!is.na(df_immune_clustered$Cluster), ]
df_immune_clustered <- df_immune_clustered[, -c(24:27, 29:30, 43, 46, 52:58, 63, 64, 66, 75, 77, 81, 83, 85:88, 92:94)]

# Sort the data frame by cluster before transposing:
df_immune_clustered <- df_immune_clustered[order(df[[ncol(df)]]), ]

# Prepare data
heat_data <- df_immune_clustered[, 2:(ncol(df_immune_clustered) - 1)]  # Select only the heatmap data
rownames(heat_data) <- df_immune_clustered$patient_id  # Set patient IDs as rownames

# Transpose so that features become rows (Y-axis) and patients are columns (X-axis)
heat_matrix <- t(as.matrix(heat_data))

# Extract cluster information for column annotation
clusters <- df_immune_clustered[[ncol(df_immune_clustered)]]
names(clusters) <- df_immune_clustered$patient_id

# Define custom colors
custom_colors <- c(
  "1" = "#AF2FD0",
  "2" = "#2FA0D0",
  "3" = "#50D02F",
  "4" = "#D05F2F"
)

# Create column annotation
ha <- HeatmapAnnotation(
  cluster = as.factor(clusters),
  col = list(cluster = custom_colors),
  annotation_name_side = "left"
)

# Generate heatmap
p <- Heatmap(
  heat_matrix,
  name = "Expression",
  top_annotation = ha,
  column_split = as.factor(clusters),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE
)

jpeg(filename = paste0(filepath, "Immune_heatmap.jpeg"), res = 300, height = 25, width = 16, units = "in")

print(p)
if (dev.cur() != 1) dev.off()


#===============================================================================
# Driver Gene mutation box plot
#===============================================================================

# Read input data from CSV file
data <- read.csv("Immune_Whole_data.csv")

# Ensure 'Cluster' is treated as a categorical variable
data$Cluster <- as.factor(data$Cluster)

# Define custom color palette for each cluster
custom_colors <- c(
  "1" = "#AF2FD0",
  "2" = "#2FA0D0",
  "3" = "#50D02F",
  "4" = "#D05F2F"
)

# Create boxplot with overlaid jittered points and statistical comparisons
p <- ggplot(data, aes(x = Cluster, y = Mutated.Driver.Genes, fill = Cluster)) +
  
  # Draw boxplots with reduced edge thickness and transparent fill
  geom_boxplot(size = 0.3,
               alpha = 0.7,
               outlier.shape = NA) +
  
  # Add individual data points with jitter for visibility
  geom_jitter(
    aes(color = Cluster),
    position = position_jitter(width = 0),
    size = 0.5,
    alpha = 0.6
  ) +
  
  # Apply custom colors to both fill and jitter points
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
  # Set plot labels
  labs(title = "Driver Gene Mutations per Cluster", x = "Cluster", y = "Number of Driver Genes Mutated") +
  
  # Adjust overall theme and styling
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 3),
    legend.position = "right"
  ) +
  
  # Perform and display pairwise t-tests with p-values between clusters
  stat_compare_means(
    comparisons = list(
      c("1", "2"),
      c("1", "3"),
      c("1", "4"),
      c("2", "3"),
      c("2", "4"),
      c("3", "4")
    ),
    method = "t.test",
    label = "p",
    size = 2,
    step.increase = 0.1,
    hide.ns = TRUE
  )

# Save the resulting plot as a high-resolution JPEG file
jpeg(filename = paste0(filepath, "DriverGeneMutations_Boxplot_Pvalues.jpeg"), res = 300, height = 25, width = 16, units = "in")

print(p)
if (dev.cur() != 1) dev.off()
