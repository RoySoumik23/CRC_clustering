#===============================================================================
# CRC Clustering - Script 1
#===============================================================================
# Input matrix: Mutation Signatures (rows) × Patients (columns)
df_t <- df_scaled  # Now: Patients × Features

# K-means-----------------------------------------------------------------------
set.seed(123)
kmeans_res <- kmeans(df_t, centers = 4, nstart = 25)
clusters_kmeans <- kmeans_res$cluster

# Hierarchical Clustering-------------------------------------------------------
dist_mat <- dist(df_t, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
clusters_hc <- cutree(hc, k = 4)

# Mclust------------------------------------------------------------------------
mclust_res <- Mclust(df_t, G = 1:6)
clusters_mclust <- mclust_res$classification

# Consensus Clustering----------------------------------------------------------
results_cc <- ConsensusClusterPlus(
  t(df_scaled),
  maxK = 6,
  reps = 500,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "pearson",
  seed = 123,
  plot = "png",  # or "pdf"
  title = "ConsensusPlots"
)

clusters_consensus <- results_cc[[4]]$consensusClass  # K = 4

#-------------------------------------------------------------------------------
# Combine cluster outputs into a list
#-------------------------------------------------------------------------------
cluster_outputs <- list(
  kmeans = clusters_kmeans,
  hierarchical = clusters_hc,
  mclust = clusters_mclust,
  consensus = clusters_consensus
)

# Save for next steps
saveRDS(cluster_outputs, "cluster_outputs.rds")

#===============================================================================
# KM Plot - Script 2
#===============================================================================
# Loop through all cluster outputs
for (method_name in names(cluster_outputs)) {
  
  cat("Running KM plot for:", method_name, "\n")
  
  # Assign current clustering result to km_df
  km_df <- as.data.frame(cluster_outputs[[method_name]])
  
  # Assign rownames as a new column for patient/sample identification
  km_df$pnts <- rownames(km_df)
  
  # Rename the first two columns for clarity and consistency
  if (ncol(km_df) >= 2) {
    colnames(km_df)[1] <- "class"
    colnames(km_df)[2] <- "DNA Tumor Sample Barcode"
  }
  
  # Reset row names to a numerical index
  rownames(km_df) <- seq_len(nrow(km_df))
  
  # Merge clustering data with full dataset based on sample barcode
  km_df <- left_join(km_df, df_whole_data, by = "DNA Tumor Sample Barcode")
  
  # Select relevant columns for survival analysis
  km_df <- km_df[, c(1, 31, 32, 2)]
  
  # Check and clean non-numeric values before conversion
  km_df$`Overall survival days` <- suppressWarnings(as.numeric(km_df$`Overall survival days`))
  
  # Filter out rows where either survival time or vital status is missing
  km_df <- subset(km_df, !is.na(km_df$`Overall survival days`) & !is.na(km_df$`Vital Status`))
  
  # Convert survival days to months for better interpretation
  km_df$`Overall survival days` <- km_df$`Overall survival days` / 30
  
  # Convert vital status to binary format (1 = Dead, 0 = Alive)
  km_df$`Vital Status` <- ifelse(km_df$`Vital Status` == "Dead", 1, 0)
  
  # Rename final columns for consistency with survival analysis functions
  colnames(km_df) <- c("class", "OS_MONTHS", "OS_STATUS", "PATIENT_IDs")
  
  # Dynamically generate risk labels
  cluster_levels <- sort(unique(as.character(km_df$class)))
  risk_labels <- setNames(paste0("CRC_cluster ", cluster_levels), cluster_levels)
  
  # Assign labels to the factor
  km_df$class <- factor(km_df$class, levels = cluster_levels, labels = risk_labels)
  
  # Fit Kaplan-Meier survival model
  kmcurve <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ class, data = km_df)
  
  # Fit Cox proportional hazards model
  coxph_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ class, data = km_df)
  hr <- summary(coxph_model)$coefficients[, "exp(coef)"]
  hr <- round(1 / hr, 2)  # Invert for consistent HR interpretation
  
  # Save KM plot
  if (dev.cur() != 1) dev.off()
  tiff(filename = file.path(filepath, paste0("survival_", method_name, ".tif")), 
       width = 8, height = 10, units = "in", res = 300)
  
  # Generate plot
  g <- ggsurvplot(kmcurve,
                  data = km_df,
                  conf.int = FALSE,
                  pval = TRUE,
                  risk.table = TRUE,
                  legend.labs = names(risk_labels),
                  legend.title = paste0("K-means Clusters (n = ", num_clusters, ")"),
                  title = "Kaplan-Meier Curve: CRC Clusters",
                  palette = color_palette,
                  xlim = c(0, 60)
  )
  
  # Modify x-axis breaks to intervals of 10 months
  g$plot <- g$plot + 
    ggplot2::scale_x_continuous(breaks = seq(0, 60, by = 10))
  
  # Print the plot
  print(g)
  
  # Annotate hazard ratios
  if (length(hr) > 0) {
    g$plot <- g$plot + 
      ggplot2::annotate("text", 
                        x = 0.8 * max(kmcurve$time), 
                        y = 1, 
                        hjust = 1,
                        label = paste0("HR (1/ref):\n", 
                                       paste(names(hr), hr, sep = ": ", collapse = "\n")))
  }
  
  # Print and save
  print(g)
  if (dev.cur() != 1) dev.off()
}