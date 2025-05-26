#======================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
# PART 4: Kaplan-Meier Survival Analysis (Dynamic Cluster-Safe Version)
#======================================================================

# Prepare the dataframe for survival analysis
km_rfs_df <- kmean_clusters_df

# Add sample barcode as a column from rownames for merging
km_rfs_df$`DNA Tumor Sample Barcode` <- rownames(km_rfs_df)

# Rename first column to 'class' to represent cluster assignment
if (ncol(km_rfs_df) >= 2) {
  colnames(km_rfs_df)[1] <- "class"
  rownames(km_rfs_df) <- NULL
}

# Merge clustering info with clinical and survival data by sample barcode
km_rfs_df <- left_join(km_rfs_df, df_whole_data, by = "DNA Tumor Sample Barcode")

# Select columns relevant for overall survival analysis and covariates
km_rfs_df <- km_rfs_df[, c("class",
                           "Recurrence free survival days",
                           "Recurrence",
                           "DNA Tumor Sample Barcode",
                           "Age at diagnosis",
                           "Sex",
                           "Tumour Stage",
                           "Tumour Grade",
                           "MSI Status",
                           "CMS Tumour")]

# Ensure survival days is numeric and filter out missing values
km_rfs_df$`Recurrence free survival days` <- suppressWarnings(as.numeric(km_rfs_df$`Recurrence free survival days`))
km_rfs_df <- subset(km_rfs_df, !is.na(km_rfs_df$`Recurrence free survival days`) & !is.na(km_rfs_df$`Recurrence`))

# Convert survival days to months for easier clinical interpretation
km_rfs_df$`Recurrence free survival days` <- km_rfs_df$`Recurrence free survival days` / 30

# Encode vital status as binary: 1 = Dead (event), 0 = Alive (censored)
km_rfs_df$`Recurrence` <- ifelse(km_rfs_df$`Recurrence` == "Yes", 1, 0)

# Rename key columns for survival analysis function compatibility
colnames(km_rfs_df)[colnames(km_rfs_df) == "Recurrence free survival days"] <- "RFS_MONTHS"
colnames(km_rfs_df)[colnames(km_rfs_df) == "Recurrence"] <- "RFS_STATUS"
colnames(km_rfs_df)[colnames(km_rfs_df) == "DNA Tumor Sample Barcode"] <- "PATIENT_IDs"

# Create readable cluster labels dynamically
cluster_levels <- sort(unique(as.character(km_rfs_df$class)))
risk_labels <- setNames(paste0("CRC_cluster ", cluster_levels), cluster_levels)

# Assign factor levels and labels to cluster variable
km_rfs_df$class <- factor(km_rfs_df$class, levels = cluster_levels, labels = risk_labels)

# Convert clinical covariates to factors as needed
km_rfs_df$Sex <- factor(km_rfs_df$Sex)
km_rfs_df$`Tumour Stage` <- factor(km_rfs_df$`Tumour Stage`)
km_rfs_df$`Tumour Grade` <- factor(km_rfs_df$`Tumour Grade`)
km_rfs_df$`MSI Status` <- factor(km_rfs_df$`MSI Status`)
km_rfs_df$`CMS Tumour` <- factor(km_rfs_df$`CMS Tumour`)

# List of cluster pair combinations
cluster_pairs <- list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))

# Loop through each pair and generate KM plot
for (pair in cluster_pairs) {
  # Extract cluster labels
  cluster_a <- pair[1]
  cluster_b <- pair[2]
  
  # Subset the data for the two clusters
  pair_df <- subset(km_rfs_df, class %in% c(paste0("CRC_cluster ", cluster_a), paste0("CRC_cluster ", cluster_b)))
  
  # Drop unused factor levels
  pair_df$class <- droplevels(pair_df$class)
  
  # Fit KM curve
  kmcurve_pair <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ class, data = pair_df)
  
  # Fit Cox PH model for this pair
  coxph_model_pair <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ class + `Age at diagnosis` + Sex + 
                              `Tumour Stage` + `Tumour Grade` + `MSI Status` + `CMS Tumour`, data = pair_df)
  
  # Extract HR and 95% CI
  cox_summary <- summary(coxph_model_pair)
  hr <- round(cox_summary$coefficients[, "exp(coef)"], 2)
  hr_lower <- round(cox_summary$conf.int[, "lower .95"], 2)
  hr_upper <- round(cox_summary$conf.int[, "upper .95"], 2)
  hr_text <- paste0(names(hr), ": HR=", hr, " (", hr_lower, "-", hr_upper, ")")
  
  # File name for saving
  file_name <- paste0("KM_RFS_cluster_", cluster_a, "_vs_", cluster_b, ".tif")
  
  # Close any open graphics devices
  if (dev.cur() != 1) dev.off()
  
  # Save KM plot
  tiff(filename = file.path(filepath, file_name), width = 10, height = 10, units = "in", res = 300)
  
  # Generate KM plot
  g <- ggsurvplot(kmcurve_pair,
                  data = pair_df,
                  conf.int = FALSE,
                  pval = TRUE,
                  risk.table = "abs_pct",
                  risk.table.title = "Number at Risk",
                  risk.table.y.text.col = TRUE,
                  risk.table.y.text = FALSE,
                  surv.median.line = "hv",
                  legend.title = paste0("Clusters ", cluster_a, " vs ", cluster_b),
                  palette = color_palette,
                  title = paste0("Kaplan-Meier: CRC Cluster ", cluster_a, " vs ", cluster_b),
                  xlim = c(0, 60))
  
  # Customize x-axis
  g$plot <- g$plot + 
    ggplot2::scale_x_continuous(breaks = seq(0, 60, by = 10))
  
  # Annotate HR
  g$plot <- g$plot + 
    ggplot2::annotate("text", 
                      x = 0.8 * max(kmcurve_pair$time), 
                      y = 1, 
                      hjust = 1,
                      label = paste0("HR (1/ref):\n", 
                                     paste(names(hr), hr, sep = ": ", collapse = "\n")))
  
  # Print and save
  print(g)
  if (dev.cur() != 1) dev.off()
}