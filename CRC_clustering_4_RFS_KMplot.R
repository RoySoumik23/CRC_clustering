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

# Fit Kaplan-Meier survival model stratified by cluster
kmcurve <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ class, data = km_rfs_df)

# Fit multivariate Cox proportional hazards model adjusting for clinical covariates
coxph_model <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ class + `Age at diagnosis` + Sex + 
                       `Tumour Stage` + `Tumour Grade` + `MSI Status` + `CMS Tumour`, data = km_rfs_df)

# Extract hazard ratios (HR) and invert for consistent interpretation
hr <- summary(coxph_model)$coefficients[, "exp(coef)"]
hr <- round(1 / hr, 2)  

# Prepare for plot output: close any open graphics device
if (dev.cur() != 1) dev.off()

# Calculate median survival times per cluster
medians <- summary(kmcurve)$table[, "median"]
median_text <- paste0(names(medians), ": ", round(medians, 1), " months")

# Calculate HR and 95% CI from Cox model
cox_summary <- summary(coxph_model)
hr <- round(cox_summary$coefficients[, "exp(coef)"], 2)
hr_lower <- round(cox_summary$conf.int[, "lower .95"], 2)
hr_upper <- round(cox_summary$conf.int[, "upper .95"], 2)
hr_text <- paste0(names(hr), ": HR=", hr, " (", hr_lower, "-", hr_upper, ")")

# Save KM plot
if (dev.cur() != 1) dev.off()
tiff(filename = file.path(filepath, "OS_RFS_survival.tif"), width = 10, height = 10, units = "in", res = 300)

# Generate plot
g <- ggsurvplot(kmcurve,
                data = km_rfs_df,
                conf.int = F,
                pval = TRUE,
                risk.table = "abs_pct",  # Absolute + percentage
                risk.table.title = "Number at Risk",
                risk.table.y.text.col = TRUE,
                risk.table.y.text = FALSE,
                surv.median.line = "hv",  # Show median survival line
                legend.labs = names(risk_labels),
                legend.title = paste0("K-means Clusters (n = ", num_clusters, ")"),
                title = "Kaplan-Meier Curve: CRC Clusters",
                palette = color_palette,
                xlim = c(0, 60))


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


