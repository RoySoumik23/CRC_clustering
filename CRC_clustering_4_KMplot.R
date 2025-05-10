#======================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
# PART 4: Kaplan-Meier Survival Analysis (Dynamic Cluster-Safe Version)
#======================================================================

# Prepare the dataframe for survival analysis
km_df <- kmean_clusters_df

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
tiff(filename = file.path(filepath, "survival.tif"), width = 8, height = 10, units = "in", res = 300)

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
