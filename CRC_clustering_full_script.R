#======================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
#======================================================================
# MAIN ANALYSIS SCRIPT: CRC Clustering Pipeline
#======================================================================
# Author: SOUMIK ROY (BT23MTECH11010 @CGnT Lab IITHyderabad)
# Date: [2025-02-13]
# Description: 
#   This script serves as the main entry point for clustering analyses 
#   of colorectal cancer (CRC) genomic signatures. It sequentially 
#   executes sub-scripts that process SBS (single base substitutions), 
#   DBS (double base substitutions), and ID (indels) mutation signatures.
#
#   The analysis involves:
#   1. Data pre-processing and proportion calculations
#   2. Normalization
#   3. Non-negative Matrix Factorization (NMF) for cluster estimation
#   4. K-Means clustering (n = 6) of samples
#   5. Kaplan-Meier (KM) survival analysis
#   6. Visualization (Sankey diagrams, heatmaps, etc.)
#
#======================================================================

#======================================================================
# Load Required Scripts for CRC Mutation Signature Clustering
#======================================================================

#----------------------------------------------------------------------
# Step 1: Process Single Base Substitution (SBS) Signatures
#----------------------------------------------------------------------

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# SBS_denovo Signature Analysis
message("Running SBS_denovo signature clustering...")
source("CRC_clustering_2_SBS_DENOVO.R")

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# SBS COSMIC Signature Analysis
message("Running SBS COSMIC signature clustering...")
source("CRC_Clustering_2_SBSCOSMIC.R")

#----------------------------------------------------------------------
# Step 2: Process Double Base Substitution (DBS) Signatures
#----------------------------------------------------------------------

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# DBD_denovo Signature Analysis
message("Running DBD_denovo signature clustering...")
source("CRC_clustering_2_DBS_DENOVO.R")

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# DBS COSMIC Signature Analysis
message("Running DBS COSMIC signature clustering...")
source("CRC_clustering_2_DBSCOSMIC.R")

#----------------------------------------------------------------------
# Step 3: Process Indel (ID) Signatures
#----------------------------------------------------------------------

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# ID_DENOVO Signature Analysis
message("Running ID DENOVO signature clustering...")
source("CRC_clustering_2_ID_DENOVO.R")

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# ID COSMIC Signature Analysis
message("Running ID COSMIC signature clustering...")
source("CRC_clustering_2_IDCOSMIC.R")

#----------------------------------------------------------------------
# Step 4: Process clubbed Signatures
#----------------------------------------------------------------------

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# All COSMIC Signatures Analysis
message("Running all COSMIC signature clustering...")
source("CRC_clustering_2_ALL_COSMIC.R")

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# All denovo Signatures Analysis
message("Running all denovo signature clustering...")
source("CRC_clustering_2_ALL_DENOVO.R")

# Set working directory
setwd("/home/shivalik/Downloads/R_codes/")
# All Signatures Analysis
message("Running all signature clustering...")
source("CRC_clustering_2_ALL_SIGS.R")
#======================================================================
# End of CRC Clustering Analysis
#======================================================================

message("CRC clustering analyses completed successfully!")
