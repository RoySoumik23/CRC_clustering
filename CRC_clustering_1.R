#======================================================================
# Prognostic Genome and Transcriptome Signatures in Colorectal Cancers
#======================================================================
# PART 1: Data Preprocessing and Preparation for Clustering
#======================================================================

#======================================================================
# Number ocluster_dfs_merged
#======================================================================
num_clusters = 4
#======================================================================
# Load required libraries and set up the environment
source("https://raw.githubusercontent.com/RoySoumik23/R_scripts/refs/heads/main/libraries_and_setup.R")

# Set working directory
setwd("/home/shivalik2/Soumik/my_codes/Avik_bhaiya_project")

#======================================================================
# Load Data
#======================================================================
# Load the main dataset, skipping metadata rows
df <- read_excel("Supplementary_Table_13.xlsx", skip = 4)  # Mutation data
df_whole_data <- read_excel("Supplementary_Table_01.xlsx", skip = 2)  # Full dataset

#======================================================================
# Data Preprocessing
#======================================================================
# Rename columns for clarity
column_names <- c(
  "Sample_IDs", "Total Mutations", "Cosine Similarity", "L1 Norm", "L1 Norm %", "L2 Norm", "L2 Norm %", 
  "KL Divergence", "Correlation", "SBS96A", "SBS96B", "SBS96C", "SBS96D", "SBS96E", "SBS96F", "SBS96G", 
  "SBS96H", "SBS96I", "SBS96J", "SBS96K", "SBS96L", "SBS96M", "SBS96N", "SBS96O", "SBS96P", "SBS96Q", 
  "SBS96R", "SBS96S", "SBS96T", "SBS96U", "SBS96V", "SBS96W", "SBS96X", "SBS96Y", "SBS96Z", "SBS96AA", 
  "SBS1", "SBS2", "SBS3", "SBS5", "SBS7b", "SBS7c", "SBS8", "SBS10a", "SBS10b", "SBS13", "SBS14", "SBS15", 
  "SBS17a", "SBS17b", "SBS18", "SBS20", "SBS21", "SBS26", "SBS28", "SBS30", "SBS32", "SBS34", "SBS37", 
  "SBS42", "SBS44", "SBS54", "SBS57", "SBS60", "SBS86", "SBS88", "SBS89", "SBS90", "SBS93", "SBS94", 
  "SBS-CRC1", "SBS-CRC2", "Total Mutations", "Cosine Similarity", "L1 Norm", "L1 Norm %", "L2 Norm", "L2 Norm %", 
  "KL Divergence", "Correlation", "DBS78A", "DBS78B", "DBS78C", "DBS78D", "DBS78E", "DBS78F", "DBS78G", 
  "DBS78H", "DBS2", "DBS4", "DBS5", "DBS6", "DBS7", "DBS8", "DBS-CRC1", "DBS-CRC2", "DBS-CRC3", "DBS-CRC4", 
  "DBS-CRC5", "Total Mutations", "Cosine Similarity", "L1 Norm", "L1 Norm %", "L2 Norm", "L2 Norm %", "KL Divergence", 
  "Correlation", "ID83A", "ID83B", "ID83C", "ID83D", "ID83E", "ID83F", "ID83G", "ID83H", "ID83I", "ID83J", "ID83K", 
  "ID1", "ID2", "ID4", "ID6", "ID7", "ID8", "ID9", "ID14", "ID18", "ID-CRC1", "ID-CRC2"
)
colnames(df) <- column_names
rm(column_names)

# Convert the dataset into a standard dataframe
df <- as.data.frame(df)

# Assign first column as row names and remove it from the dataset
rownames(df) <- df[, 1]
df <- df[, -c(1,2)]  # Remove first two columns if no longer needed

# Convert all elements to numeric format
df[] <- lapply(df, as.numeric)

#======================================================================
# Generating Dataframes for Clustering Analysis
#======================================================================

# Function to create directories dynamically (optional)
create_dir <- function(folder_name) {
  if (!dir.exists(folder_name)) dir.create(folder_name)
  return(folder_name)
}


# Function to create directories dynamically
create_dir <- function(folder_name) {
  if (!dir.exists(folder_name)) dir.create(folder_name)
  return(folder_name)
}

# COSMIC Variants
df_COSMIC <- df[, c(8:68, 79:92, 106:125)]
dir_COSMIC <- create_dir("all_COSMIC_variant")

# de-novo Signatures
df_denovo <- df[, c(69:70, 93:97, 126:127)]
dir_denovo <- create_dir("all_de_novo_variants")

# ALL Signatures
df_all <- df[, c(8:70, 79:97, 106:127)]
dir_all <- create_dir("all_signatures")

# SBS_denovo Variants
df_sbs1 <- df[, 69:70]  
dir_sbs1 <- create_dir("sbs_denovo_variant")

# SBS COSMIC Signatures
df_sbs2 <- df[, 8:68]  
dir_sbs2 <- create_dir("sbs2_COSMIC_variant")

# DBS_denovo Variants
df_dbs1 <- df[, 93:97]  
dir_dbs1 <- create_dir("dbs1_denovo_variant")

# DBS COSMIC Signatures
df_dbs2 <- df[, 79:92]  
dir_dbs2 <- create_dir("dbs2_COSMIC_variant")

# ID_denovo Variants
df_id1 <- df[, 126:127]  
dir_id1 <- create_dir("id1_denovo_variant")

# ID COSMIC Signatures
df_id2 <- df[, 106:125]  
dir_id2 <- create_dir("id2_COSMIC_variant")

#======================================================================
message("End of Preprocessing")
#======================================================================