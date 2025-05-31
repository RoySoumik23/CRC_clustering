# CRC_clustering

## Overview
This repository contains R scripts and supplementary data for clustering analysis in colorectal cancer (CRC) research. The primary goal is to identify molecular subtypes, mutational signatures, and prognostic patterns within CRC datasets using various clustering algorithms and survival analysis techniques.

## Features
- Identification of CRC molecular subtypes using clustering methods such as K-means and DBSCAN.
- Extraction and analysis of COSMIC mutational signatures and de novo mutational signatures.
- Survival analysis including Kaplan-Meier plots to assess prognostic significance of subtypes.
- Comprehensive visualization tools including heatmaps and cluster plots.
- Supplementary data files with clinical and genomic annotations.

## Repository Structure
- `CRC_clustering_1.R` to `CRC_clustering_5_shanky.R`: Core R scripts for clustering, signature analysis, and survival analysis.
- `ComplexHeatmap_CRC.R`: Script for generating heatmaps to visualize clustering results.
- Supplementary data files (`Supplementary_Table_01.xlsx`, `Supplementary_Table_13.xlsx`): Clinical and genomic data supporting the analysis.

## Requirements
- R version >= 4.0
- R packages: `ComplexHeatmap`, `survival`, `dbscan`, `ggplot2`, `ConsensusClusterPlus`, and other dependencies as used in the scripts.

## Usage
1. Clone the repository:
git clone [https://github.com/RoySoumik23/CRC_clustering.git]
2. Open the R scripts in your R environment.
3. Install any missing packages using `install.packages()` or Bioconductor if needed.
4. Load the supplementary data files as required by the scripts.
5. Run the scripts sequentially or individually depending on your analysis needs.

## Notes
- The scripts assume familiarity with colorectal cancer molecular subtyping frameworks such as the Consensus Molecular Subtypes (CMS).
- Input data should be prepared according to the formats expected by the scripts (e.g., mutation matrices, expression data).
- It is recommended to review and modify file paths and parameters within the scripts to suit your data.

## Citation
If you use this repository for your research, please cite the relevant publications or contact the author for citation details.

## Contact
For questions or support, please open an issue or contact the repository owner.

---

