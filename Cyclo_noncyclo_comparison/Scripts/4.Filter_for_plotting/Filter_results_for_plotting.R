# Hank Cheng

# Load necessary libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(stats) # For chisq.test
library(testthat)

#Usage: Rscript Filter_results_for_plotting.R <isoform_file> <gene_file>


############################################################
# Define inputs
############################################################

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  print("Usage: Rscript Filter_results_for_plotting.R <isoform_file> <gene_file>")
  stop("Error: Incorrect number of arguments provided.")
}

# Assign arguments to variables
isoform_file <- args[1]
gene_file <- args[2]

# Read the CSV files into datatables using the input file paths
data_combined_full_isoform <- fread(isoform_file)
data_combined_full_gene <- fread(gene_file)


############################################################
# Isoform level
############################################################


# Assuming data_combined_full_isoform is already a data.table
# Step 1: Filter isoforms based on the P_Value_Hyp1 condition
# effectively identifies and keeps the IDs of isoforms for which any sample has a 
# P_Value_Hyp1 less than 0.01. This line filters data_combined_full_isoform to only those 
# rows where P_Value_Hyp1 is less than 0.01 and then applies unique(Isoform_PBid) 
# to get a list of unique isoform IDs that meet this condition.
isoforms_to_keep_hyp1 <- data_combined_full_isoform[P_Value_Hyp1 < 0.1 & NormalizedFractionDifference > 0, unique(Isoform_PBid)]
isoforms_to_keep_hyp2_below_median <- data_combined_full_isoform[Max_P_Value_Hyp2_below_median < 0.1, unique(Isoform_PBid)]
isoforms_to_keep_hyp2_above_median <- data_combined_full_isoform[Max_P_Value_Hyp2_above_median < 0.1, unique(Isoform_PBid)]
isoforms_to_keep_hyp3_below_median <- data_combined_full_isoform[Max_P_Value_Hyp3_below_median < 0.1, unique(Isoform_PBid)]
isoforms_to_keep_hyp3_above_median <- data_combined_full_isoform[Max_P_Value_Hyp3_above_median < 0.1, unique(Isoform_PBid)]

# Step 2: Select rows with the isoforms of interest and only keep specified columns
isoforms_filtered_hyp1 <- data_combined_full_isoform[Isoform_PBid %in% isoforms_to_keep_hyp1, .(
  Isoform_PBid,
  Sample,
  associated_gene,
  P_Value_Hyp1,
  Cyclo_TPM,
  Noncyclo_TPM,
  NormalizedFractionDifference,
  Phenotypes
)]

# Hypothesis 2 Below Median
isoforms_filtered_hyp2_below_median <- data_combined_full_isoform[Isoform_PBid %in% isoforms_to_keep_hyp2_below_median, .(
  Isoform_PBid,
  Sample,
  associated_gene,
  Max_P_Value_Hyp2_below_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Noncyclo_Z_Score, #Take this out and replace with appropriate x-a
  Phenotypes
)]

# Hypothesis 2 Above Median
isoforms_filtered_hyp2_above_median <- data_combined_full_isoform[Isoform_PBid %in% isoforms_to_keep_hyp2_above_median, .(
  Isoform_PBid,
  Sample,
  associated_gene,
  Max_P_Value_Hyp2_above_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Noncyclo_Z_Score,
  Phenotypes
)]

# Hypothesis 3 Below Median
isoforms_filtered_hyp3_below_median <- data_combined_full_isoform[Isoform_PBid %in% isoforms_to_keep_hyp3_below_median, .(
  Isoform_PBid,
  Sample,
  associated_gene,
  Max_P_Value_Hyp3_below_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Cyclo_Z_Score,
  Phenotypes
)]

# Hypothesis 3 Above Median
isoforms_filtered_hyp3_above_median <- data_combined_full_isoform[Isoform_PBid %in% isoforms_to_keep_hyp3_above_median, .(
  Isoform_PBid,
  Sample,
  associated_gene,
  Max_P_Value_Hyp3_above_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Cyclo_Z_Score,
  Phenotypes
)]

# Remove rows with NormalizedFractionDifference <0
isoforms_filtered_hyp1 <- isoforms_filtered_hyp1[NormalizedFractionDifference > 0]

############################################################
# Gene level
############################################################

# Assuming data_combined_full_gene is already a data.table
# Step 1: Filter isoforms based on the P_Value_Hyp1 condition
# effectively identifies and keeps the IDs of isoforms for which any sample has a 
# P_Value_Hyp1 less than 0.01. This line filters data_combined_full_gene to only those 
# rows where P_Value_Hyp1 is less than 0.01 and then applies unique(Isoform_PBid) 
# to get a list of unique isoform IDs that meet this condition.
genes_to_keep_hyp1 <- data_combined_full_gene[P_Value_Hyp1 < 0.1 & NormalizedFractionDifference > 0, unique(Isoform_PBid)]
genes_to_keep_hyp2_below_median <- data_combined_full_gene[Max_P_Value_Hyp2_below_median < 0.1, unique(Isoform_PBid)]
genes_to_keep_hyp2_above_median <- data_combined_full_gene[Max_P_Value_Hyp2_above_median < 0.1, unique(Isoform_PBid)]
genes_to_keep_hyp3_below_median <- data_combined_full_gene[Max_P_Value_Hyp3_below_median < 0.1, unique(Isoform_PBid)]
genes_to_keep_hyp3_above_median <- data_combined_full_gene[Max_P_Value_Hyp3_above_median < 0.1, unique(Isoform_PBid)]
#genes_to_keep_hyp5 <- data_combined_full_gene[Max_P_Value_Hyp3_above_median < 0.1, unique(Isoform_PBid)]

# Step 2: Select rows with the isoforms of interest and only keep specified columns
genes_filtered_hyp1 <- data_combined_full_gene[Isoform_PBid %in% genes_to_keep_hyp1, .(
  Sample,
  associated_gene,
  P_Value_Hyp1,
  Cyclo_TPM,
  Noncyclo_TPM,
  NormalizedFractionDifference,
  Phenotypes
)]

# Hypothesis 2 Below Median
genes_filtered_hyp2_below_median <- data_combined_full_gene[Isoform_PBid %in% genes_to_keep_hyp2_below_median, .(
  Sample,
  associated_gene,
  Max_P_Value_Hyp2_below_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Noncyclo_Z_Score,
  Phenotypes
)]

# Hypothesis 2 Above Median
genes_filtered_hyp2_above_median <- data_combined_full_gene[Isoform_PBid %in% genes_to_keep_hyp2_above_median, .(
  Sample,
  associated_gene,
  Max_P_Value_Hyp2_above_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Noncyclo_Z_Score,
  Phenotypes
)]

# Hypothesis 3 Below Median
genes_filtered_hyp3_below_median <- data_combined_full_gene[Isoform_PBid %in% genes_to_keep_hyp3_below_median, .(
  Sample,
  associated_gene,
  Max_P_Value_Hyp3_below_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Cyclo_Z_Score,
  Phenotypes
)]

# Hypothesis 3 Above Median
genes_filtered_hyp3_above_median <- data_combined_full_gene[Isoform_PBid %in% genes_to_keep_hyp3_above_median, .(
  Sample,
  associated_gene,
  Max_P_Value_Hyp3_above_median,
  Cyclo_TPM,
  Noncyclo_TPM,
  Cyclo_Z_Score,
  Phenotypes
)]

# Remove rows with NormalizedFractionDifference <0
genes_filtered_hyp1 <- genes_filtered_hyp1[NormalizedFractionDifference > 0]


############################################################
# Gene level hypothesis 5
############################################################

# Step 1: Identify genes to keep based on P_Value_Hyp5 threshold
genes_to_keep_hyp5 <- data_combined_full_gene[P_Value_Hyp5 < 0.1, unique(Isoform_PBid)]

# Step 2: Filter data for selected isoforms and retain specified columns
genes_filtered_hyp5 <- data_combined_full_gene[Isoform_PBid %in% genes_to_keep_hyp5, .(
  Sample,
  associated_gene,
  P_Value_Hyp5,
  proportion_in_Bin1_cyclo,
  proportion_in_Bin1_noncyclo,
  Cyclo_TPM,
  Noncyclo_TPM,
  NormalizedFractionDifference,
  Phenotypes,
  P_Value_Hyp1
)]


############################################################
# Make into R object for easy loading
############################################################

data <- c()

data$isoforms_filtered_hyp1 <- isoforms_filtered_hyp1
data$isoforms_filtered_hyp2_below_median <- isoforms_filtered_hyp2_below_median
data$isoforms_filtered_hyp2_above_median <- isoforms_filtered_hyp2_above_median
data$isoforms_filtered_hyp3_below_median <- isoforms_filtered_hyp3_below_median
data$isoforms_filtered_hyp3_above_median <- isoforms_filtered_hyp3_above_median

data$genes_filtered_hyp1 <- genes_filtered_hyp1
data$genes_filtered_hyp2_below_median <- genes_filtered_hyp2_below_median
data$genes_filtered_hyp2_above_median <- genes_filtered_hyp2_above_median
data$genes_filtered_hyp3_below_median <- genes_filtered_hyp3_below_median
data$genes_filtered_hyp3_above_median <- genes_filtered_hyp3_above_median
data$genes_filtered_hyp5 <- genes_filtered_hyp5


saveRDS(data, "data_for_plotting.rds")