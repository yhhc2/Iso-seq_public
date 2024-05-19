
# Load necessary libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(stats) # For chisq.test
library(testthat)

# Usage: 
# conda activate r_env_per_isoform
# Rscript Filter_results_for_shiny.R > Filter_results_for_shiny_output.txt 2>&1

############################################################
# Isoform level
############################################################

# Read the CSV file into a datatable
data_combined_full_isoform <- fread("/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/4.24.24_merge_aligned_bams/3.Comparison_between_samples/Isoform/data_combined_full.csv")

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

# Read the CSV file into a datatable
data_combined_full_gene <- fread("/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/4.24.24_merge_aligned_bams/3.Comparison_between_samples/Gene/data_combined_full.csv")

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

saveRDS(data, "hank_data.rds")