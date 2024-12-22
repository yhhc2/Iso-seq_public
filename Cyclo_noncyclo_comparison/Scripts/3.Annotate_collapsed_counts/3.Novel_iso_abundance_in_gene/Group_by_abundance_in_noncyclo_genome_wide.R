# Hank Cheng
# yhhc@uw.edu

# Usage: 
# conda activate r_env_per_isoform
# Rscript 1.Group_by_abundance_in_noncyclo_genome_wide.R "/path/to/isoform_level.csv" "/path/to/gene_level.csv" 0.5 0.005 0.01 0.00000001


# Load necessary libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(stats) # For chisq.test
library(testthat)

# Parsing Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if (length(args) != 6) {
  stop("Usage: Rscript 1.Group_by_abundance_in_noncyclo_genome_wide.R <dt_isoform_level_file> <dt_gene_level_file> <hyp1_sig_proportion_masking_threshold> <bin_proportion> <significance_threshold> <masking_threshold>", call. = FALSE)
}

# Assign variables from arguments
dt_isoform_level <- fread(args[1])
dt_gene_level <- fread(args[2])
hyp1_sig_proportion_masking_threshold <- as.numeric(args[3]) # Use 10 if you want no masking. Use 1 if you want masking only if hyp1 is sig in ALL controls
bin_proportion <- as.numeric(args[4]) # Bin proportion for minor bin used in hypothesis 5
significance_threshold <- as.numeric(args[5])
masking_threshold <- as.numeric(args[6])


# This is used for unit testing
number_of_rows_dt_gene <- nrow(dt_gene_level)


#############################################################
# Perform bin calculations
#############################################################

print(paste("Starting bin calculations at:", Sys.time()))

#isoform_noncyclo_proportion can have a value of NaN
#This is because 0/0 = NaN. So if a gene has 0 counts,
#then the denominator is 0. Take this into account by replacing
#NA with zero.
dt_isoform_level[is.na(isoform_noncyclo_proportion), isoform_noncyclo_proportion := 0]

# isoforms that are not seen at all in the noncyclo sample are removed.
# If this is commented out then isoforms that are not seen at all in the noncyclo sample (but are seen in the cyclo sample) are considered low abundance isoforms
#dt <- dt[isoform_noncyclo_proportion != 0]

# isoforms that are not seen at all in the noncyclo sample (but are seen in the cyclo sample) are considered low abundance isoforms
# Step 2: Create a new column for binning for all genes
dt_isoform_level[, bin := ifelse(isoform_noncyclo_proportion <= bin_proportion, "Bin1_le", "Bin2_g")]


# Step 3: Aggregate the Cyclo_TPM and Noncyclo_TPM, and list Isoform_PBid for each sample and gene
result <- dt_isoform_level[, .(Isoforms = list(Isoform_PBid),
                 Total_bin_cyclo_count = sum(cyclo_count),
                 Total_bin_noncyclo_count = sum(noncyclo_count)),
             by = .(Sample, associated_gene, bin)]

rm(dt_isoform_level)

# Step 4: Pivot the data to wide format
wide_result <- dcast(result, Sample + associated_gene ~ bin, value.var = c("Total_bin_cyclo_count", "Total_bin_noncyclo_count"))

rm(result)

# Replace NA values with zero. NA values appear in the Total_bin_cyclo_count_Bin2_g and Total_bin_noncyclo_count_Bin2_g
# bins if the gene has zero counts or happens to have no isoforms that make up more than 0.01 of the gene.
# Similar thing can happen to the le bins.
wide_result[is.na(Total_bin_cyclo_count_Bin2_g), Total_bin_cyclo_count_Bin2_g := 0]
wide_result[is.na(Total_bin_noncyclo_count_Bin2_g), Total_bin_noncyclo_count_Bin2_g := 0]
wide_result[is.na(Total_bin_cyclo_count_Bin1_le), Total_bin_cyclo_count_Bin1_le := 0]
wide_result[is.na(Total_bin_noncyclo_count_Bin1_le), Total_bin_noncyclo_count_Bin1_le := 0]

print(paste("Finished bin calculations at:", Sys.time()))


#############################################################
# Perform chi-square
#############################################################

print(paste("Starting chi-square calculations at:", Sys.time()))


# Step 5: Apply the chi-square test for each Sample and gene
wide_result[, c("p_value") := {
  matrix_data <- matrix(c(Total_bin_cyclo_count_Bin1_le, Total_bin_cyclo_count_Bin2_g,
                          Total_bin_noncyclo_count_Bin1_le, Total_bin_noncyclo_count_Bin2_g),
                        nrow = 2, byrow = TRUE)
  if (any(matrix_data != 0)){
    test <- chisq.test(matrix_data)
    list(test$p.value)
  }
}, by = .(Sample, associated_gene)]
# Rename p-value to P_Value_Hyp5
setnames(wide_result, "p_value", "P_Value_Hyp5")

# Step 6: Calculate proportions for cyclo and noncyclo counts across bins for each gene
wide_result[, proportion_in_Bin1_cyclo := Total_bin_cyclo_count_Bin1_le / (Total_bin_cyclo_count_Bin2_g + Total_bin_cyclo_count_Bin1_le)]
wide_result[, proportion_in_Bin1_noncyclo := Total_bin_noncyclo_count_Bin1_le / (Total_bin_noncyclo_count_Bin2_g + Total_bin_noncyclo_count_Bin1_le)]

print(paste("Finished chi-square calculations at:", Sys.time()))

#############################################################
# Add all the columns from dt_gene_level
#############################################################

print(paste("Started add dt_gene_level columns at:", Sys.time()))


setDT(wide_result)

#names(dt_gene_level)[names(dt_gene_level) == "Isoform_PBid"] <- "associated_gene"

# Step 1: Identify the columns in dt_gene_level not in wide_result
unique_columns <- setdiff(names(dt_gene_level), names(wide_result))

# Ensure 'Sample' and 'associated_gene' are not part of the unique columns to merge
unique_columns <- setdiff(unique_columns, c("Sample", "associated_gene"))

# Step 2: Merge all unique columns from dt_gene_level to wide_result
# Specify the columns to keep from dt_gene_level (key columns + unique columns)
columns_to_merge <- c("Sample", "associated_gene", unique_columns)

# Merge operation
wide_result <- merge(wide_result, dt_gene_level[, ..columns_to_merge], by = c("Sample", "associated_gene"), all.x = TRUE)

rm(dt_gene_level)

# Create a filename using the sample name
filename <- paste0("data_combined_full_gene_with_Hyp5.csv")

# Write the data.table to a CSV file
fwrite(wide_result, file = filename)

print(paste("Finished adding dt_gene_level columns at:", Sys.time()))


#######################################################################
# Top genes based on Hyp5 p-values.
#######################################################################

print(paste("Started filtering for top genes at:", Sys.time()))


# Get a vector of unique samples
samples <- unique(wide_result$Sample)

for (sample in samples) {
  
  print(sample)
  
  # Filter current sample data
  sample_data <- wide_result[Sample == sample]
  
  # Filter control data (all samples except the current one)
  control_data <- wide_result[Sample != sample]
  
  # Find isoforms significant in the sample but not in controls. One sided significance.
  significant_isoforms <- sample_data[P_Value_Hyp5 < significance_threshold]
  
  
  # Filter isoforms that are non-significant across all controls. One sided significance.
  # na.rm because nan values are reported as minimum even though they should be considered the same as 
  # a p-value of 1.
  # Cyclo_TPM > Noncyclo_TPM.
  isoforms_significant_in_controls <- control_data[
    Cyclo_TPM > Noncyclo_TPM, 
    .(MinPValue = min(P_Value_Hyp5, na.rm = TRUE)), 
    by = .(Isoform_PBid)
  ][MinPValue < masking_threshold, .(Isoform_PBid)]
  
  
  # To remove isoforms significant in controls from significant_isoforms, perform an anti-join
  # This uses the `data.table` syntax for anti-join: `!in`
  top_isoforms <- significant_isoforms[!Isoform_PBid %in% isoforms_significant_in_controls$Isoform_PBid]

  # Using hyp1 to further remove genes using controls
  ###############
  # Check the proportion of controls where P_Values_Hyp1 for a gene are below the masking threshold and NormalizedFractionDifference > 0
  gene_significant_in_controls_hyp1 <- control_data[, .(
    ProportionSignificant = sum(P_Value_Hyp1 < masking_threshold & NormalizedFractionDifference > 0) / .N
  ), by = .(associated_gene)]

  # Filter to keep only those genes where ProportionSignificant is greater than or equal to 0.5
  significant_genes <- gene_significant_in_controls_hyp1[ProportionSignificant >= hyp1_sig_proportion_masking_threshold, associated_gene]

  # Remove these genes from wide_result_sample
  top_isoforms <- top_isoforms[!associated_gene %in% significant_genes]
  ###############
  
  
  # Sort top_isoforms by P_Value_Hyp1. Ascending order.
  top_isoforms <- top_isoforms[order(P_Value_Hyp5)]
    
  # Create a filename using the sample name
  filename <- paste0("Hyp5_Top_Isoforms_", sample, ".csv")
  
  # Write the data.table to a CSV file
  fwrite(top_isoforms, file = filename)
  
  rm(sample_data)
  rm(control_data)
  rm(significant_isoforms)
  rm(isoforms_significant_in_controls)
  rm(top_isoforms)
  
}

print(paste("Finished filtering for top genes at:", Sys.time()))

#######################################################################
# Unit tests
#######################################################################

print("Unit test results:")

# Unit tests on the final wide_result
test_that("wide_result has the correct columns", {
  expect_true(all(c("Sample", "associated_gene", "Total_bin_cyclo_count_Bin1_le", 
                    "Total_bin_noncyclo_count_Bin1_le", "Total_bin_cyclo_count_Bin2_g", 
                    "Total_bin_noncyclo_count_Bin2_g", "P_Value_Hyp5", 
                    "proportion_in_Bin1_cyclo", "proportion_in_Bin1_noncyclo") %in% names(wide_result)))
})

test_that("no NA values in key columns of wide_result", {
  key_columns <- c("Total_bin_cyclo_count_Bin1_le", "Total_bin_noncyclo_count_Bin1_le", 
                   "Total_bin_cyclo_count_Bin2_g", "Total_bin_noncyclo_count_Bin2_g")
  for(col in key_columns) {
    expect_true(all(!is.na(wide_result[[col]])))
  }
})

test_that("P_Value_Hyp5 are valid probabilities", {
  expect_true(all(is.na(wide_result$P_Value_Hyp5) | wide_result$P_Value_Hyp5 >= 0 & wide_result$P_Value_Hyp5 <= 1))
})

test_that("dt_gene_level and wide_result have the same number of rows", {
  # Assume dt_gene_level and wide_result are already defined and loaded
  
  num_rows_wide_result <- nrow(wide_result)
  
  # Test if the number of rows are the same
  expect_equal(number_of_rows_dt_gene, num_rows_wide_result, 
               info = "The number of rows in dt_gene_level should match the number of rows in wide_result.")
})