# Hank Cheng
# yhhc@uw.edu

# Start timer
start_time <- Sys.time()
print(paste("Start time:", start_time))

# Usage: 
# conda activate r_env_per_isoform
# Rscript 1.Isoform_analysis_p_value_generation.R "/path/to/sample_info.tsv" "/path/to/omim_file.txt" FALSE 0.01 0.01 6
#Note: You might need an ultramem compute node. salloc -A stergachislab --mem=500G --time="1-0" -p compute-ultramem -c 18

# Parsing Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if (length(args) != 6) {
  stop("Usage: Rscript 1.Isoform_analysis_p_value_generation.R <sample_info_path> <omim_file_path> <gene_level> <significance_threshold> <masking_threshold> <count_threshold>", call. = FALSE)
}

# Assign variables from arguments
sample_info_path <- args[1]
omim_file_path <- args[2]
gene_level <- as.logical(args[3])  # Expecting 'TRUE' or 'FALSE'

# Load necessary libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(stats) # For chisq.test
library(testthat)

sample_info <- read.delim(sample_info_path, header = TRUE, sep = "\t")

# Function to check if all values in a row are NA or empty strings
is_all_empty <- function(row) {
  all(is.na(row) | row == "")
}
# Remove rows where all columns are NA or empty strings. Sometimes empty rows are added to the end by accident.
# This will remove these empty rows.
sample_info <- sample_info[!apply(sample_info, 1, is_all_empty), ]

sample_info_unique <- distinct(sample_info, collapsed_by_isoform_file_cyclo_noncyclo_counts_classified, .keep_all = TRUE)

print(sample_info_unique$collapsed_by_isoform_file_cyclo_noncyclo_counts_classified)

# Load data for each sample
sample_data_list <- lapply(1:nrow(sample_info_unique), function(i) {
  df <- read.table(sample_info_unique$collapsed_by_isoform_file_cyclo_noncyclo_counts_classified[i], 
                   fill = TRUE, header = TRUE, sep = " ", check.names = FALSE)
  df$Sample <- sample_info_unique$id[i]
  return(df)
})
data_combined <- bind_rows(sample_data_list)

# Perform the Gene-Level aggregation if analysis is performed per-gene instead
# of per-isoform.
if(gene_level){
  data_combined <- data_combined %>%
    group_by(Sample, associated_gene) %>%
    summarize(
      cyclo_count = sum(cyclo_count, na.rm = TRUE),
      noncyclo_count = sum(noncyclo_count, na.rm = TRUE),
      .groups = 'drop'  # This option removes the grouping structure after summarizing
    )
  data_combined <- data_combined %>%
    rename(Isoform_PBid = associated_gene)
}

# Handle missing isoforms
all_combinations <- expand.grid(Sample = unique(data_combined$Sample), 
                                Isoform_PBid = unique(data_combined$Isoform_PBid))
data_combined_full <- merge(all_combinations, data_combined, all = TRUE)
data_combined_full[is.na(data_combined_full)] <- 0
data_combined_full <- data_combined_full[, 1:4] 

# Test: Data combined correctly
test_that("Data combined correctly", {
  expect_true(nrow(data_combined) > 0)
  expect_true("Sample" %in% colnames(data_combined))
})


#######################################################################
# Populating datatable
#######################################################################

#######################################################################
# Add pigeon annotations
#######################################################################

# Only perform if analysis is at isoform level
if(!gene_level){
  
  print(paste("Starting adding pigeon annotations at:", Sys.time()))
  
  # Handle missing isoforms
  all_combinations <- expand.grid(Sample = unique(data_combined$Sample), 
                                  Isoform_PBid = unique(data_combined$Isoform_PBid))
  data_combined_full <- merge(all_combinations, data_combined, all = TRUE)
  data_combined_full[is.na(data_combined_full)] <- 0
  
  rm(data_combined)
  
  # Step 1: Identify non-zero rows for each Isoform
  non_zero_rows <- data_combined_full %>%
    filter(structural_category != 0 | associated_gene != 0 | associated_transcript != 0 | subcategory != 0) %>%
    group_by(Isoform_PBid) %>%
    slice(1) %>% # Assuming the first non-zero row per group is representative
    select(Isoform_PBid, structural_category, associated_gene, associated_transcript, subcategory)
  
  # Step 2: Fill in zeros in the main dataframe
  data_combined_full <- data_combined_full %>%
    left_join(non_zero_rows, by = "Isoform_PBid", suffix = c("", ".nonzero")) %>%
    mutate(structural_category = if_else(structural_category == 0, structural_category.nonzero, structural_category),
           associated_gene = if_else(associated_gene == 0, associated_gene.nonzero, associated_gene),
           associated_transcript = if_else(associated_transcript == 0, associated_transcript.nonzero, associated_transcript),
           subcategory = if_else(subcategory == 0, subcategory.nonzero, subcategory)) %>%
    select(-ends_with(".nonzero")) # Remove the temporary non-zero columns
  
  print(paste("Finished adding pigeon annotations at:", Sys.time()))

  rm(non_zero_rows)
  
}

#######################################################################
# Add total sums and TPM
#######################################################################

print(paste("Started adding total sums, TPM, etc. at:", Sys.time()))

# Convert data_combined_full to a data.table
setDT(data_combined_full)

# Gene-level change.If doing analysis gene-level instead of isoform-level
if(gene_level){
  data_combined_full[, associated_gene := Isoform_PBid]
}

# Group by Sample and calculate sums
data_combined_full[, `:=` (cyclo_total_sum = sum(cyclo_count),
                           noncyclo_total_sum = sum(noncyclo_count)), by = .(Sample)]

# Step 1: Add fractions and TPM calculations
data_combined_full[, `:=` (
  CycloFraction = cyclo_count / cyclo_total_sum,
  NoncycloFraction = noncyclo_count / noncyclo_total_sum,
  Cyclo_TPM = (cyclo_count / cyclo_total_sum) * 1000000,
  Noncyclo_TPM = (noncyclo_count / noncyclo_total_sum) * 1000000
)]

# Step 2: Calculate normalized fractions and their difference
data_combined_full[, `:=` (
  NormalizedCycloFraction = CycloFraction / (CycloFraction + NoncycloFraction),
  NormalizedNoncycloFraction = NoncycloFraction / (CycloFraction + NoncycloFraction)
)]

data_combined_full[, `:=` (
  NormalizedFractionDifference = NormalizedCycloFraction - NormalizedNoncycloFraction
)]

# Calculate the average Cyclo_TPM and Noncyclo_TPM for each isoform
isoform_averages <- data_combined_full[, .(Avg_Cyclo_TPM = mean(Cyclo_TPM, na.rm = TRUE),
                                           Avg_Noncyclo_TPM = mean(Noncyclo_TPM, na.rm = TRUE)), by = .(Isoform_PBid)]

# Join this back to the original data.table by creating a key column for a faster merge
data_combined_full[isoform_averages, `:=` (Avg_Cyclo_TPM = i.Avg_Cyclo_TPM, 
                                           Avg_Noncyclo_TPM = i.Avg_Noncyclo_TPM), on = .(Isoform_PBid)]



# Calculate gene-level total counts for cyclo and noncyclo treatments
data_combined_full[, `:=` (gene_cyclo_total = sum(cyclo_count),
                           gene_noncyclo_total = sum(noncyclo_count)), by = .(Sample, associated_gene)]

# Calculate the proportion of reads the isoform has for the gene from the cyclo and noncyclo treatments
data_combined_full[, `:=` (isoform_cyclo_proportion = cyclo_count / gene_cyclo_total,
                           isoform_noncyclo_proportion = noncyclo_count / gene_noncyclo_total)]


# Add Z scores
data_combined_full[, `:=` (
  Avg_Cyclo_TPM = mean(Cyclo_TPM, na.rm = TRUE),
  Avg_Noncyclo_TPM = mean(Noncyclo_TPM, na.rm = TRUE),
  SD_Cyclo_TPM = sd(Cyclo_TPM, na.rm = TRUE),
  SD_Noncyclo_TPM = sd(Noncyclo_TPM, na.rm = TRUE)
), by = Isoform_PBid][
  , `:=` (
    Cyclo_Z_Score = (Cyclo_TPM - Avg_Cyclo_TPM) / SD_Cyclo_TPM,
    Noncyclo_Z_Score = (Noncyclo_TPM - Avg_Noncyclo_TPM) / SD_Noncyclo_TPM
  )
]

# Calculate the median TPMs for each Isoform_PBid
median_tpm <- data_combined_full[, .(
  MedianCyclo_TPM = median(Cyclo_TPM, na.rm = TRUE),
  MedianNoncyclo_TPM = median(Noncyclo_TPM, na.rm = TRUE)
), by = Isoform_PBid]
# Merge the median values back into the original datatable
data_combined_full <- merge(data_combined_full, median_tpm, by = "Isoform_PBid", all.x = TRUE)

# Add the rank of the TPMs
# I have a datatable with these columns Isoform_PBid, Sample, Cyclo_TPM, Noncyclo_TPM. I want to add
# Two more columns: Cyclo_TPM_Rank and Noncyclo_TPM_Rank. Cyclo_TPM_Rank and Noncyclo_TPM_Rank are calculated
# for each unique Isoform_PBid.
# Add ranking columns based on Cyclo_TPM and Noncyclo_TPM within each Isoform_PBid
data_combined_full[, `:=` (
  Cyclo_TPM_Rank = frank(Cyclo_TPM, ties.method = "average"),
  Noncyclo_TPM_Rank = frank(Noncyclo_TPM, ties.method = "average")
), by = Isoform_PBid]

print(paste("Finished adding total sums, TPM, etc. at:", Sys.time()))



#######################################################################
# Add OMIM info
#######################################################################

print(paste("Started adding omim info at:", Sys.time()))


# some genes are not included when I read in the data for some reason. Like PDZK1
# This is fixed by including the quote="" argument. https://www.biostars.org/p/221983/
omim_data <- read.table(omim_file_path, header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE, quote="")

# Remove any rows where the approved gene symbol is empty. 
omim_data <- omim_data %>%
  filter(`Approved Gene Symbol` != "")

# Handle pseudoautosomal genes (or any gene with multiple rows in omim file). These genes have multiple rows in the omim data. For each gene, I want to collapse them into a single row.
omim_data <- omim_data %>%
  group_by(`Approved Gene Symbol`) %>%  # Ensure the column name is correct
  summarize(
    across(
      .cols = everything(),  # Apply the function to all character columns except "Approved Gene Symbol"
      .fns = ~paste(na.omit(.x), collapse = "; "),  # Concatenate, omitting NA values
      .names = "{.col}"  # Keeps the original column names
    ),
    .groups = 'drop'  # Drop grouping structure after summarization
  )

# Handle the Phenotypes column of pseudoautosomal genes. Don't want ";" to be mistaken as
# a valid phenotype.
# Assuming omim_data is your dataframe
# Replace ';' with '' in the Phenotypes column
omim_data$Phenotypes[omim_data$Phenotypes == "; "] <- ""

# Merge data frames based on the gene name column
data_combined_full <- left_join(data_combined_full, rename(omim_data, associated_gene = "Approved Gene Symbol"), by = "associated_gene")

# Add new column that indicates if omim phenotype is associated with gene
data_combined_full <- data_combined_full %>%
  mutate(PhenotypesNotEmpty = !is.na(Phenotypes) & Phenotypes != "")

setDT(data_combined_full)

# Create a filename
filename <- paste0("data_combined_full.csv")
# Write the data.table to a CSV file
fwrite(data_combined_full, file = filename)

print(paste("Finished adding omim info at:", Sys.time()))


#######################################################################
# Unit tests.
#######################################################################

print("Unit test results for final table construction and filtering:")

# Test: Sample information is loaded correctly
test_that("Sample information is loaded correctly", {
  expect_true(nrow(sample_info) > 0)
  expect_true("id" %in% colnames(sample_info))
  expect_true("collapsed_by_isoform_file_cyclo_noncyclo_counts_classified" %in% colnames(sample_info))
})


# Test: Files are read correctly
test_that("Files are read correctly", {
  expect_true(nrow(sample_info) > 0)
  expect_true(nrow(omim_data) > 0)
})

test_that("data_combined_full is correctly populated", {
  expect_true(nrow(data_combined_full) > 0)
  expect_true(all(c("Sample", "Isoform_PBid") %in% names(data_combined_full)))
  # Check for no NAs in key columns
  expect_equal(sum(is.na(data_combined_full$Sample)), 0)
  expect_equal(sum(is.na(data_combined_full$Isoform_PBid)), 0)
})

# Only perform if doing isoform-level analysis
if(!gene_level){
  test_that("pigeon annotations are correctly added", {
    # Assuming there's an original state of data_combined_full for comparison
    expect_true(all(c("structural_category", "associated_gene", "associated_transcript", "subcategory") %in% names(data_combined_full)))
    # More specific tests can be added depending on known expected outcomes
  })
}


# Test to check if there are any unannotated isoforms. 
# Unannotated isoforms. Missing annotation in the pigeon file. 
if(!gene_level){
  data_combined_fill_unannotated <- data_combined_full[associated_gene == "" | is.na(associated_gene)]
  
  fwrite(data_combined_fill_unannotated, file = "unannotated_isoforms.csv")
  
  if(nrow(data_combined_fill_unannotated) > 0){
    print("Warning: There are unannotated isoforms. See the unannotated_isoforms.csv file")
  }
}
