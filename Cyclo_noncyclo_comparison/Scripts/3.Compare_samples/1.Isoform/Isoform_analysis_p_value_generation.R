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
significance_threshold <- as.numeric(args[4])
masking_threshold <- as.numeric(args[5])
count_threshold <- as.numeric(args[6])


# Load necessary libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(stats) # For chisq.test
library(testthat)

sample_info <- read.delim(sample_info_path, header = TRUE, sep = "\t")
sample_info_unique <- distinct(sample_info, collapsed_by_isoform_file_cyclo_noncyclo_counts_classified, .keep_all = TRUE)


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

# Calculate total counts
total_counts_by_sample <- data_combined_full %>%
  group_by(Sample) %>%
  summarise(TotalCyclo = sum(cyclo_count), TotalNoncyclo = sum(noncyclo_count))

# Join `total_counts_by_sample` with `data_combined_full` to have total counts available for each row
data_with_totals <- left_join(data_combined_full, total_counts_by_sample, by = "Sample")

# Convert to data.table
setDT(data_with_totals)

# Add TPM calculations to the data.table
data_with_totals[, `:=` (
  Cyclo_TPM = (cyclo_count / TotalCyclo) * 1000000,
  Noncyclo_TPM = (noncyclo_count / TotalNoncyclo) * 1000000
)]

# For each isoform, add the median TPM. 
# I have a datatable with columns Sample, Isoform_PBid, Cyclo_TPM, and Noncyclo_TPM.
# For each unique value of Isoform_PBid, I want to calculate the median Cyclo_TPM and
# the median Noncyclo_TPM. I want to add these median values for the Isoforms as two
# new columns to the datatable.
# Calculate the median TPMs for each Isoform_PBid
median_tpm <- data_with_totals[, .(
  MedianCyclo_TPM = median(Cyclo_TPM, na.rm = TRUE),
  MedianNoncyclo_TPM = median(Noncyclo_TPM, na.rm = TRUE)
), by = Isoform_PBid]

# Merge the median values back into the original datatable
data_with_totals <- merge(data_with_totals, median_tpm, by = "Isoform_PBid", all.x = TRUE)

rm(sample_data_list)
rm(all_combinations)

#######################################################################
# Hypothesis 1
#######################################################################

print(paste("Starting hypothesis 1 p-value calculations at:", Sys.time()))

# Print helpful statements
print("Number of unique isoforms:")
print(length(unique(data_combined$Isoform_PBid)))

# Setup parallel backend
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Calculate chunk size and split data into chunks
chunk_size <- ceiling(nrow(data_with_totals) / no_cores)
chunks <- split(data_with_totals, (seq_len(nrow(data_with_totals)) - 1) %/% chunk_size + 1)


# Parallel computation of p-values for each chunk
p_values_list <- foreach(chunk = chunks, .packages = c("data.table", "stats")) %dopar% {
  # First p-value calculation using cyclo_count
  chunk[, P_Value_Hyp1 := mapply(function(cyclo_count, TotalCyclo, noncyclo_count, TotalNoncyclo) {
    chisq.test(matrix(c(cyclo_count, TotalCyclo, noncyclo_count, TotalNoncyclo), nrow = 2))$p.value
  }, cyclo_count, TotalCyclo, noncyclo_count, TotalNoncyclo, SIMPLIFY = TRUE)]
  
  return(chunk)
}

# Combine chunks back into the full dataset
data_with_totals_Hyp1 <- rbindlist(p_values_list)

# Stop the cluster
stopCluster(cl)

# Select the desired columns. In data.table, you can directly assign the result to a new object
results_hypothesis_1 <- data_with_totals_Hyp1[, .(Isoform_PBid, Sample, P_Value_Hyp1)]

# Write to CSV
write.csv(results_hypothesis_1, "Results_hypothesis_1.csv", row.names = FALSE)

# Remove to save space in environment
rm(p_values_list)
rm(data_with_totals_Hyp1)
rm(chunk_size)
rm(chunks)

print(paste("Completed hypothesis 1 p-value calculations at:", Sys.time()))


#######################################################################
# Hypotheses 2 and 3
#######################################################################

print(paste("Starting hypothesis 2 and 3 p-value calculations at:", Sys.time()))

# To save memory, remove any isoforms that have <count_threshold for all samples.
isoforms_to_keep <- data_with_totals[, .(
  keep = any(cyclo_count >= count_threshold | noncyclo_count >= count_threshold)
), by = Isoform_PBid][keep == TRUE, Isoform_PBid]

# Print helpful statements
print("Number of isoforms to keep for hyp2/3 after filtering based on counts:")
print(length(isoforms_to_keep))

# Filter out the data_with_totals to keep only the desired Isoform_PBid
data_with_totals <- data_with_totals[Isoform_PBid %in% isoforms_to_keep]

# Cross join data_with_totals by itself to get all combinations
data_with_totals_all_combos <- data_with_totals[data_with_totals, on = .(Isoform_PBid = Isoform_PBid), allow.cartesian = TRUE]

# Rename columns to reflect the combination nature
setnames(data_with_totals_all_combos, old = c("Sample", "Isoform_PBid", "cyclo_count", "noncyclo_count", "TotalCyclo", "TotalNoncyclo", "Cyclo_TPM", "Noncyclo_TPM", "MedianCyclo_TPM", "MedianNoncyclo_TPM", 
                                              "i.Sample", "i.cyclo_count", "i.noncyclo_count", "i.TotalCyclo", "i.TotalNoncyclo", "i.Cyclo_TPM", "i.Noncyclo_TPM", "i.MedianCyclo_TPM", "i.MedianNoncyclo_TPM"),
         new = c("Sample1", "Isoform_PBid", "cyclo_count1", "noncyclo_count1", "TotalCyclo1", "TotalNoncyclo1", "Cyclo_TPM1", "Noncyclo_TPM1", "MedianCyclo_TPM1", "MedianNoncyclo_TPM1", 
                 "Sample2", "cyclo_count2", "noncyclo_count2", "TotalCyclo2", "TotalNoncyclo2", "Cyclo_TPM2", "Noncyclo_TPM2", "MedianCyclo_TPM2", "MedianNoncyclo_TPM2"))

# Convert Sample1 and Sample2 columns to character type
set(data_with_totals_all_combos, j = "Sample1", value = as.character(data_with_totals_all_combos$Sample1))
set(data_with_totals_all_combos, j = "Sample2", value = as.character(data_with_totals_all_combos$Sample2))

# Filter out the same sample combinations to ensure Sample1 != Sample2
data_with_totals_all_combos <- data_with_totals_all_combos[Sample1 != Sample2]

# Remove duplicate combinations where (Sample1, Sample2) is considered the same as (Sample2, Sample1)
# This step assumes that the order of samples in a pair does not matter for your analysis
data_with_totals_all_combos <- data_with_totals_all_combos[Sample1 < Sample2]

print("Number of comparisons:")
print(nrow(data_with_totals_all_combos))

original_nrow <- nrow(data_with_totals_all_combos)

# Setup parallel backend
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Calculate chunk size and split data into chunks
chunk_size <- ceiling(nrow(data_with_totals_all_combos) / no_cores)
chunks <- split(data_with_totals_all_combos, (seq_len(nrow(data_with_totals_all_combos)) - 1) %/% chunk_size + 1)

# Parallel computation of p-values for each chunk
p_values_list <- foreach(chunk = chunks, .packages = c("data.table", "stats")) %dopar% {
  
  # Below median
  
  # Add a new column for P_Value_Hyp3, initialized to NULL or NA
  chunk[, P_Value_Hyp3_below_median := NA_character_]
  # Apply the chi-squared test only if CycloTPM1 and CycloTPM2 are below their respective MedianCycloTPM
  chunk[, P_Value_Hyp3_below_median := ifelse(
    Cyclo_TPM1 < MedianCyclo_TPM1 & Cyclo_TPM2 < MedianCyclo_TPM1,
    mapply(function(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2) {
      test_matrix <- matrix(c(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2), nrow = 2)
      chisq.test(test_matrix)$p.value
    }, cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2, SIMPLIFY = TRUE),
    NA_real_ #"Both_TPM_are_not_below_median"
  )]
  
  # Adding initialization of P_Value_Hyp2 to a default message or NA
  chunk[, P_Value_Hyp2_below_median := NA_character_]
  # Apply the chi-squared test conditionally based on NonCycloTPM values
  chunk[, P_Value_Hyp2_below_median := ifelse(
    Noncyclo_TPM1 < MedianNoncyclo_TPM1 & Noncyclo_TPM2 < MedianNoncyclo_TPM1,
    mapply(function(noncyclo_count1, TotalNoncyclo1, noncyclo_count2, TotalNoncyclo2) {
      test_matrix <- matrix(c(noncyclo_count1, TotalNoncyclo1, noncyclo_count2, TotalNoncyclo2), nrow = 2)
      chisq.test(test_matrix)$p.value
    }, noncyclo_count1, TotalNoncyclo1, noncyclo_count2, TotalNoncyclo2, SIMPLIFY = TRUE),
    NA_real_ #"Both_TPM_are_not_below_median"
  )]
  
  
  # Above median
  
  # Add a new column for P_Value_Hyp3, initialized to NULL or NA
  chunk[, P_Value_Hyp3_above_median := NA_character_]
  # Apply the chi-squared test only if CycloTPM1 and CycloTPM2 are below their respective MedianCycloTPM
  chunk[, P_Value_Hyp3_above_median := ifelse(
    Cyclo_TPM1 > MedianCyclo_TPM1 & Cyclo_TPM2 > MedianCyclo_TPM1,
    mapply(function(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2) {
      test_matrix <- matrix(c(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2), nrow = 2)
      chisq.test(test_matrix)$p.value
    }, cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2, SIMPLIFY = TRUE),
    NA_real_ #"Both_TPM_are_not_above_median"
  )]
  
  # Adding initialization of P_Value_Hyp2 to a default message or NA
  chunk[, P_Value_Hyp2_above_median := NA_character_]
  # Apply the chi-squared test conditionally based on NonCycloTPM values
  chunk[, P_Value_Hyp2_above_median := ifelse(
    Noncyclo_TPM1 > MedianNoncyclo_TPM1 & Noncyclo_TPM2 > MedianNoncyclo_TPM1,
    mapply(function(noncyclo_count1, TotalNoncyclo1, noncyclo_count2, TotalNoncyclo2) {
      test_matrix <- matrix(c(noncyclo_count1, TotalNoncyclo1, noncyclo_count2, TotalNoncyclo2), nrow = 2)
      chisq.test(test_matrix)$p.value
    }, noncyclo_count1, TotalNoncyclo1, noncyclo_count2, TotalNoncyclo2, SIMPLIFY = TRUE),
    NA_real_ #"Both_TPM_are_not_above_median"
  )]
  
  return(chunk)
}

# Combine chunks back into the full dataset
data_with_totals_all_combos <- rbindlist(p_values_list)
rm(p_values_list)

# Stop the cluster
stopCluster(cl)

# Select the desired columns. In data.table, you can directly assign the result to a new object
results_hypothesis_2 <- data_with_totals_all_combos[, .(Isoform_PBid, Sample1, Sample2, P_Value_Hyp2_below_median, P_Value_Hyp2_above_median)]
# Write to CSV
write.csv(results_hypothesis_2, "Results_hypothesis_2.csv", row.names = FALSE)

# Select the desired columns. In data.table, you can directly assign the result to a new object
results_hypothesis_3 <- data_with_totals_all_combos[, .(Isoform_PBid, Sample1, Sample2, P_Value_Hyp3_below_median, P_Value_Hyp3_above_median)]
# Write to CSV
write.csv(results_hypothesis_3, "Results_hypothesis_3.csv", row.names = FALSE)

end_time <- Sys.time()
print(paste("End time:", end_time))

#######################################################################
# Unit tests for p-value generation
#######################################################################

print("Unit test results for p-value generation:")

# Test: Data transformations are correct
test_that("Data transformations are correct", {
  # Test median calculations
  expect_true("MedianCyclo_TPM" %in% names(data_with_totals))
  expect_true("MedianNoncyclo_TPM" %in% names(data_with_totals))
  
  # Ensure that after joining, no key columns are filled with NA
  expect_equal(sum(is.na(data_with_totals$Sample)), 0)
  expect_equal(sum(is.na(data_with_totals$Isoform_PBid)), 0)
})

# Data joins do not create duplicate rows
test_that("Data joins do not create duplicate rows", {
  original_n <- nrow(data_with_totals)
  data_with_totals <- merge(data_with_totals, median_tpm, by = "Isoform_PBid", all.x = TRUE)
  expect_equal(nrow(data_with_totals), original_n)
})




# Test: Sample information is loaded correctly
test_that("Sample information is loaded correctly", {
  expect_true(nrow(sample_info) > 0)
  expect_true("id" %in% colnames(sample_info))
  expect_true("collapsed_by_isoform_file_cyclo_noncyclo_counts_classified" %in% colnames(sample_info))
})

# Test: Data combined correctly
test_that("Data combined correctly", {
  expect_true(nrow(data_combined) > 0)
  expect_true("Sample" %in% colnames(data_combined))
})

# Test: Verify the number of rows in data_with_totals_all_combos is as expected
test_that("Number of rows in data_with_totals_all_combos is correct", {
  num_samples <- length(unique(data_with_totals$Sample))
  num_isoforms <- length(unique(data_with_totals$Isoform_PBid))
  expected_num_combinations_per_isoform <- choose(num_samples, 2)
  expected_total_rows <- num_isoforms * expected_num_combinations_per_isoform
  
  actual_total_rows <- nrow(data_with_totals_all_combos)
  
  expect_equal(actual_total_rows, expected_total_rows)
})

# Test: Total counts are calculated correctly
test_that("Total counts are calculated correctly", {
  total_counts_check <- data_combined_full %>%
    group_by(Sample) %>%
    summarise(check_TotalCyclo = sum(cyclo_count), check_TotalNoncyclo = sum(noncyclo_count))
  
  expect_equal(total_counts_by_sample$TotalCyclo, total_counts_check$check_TotalCyclo)
  expect_equal(total_counts_by_sample$TotalNoncyclo, total_counts_check$check_TotalNoncyclo)
})


# Test: P-value calculation is correct for each sample combination
test_that("P-value calculation is correct for each sample combination", {
  suppressWarnings({
    # Create a small mock dataset to mimic the expected structure
    mock_data <- data.table(
      Sample1 = c("Sample1", "Sample2"),
      Isoform_PBid = c("Iso1", "Iso1"),
      cyclo_count1 = c(10, 0),
      TotalCyclo1 = c(100, 50),
      Sample2 = c("Sample2", "Sample3"),
      cyclo_count2 = c(20, 15),
      TotalCyclo2 = c(200, 150)
    )
    # Manually compute expected p-values
    expected_p_values <- c(
      chisq.test(matrix(c(10, 100, 20, 200), nrow = 2))$p.value,
      chisq.test(matrix(c(0, 50, 15, 150), nrow = 2))$p.value
    )
    # Calculate p-values for mock data
    mock_data[, P_Value := mapply(function(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2) {
      chisq.test(matrix(c(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2), nrow = 2))$p.value
    }, cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2, SIMPLIFY = TRUE), by = .(Isoform_PBid, Sample1, Sample2)]
    
    # Test
    expect_equal(mock_data$P_Value, expected_p_values)
    
    # Expecting the first p-value to be non-significant.
    expect_gt(mock_data$P_Value[1], 0.1)
    # Expecting the last p-value to be significant.
    expect_lt(mock_data$P_Value[2], 0.1)
    
  })
})


# Test the parallel processing. Do chunks add up to the original data
# Unit Test: Verify that recombined chunks match the original data table size and content
test_that("Recombined chunks match original data table", {
  # Check if the number of rows matches
  expect_equal(original_nrow, nrow(data_with_totals_all_combos))
  
})


# Unit test: Additional tests to ensure that hypothesis corrections are done correctly."
test_that("Additional tests to ensure that hypothesis corrections are done correctly.", {
  suppressWarnings({
    # Create sample data
    data_with_totals_sample_hyp1_3 <- data.frame(
      Isoform_PBid = rep("PB.29164.115", 9),
      Sample = c("BCH_1877-01", "BCH_2275-01", "BCH_2410-01", "BCH_2443-01", "UDN052264", 
                 "UDN212054", "UDN215640", "UDN633333", "UDN687128"),
      cyclo_count = c(10, 0, 0, 0, 0, 74, 0, 0, 0),
      noncyclo_count = c(10, 0, 0, 0, 0, 0, 0, 0, 0),
      TotalCyclo = c(8514101, 5298875, 5115261, 7983565, 4388293, 
                     5582703, 5138817, 4965770, 4951492),
      TotalNoncyclo = c(8390416, 6226664, 5604044, 8045115, 5339570, 
                        5443395, 5563024, 4950758, 4938020)
    )
    
    data_with_totals_sample_hyp1_2 <- data.frame(
      Isoform_PBid = rep("PB.29164.115", 9),
      Sample = c("BCH_1877-01", "BCH_2275-01", "BCH_2410-01", "BCH_2443-01", "UDN052264", 
                 "UDN212054", "UDN215640", "UDN633333", "UDN687128"),
      cyclo_count = c(10, 0, 0, 0, 0, 0, 0, 0, 0),
      noncyclo_count = c(10, 0, 0, 0, 0, 74, 0, 0, 0),
      cyclo_total_sum = c(8514101, 5298875, 5115261, 7983565, 4388293, 
                          5582703, 5138817, 4965770, 4951492),
      noncyclo_total_sum = c(8390416, 6226664, 5604044, 8045115, 5339570, 
                             5443395, 5563024, 4950758, 4938020)
    )
    
    setDT(data_with_totals_sample_hyp1_3)
    setDT(data_with_totals_sample_hyp1_2)
    
    # Calculate hypothesis 1 p-values.
    data_with_totals_sample_hyp1_3[, P_Value_Hyp1 := mapply(function(cyclo_count, TotalCyclo, noncyclo_count, TotalNoncyclo) {
      chisq.test(matrix(c(cyclo_count, TotalCyclo, noncyclo_count, TotalNoncyclo), nrow = 2))$p.value
    }, cyclo_count, TotalCyclo, noncyclo_count, TotalNoncyclo, SIMPLIFY = TRUE)]
    
    expect_lt(data_with_totals_sample_hyp1_3[6,P_Value_Hyp1], 0.01)
    expect_gt(data_with_totals_sample_hyp1_3[1,P_Value_Hyp1], 0.01)
    expect_true(is.nan(data_with_totals_sample_hyp1_3[2,P_Value_Hyp1]))
    
    # Cross join data_with_totals by itself to get all combinations
    data_with_totals_sample_hyp1_3_all_combos <- data_with_totals_sample_hyp1_3[data_with_totals_sample_hyp1_3, on = .(Isoform_PBid = Isoform_PBid), allow.cartesian = TRUE]
    # Rename columns to reflect the combination nature
    setnames(data_with_totals_sample_hyp1_3_all_combos, old = c("Sample", "Isoform_PBid", "cyclo_count", "noncyclo_count", "TotalCyclo", "TotalNoncyclo", "i.Sample", "i.cyclo_count", "i.noncyclo_count", "i.TotalCyclo", "i.TotalNoncyclo"),
             new = c("Sample1", "Isoform_PBid", "cyclo_count1", "noncyclo_count1", "TotalCyclo1", "TotalNoncyclo1", "Sample2", "cyclo_count2", "noncyclo_count2", "TotalCyclo2", "TotalNoncyclo2"))
    # Convert Sample1 and Sample2 columns to character type
    set(data_with_totals_sample_hyp1_3_all_combos, j = "Sample1", value = as.character(data_with_totals_sample_hyp1_3_all_combos$Sample1))
    set(data_with_totals_sample_hyp1_3_all_combos, j = "Sample2", value = as.character(data_with_totals_sample_hyp1_3_all_combos$Sample2))
    # Filter out the same sample combinations to ensure Sample1 != Sample2
    data_with_totals_sample_hyp1_3_all_combos <- data_with_totals_sample_hyp1_3_all_combos[Sample1 != Sample2]
    # Remove duplicate combinations where (Sample1, Sample2) is considered the same as (Sample2, Sample1)
    # This step assumes that the order of samples in a pair does not matter for your analysis
    data_with_totals_sample_hyp1_3_all_combos <- data_with_totals_sample_hyp1_3_all_combos[Sample1 < Sample2]
    
    # Calculate hypothesis 3 p-values
    data_with_totals_sample_hyp1_3_all_combos[, P_Value_Hyp3 := mapply(function(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2) {
      chisq.test(matrix(c(cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2), nrow = 2))$p.value
    }, cyclo_count1, TotalCyclo1, cyclo_count2, TotalCyclo2, SIMPLIFY = TRUE)]
    
    
    # Filter rows where Sample1 or Sample2 is UDN212054
    relevant_rows <- data_with_totals_sample_hyp1_3_all_combos[Sample1 == "UDN212054" | Sample2 == "UDN212054", ]
    # Test if all P_Value_Hyp3 values in these rows are < 0.00001
    expect_true(all(relevant_rows$P_Value_Hyp3 < 0.00001))
    
    relevant_rows <- data_with_totals_sample_hyp1_3_all_combos[(Sample1 == "BCH_1877-01" | Sample2 == "BCH_1877-01") & Sample1 != "UDN212054" & Sample2 != "UDN212054", ]
    # Test if all P_Value_Hyp3 values in these rows are > 0.00001
    expect_true(all(relevant_rows$P_Value_Hyp3 > 0.00001))
    
    relevant_rows <- data_with_totals_sample_hyp1_3_all_combos[Sample1 != "BCH_1877-01" & Sample2 != "BCH_1877-01" & Sample1 != "UDN212054" & Sample2 != "UDN212054", ]
    # Test if all P_Value_Hyp3 values in these rows are > 0.00001
    expect_true(all(is.nan(relevant_rows$P_Value_Hyp3)))
    
  })
})

# Assuming data_with_totals_all_combos is your data frame or data table
test_that("MedianCycloTPM1 equals MedianCycloTPM2 in all rows", {
  # Calculate the number of rows where MedianCycloTPM1 does not equal MedianCycloTPM2
  non_equal_count <- sum(data_with_totals_all_combos$MedianCycloTPM1 != data_with_totals_all_combos$MedianCycloTPM2)
  
  # Expect that there are no rows with differing MedianCycloTPM1 and MedianCycloTPM2
  expect_equal(non_equal_count, 0)
})


rm(data_with_totals)
rm(data_with_totals_all_combos)
rm(chunks)

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
# Add hypothesis 1 p-values
#######################################################################

print(paste("Started adding hypothesis 1 p-values at:", Sys.time()))

# Merge df1 and df2
data_combined_full <- merge(data_combined_full, as.data.frame(results_hypothesis_1), by = c("Isoform_PBid", "Sample"), all.x = TRUE)

print(paste("Finished adding hypothesis 1 p-values at:", Sys.time()))

rm(results_hypothesis_1)

#######################################################################
# Add hypotheis 2 and 3 worst p-value
#######################################################################

print(paste("Started adding hypothesis 2 and 3 worst p-values at:", Sys.time()))

# Hypothesis 2 below median
# First, ensure Sample1 and Sample2 are in a long format for easier grouping
results_long_hyp2_below_median <- melt(results_hypothesis_2, id.vars = c("Isoform_PBid", "P_Value_Hyp2_below_median"), measure.vars = c("Sample1", "Sample2"), value.name = "Sample")
# Now, calculate max P-value for each Isoform_PBid and Sample. If the sample has an NaN value for a comparison, then
# it means that it has (0,0) compared to (0,0) and we want to output NaN as the p-value. 
results_max_hyp2_below_median <- results_long_hyp2_below_median[, .(Max_P_Value_Hyp2_below_median = if (any(is.nan(P_Value_Hyp2_below_median))) NaN else max(P_Value_Hyp2_below_median, na.rm = TRUE)), by = .(Isoform_PBid, Sample)]
data_combined_full <- merge(data_combined_full, results_max_hyp2_below_median, by = c("Isoform_PBid", "Sample"), all.x = TRUE)

# Hypothesis 2 above median
# First, ensure Sample1 and Sample2 are in a long format for easier grouping
results_long_hyp2_above_median <- melt(results_hypothesis_2, id.vars = c("Isoform_PBid", "P_Value_Hyp2_above_median"), measure.vars = c("Sample1", "Sample2"), value.name = "Sample")
# Now, calculate max P-value for each Isoform_PBid and Sample. If the sample has an NaN value for a comparison, then
# it means that it has (0,0) compared to (0,0) and we want to output NaN as the p-value. 
results_max_hyp2_above_median <- results_long_hyp2_above_median[, .(Max_P_Value_Hyp2_above_median = if (any(is.nan(P_Value_Hyp2_above_median))) NaN else max(P_Value_Hyp2_above_median, na.rm = TRUE)), by = .(Isoform_PBid, Sample)]
data_combined_full <- merge(data_combined_full, results_max_hyp2_above_median, by = c("Isoform_PBid", "Sample"), all.x = TRUE)


# Hypothesis 3 below median
# First, ensure Sample1 and Sample2 are in a long format for easier grouping
results_long_hyp3_below_median <- melt(results_hypothesis_3, id.vars = c("Isoform_PBid", "P_Value_Hyp3_below_median"), measure.vars = c("Sample1", "Sample2"), value.name = "Sample")
# Now, calculate max P-value for each Isoform_PBid and Sample. If the sample has an NaN value for a comparison, then
# it means that it has (0,0) compared to (0,0) and we want to output NaN as the p-value. 
results_max_hyp3_below_median <- results_long_hyp3_below_median[, .(Max_P_Value_Hyp3_below_median = if (any(is.nan(P_Value_Hyp3_below_median))) NaN else max(P_Value_Hyp3_below_median, na.rm = TRUE)), by = .(Isoform_PBid, Sample)]
data_combined_full <- merge(data_combined_full, results_max_hyp3_below_median, by = c("Isoform_PBid", "Sample"), all.x = TRUE)

# Hypothesis 3 above median
# First, ensure Sample1 and Sample2 are in a long format for easier grouping
results_long_hyp3_above_median <- melt(results_hypothesis_3, id.vars = c("Isoform_PBid", "P_Value_Hyp3_above_median"), measure.vars = c("Sample1", "Sample2"), value.name = "Sample")
# Now, calculate max P-value for each Isoform_PBid and Sample. If the sample has an NaN value for a comparison, then
# it means that it has (0,0) compared to (0,0) and we want to output NaN as the p-value. 
results_max_hyp3_above_median <- results_long_hyp3_above_median[, .(Max_P_Value_Hyp3_above_median = if (any(is.nan(P_Value_Hyp3_above_median))) NaN else max(P_Value_Hyp3_above_median, na.rm = TRUE)), by = .(Isoform_PBid, Sample)]
data_combined_full <- merge(data_combined_full, results_max_hyp3_above_median, by = c("Isoform_PBid", "Sample"), all.x = TRUE)


setDT(data_combined_full)


# List of columns to replace -Inf with NA
# -Inf occurs when na.rm removes all values.
columns_to_replace <- c("Max_P_Value_Hyp2_below_median", 
                        "Max_P_Value_Hyp2_above_median", 
                        "Max_P_Value_Hyp3_below_median", 
                        "Max_P_Value_Hyp3_above_median")
# Replace -Inf with NA in each specified column
for (col in columns_to_replace) {
  # Check if the column exists to avoid errors
  if (col %in% names(data_combined_full)) {
    # Replace -Inf values with NA
    data_combined_full[is.infinite(get(col)) & get(col) == -Inf, (col) := NA]
  }
}

# For these 4 columns for Hyp2 and Hyp3.
# NaN means there was a (0,0) comparison to (0,0)
# NA means all the values are NA. Ex. The sample has a TPM that is above the median, but the column is the below median column.

print(paste("Finished adding hypothesis 2 and 3 max fisher p-values at:", Sys.time()))

rm(results_hypothesis_2)
rm(results_hypothesis_3)

rm(results_max_hyp2_above_median)
rm(results_max_hyp2_below_median)
rm(results_max_hyp3_above_median)
rm(results_max_hyp3_below_median)

rm(results_long_hyp2_above_median)
rm(results_long_hyp2_below_median)
rm(results_long_hyp3_above_median)
rm(results_long_hyp3_below_median)


#######################################################################
# Add total sums and TPM
#######################################################################

print(paste("Started adding total sums, TPM, etc. at:", Sys.time()))

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

rm(isoform_averages)
rm(median_tpm)



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

# Handle pseudoautosomal genes. These genes have two rows in the omim data. For each gene, I want to collapse them into a single row.
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
# For each sample, select the top isoforms based on p-values.
#######################################################################

#######################################################################
# Top isoforms based on Hyp1 p-values.
#######################################################################

print(paste("Started finding top isoforms based on Hyp1 at:", Sys.time()))

# Hypothesis 1. Using controls. Mask.
# If hypothesis 1 is not significant in controls AND is significant in the sample AND cycloTPM > noncycloTPM. 
# Initialize an empty list to store results for each sample
top_isoforms_list <- list()

# Get a vector of unique samples
samples <- unique(data_combined_full$Sample)

for (sample in samples) {
  
  print(sample)
  
  # Filter current sample data
  sample_data <- data_combined_full[Sample == sample]
  
  # Filter control data (all samples except the current one)
  control_data <- data_combined_full[Sample != sample]
  
  # Find isoforms significant in the sample but not in controls. One sided significance.
  significant_isoforms <- sample_data[P_Value_Hyp1 < significance_threshold & Cyclo_TPM > Noncyclo_TPM]
  
  
  # Filter isoforms that are non-significant across all controls
  # isoforms_non_significant_in_controls <- control_data[, .(MinPValue = min(P_Value_Hyp1, na.rm = TRUE)), by = .(Isoform_PBid)
  #                                                      ][MinPValue > masking_threshold, .(Isoform_PBid)]
  
  # Filter isoforms that are non-significant across all controls. One sided significance.
  # na.rm because nan values are reported as minimum even though they should be considered the same as 
  # a p-value of 1.
  # Cyclo_TPM > Noncyclo_TPM.
  isoforms_significant_in_controls <- control_data[
    Cyclo_TPM > Noncyclo_TPM, 
    .(MinPValue = min(P_Value_Hyp1, na.rm = TRUE)), 
    by = .(Isoform_PBid)
  ][MinPValue < masking_threshold, .(Isoform_PBid)]
  
  
  # To remove isoforms significant in controls from significant_isoforms, perform an anti-join
  # This uses the `data.table` syntax for anti-join: `!in`
  top_isoforms <- significant_isoforms[!Isoform_PBid %in% isoforms_significant_in_controls$Isoform_PBid]
  
  
  # Sort top_isoforms by P_Value_Hyp1. Ascending order.
  top_isoforms <- top_isoforms[order(P_Value_Hyp1)]
  
  # Store in the list with the sample name as the list name
  # top_isoforms_list[[sample]] <- top_isoforms
  
  # Create a filename using the sample name
  filename <- paste0("Hyp1_Top_Isoforms_", sample, ".csv")
  
  # Write the data.table to a CSV file
  fwrite(top_isoforms, file = filename)

  rm(sample_data)
  rm(control_data)
  rm(significant_isoforms)
  rm(isoforms_significant_in_controls)
  rm(top_isoforms)

}


# Combine all results into one data.table
# top_isoforms_combined <- rbindlist(top_isoforms_list, idcol = "Sample")

# Create a filename using the sample name
# filename <- paste0("Hyp1_Top_Isoforms_Combined.csv")
# Write the data.table to a CSV file
# fwrite(top_isoforms_combined, file = filename)

print(paste("Finished finding top isoforms based on Hyp1 at:", Sys.time()))

#######################################################################
# Top isoforms based on Hyp2 and Hyp3 p-values.
#######################################################################

print(paste("Started finding top isoforms based on Hyp2 and Hyp3 at:", Sys.time()))

# Hypothesis 2 and 3.
# If hypothesis 2 worst p-value is significant. Use Avg_TPM to see the direction of significance. 

for (sample in samples) {
  
  print(sample)
  
  # Filter current sample data
  sample_data <- data_combined_full[Sample == sample]
  
  # Find isoforms significant in the sample
  significant_isoforms_Hyp2_below_median <- sample_data[Max_P_Value_Hyp2_below_median < significance_threshold]
  # Ascending order.
  significant_isoforms_Hyp2_below_median <- significant_isoforms_Hyp2_below_median[order(Max_P_Value_Hyp2_below_median)]
  # Create a filename using the sample name
  filename <- paste0("Hyp2_below_median_Top_Isoforms_", sample, ".csv")
  # Write the data.table to a CSV file
  fwrite(significant_isoforms_Hyp2_below_median, file = filename)
  rm(significant_isoforms_Hyp2_below_median)
  
  # Find isoforms significant in the sample
  significant_isoforms_Hyp2_above_median <- sample_data[Max_P_Value_Hyp2_above_median < significance_threshold]
  # Ascending order.
  significant_isoforms_Hyp2_above_median <- significant_isoforms_Hyp2_above_median[order(Max_P_Value_Hyp2_above_median)]
  # Create a filename using the sample name
  filename <- paste0("Hyp2_above_median_Top_Isoforms_", sample, ".csv")
  # Write the data.table to a CSV file
  fwrite(significant_isoforms_Hyp2_above_median, file = filename)
  rm(significant_isoforms_Hyp2_above_median)
  
  # Find isoforms significant in the sample
  significant_isoforms_Hyp3_below_median <- sample_data[Max_P_Value_Hyp3_below_median < significance_threshold]
  # Ascending order.
  significant_isoforms_Hyp3_below_median <- significant_isoforms_Hyp3_below_median[order(Max_P_Value_Hyp3_below_median)]
  # Create a filename using the sample name
  filename <- paste0("Hyp3_below_median_Top_Isoforms_", sample, ".csv")
  # Write the data.table to a CSV file
  fwrite(significant_isoforms_Hyp3_below_median, file = filename)
  rm(significant_isoforms_Hyp3_below_median)
  
  # Find isoforms significant in the sample
  significant_isoforms_Hyp3_above_median <- sample_data[Max_P_Value_Hyp3_above_median < significance_threshold]
  # Ascending order.
  significant_isoforms_Hyp3_above_median <- significant_isoforms_Hyp3_above_median[order(Max_P_Value_Hyp3_above_median)]
  # Create a filename using the sample name
  filename <- paste0("Hyp3_above_median_Top_Isoforms_", sample, ".csv")
  # Write the data.table to a CSV file
  fwrite(significant_isoforms_Hyp3_above_median, file = filename)
  rm(significant_isoforms_Hyp3_above_median)

  rm(sample_data)
  
}

# The p-value from the chi-square won't differentiate the direction of significance. 
# sample cyclo count = 100
# control1 cyclo count = 0
# control2 cyclo count = 200
# sample vs control1 = significant 
# sample vs control2 = significant

# If hypothesis 3 worst p-value is significant. Use Avg_TPM column to see the direction of significance. 

print(paste("Finished finding top isoforms based on Hyp2 and Hyp3 at:", Sys.time()))


#######################################################################
# Unit tests.
#######################################################################

print("Unit test results for final table construction and filtering:")

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

test_that("hypothesis 1 p-values are correctly added", {
  expect_true("P_Value_Hyp1" %in% names(data_combined_full))
  # Optionally, check for correct merging behavior (no unexpected NAs, etc.)
})

#test_that("Fisher combined p-values are correctly calculated and added", {
#  expect_true(all(c("Fisher_P_Value_Hyp2", "Fisher_P_Value_Hyp3") %in% names(data_combined_full)))
#  # You can add more tests to check the calculation accuracy by comparing with expected values
#})



# Simulate data_combined_full with sample data
data_combined_full_sample <- data.table(
  Isoform_PBid = c("Isoform1", "Isoform2", "Isoform3", "Isoform1", "Isoform2", "Isoform3"),
  Sample = c("Sample1", "Sample1", "Sample1", "Sample2", "Sample2", "Sample2"),
  P_Value_Hyp1 = c(0.005, 0.02, 0.03, 0.001, 0.15, 0.05), # p-values for hypothesis 1
  Cyclo_TPM = c(100, 150, 10, 200, 50, 20), # Cyclo counts as TPM
  Noncyclo_TPM = c(80, 140, 15, 180, 60, 25) # Noncyclo counts as TPM
)

# Simulate control data (mock)
control_data <- data.table(
  Isoform_PBid = c("Isoform1", "Isoform2", "Isoform3"),
  Sample = "Control",
  P_Value_Hyp1 = c(0.06, 0.04, 0.06),
  Cyclo_TPM = c(90, 130, 12),
  Noncyclo_TPM = c(95, 135, 17)
)


# Filtering logic
sample_data <- data_combined_full_sample[Sample == "Sample1"]
isoforms_non_significant_in_controls <- control_data[P_Value_Hyp1 > masking_threshold, unique(Isoform_PBid)]
significant_isoforms <- sample_data[Isoform_PBid %in% isoforms_non_significant_in_controls & P_Value_Hyp1 < significance_threshold & Cyclo_TPM > Noncyclo_TPM]

# Assuming the logic puts the results in a list as per your original script
top_isoforms_list <- list("Sample1" = significant_isoforms)

test_that("top isoforms are correctly selected for each sample", {
  known_top_isoforms <- c("Isoform1") # Adjusted based on our mock dataset and filtering logic
  computed_top_isoforms <- unique(top_isoforms_list[["Sample1"]]$Isoform_PBid)
  expect_equal(sort(computed_top_isoforms), sort(known_top_isoforms))
})


# Unit test to check for unique Isoform_PBid and Sample combinations
test_that("Isoform_PBid and Sample combinations are unique", {
  # Create a subset with just the key columns
  key_subset <- unique(data_combined_full[, .(Isoform_PBid, Sample)])
  
  # Check if the number of rows in the unique subset matches the original
  expect_equal(nrow(key_subset), nrow(data_combined_full))
})


# Test to check if there are any unannotated isoforms. 
# Unannotated isoforms. Missing annotation in the pigeon file. 
if(!gene_level){
  data_combined_fill_unannotated <- data_combined_full[associated_gene == "" | is.na(associated_gene)]
  
  fwrite(data_combined_fill_unannotated, file = "unannotated_isoforms.csv")
  
  if(nrow(data_combined_fill_unannotated) > 0){
    print("Warning: There are unannotated isoforms. See the unannotated_isoforms.csv file")
  }
}
