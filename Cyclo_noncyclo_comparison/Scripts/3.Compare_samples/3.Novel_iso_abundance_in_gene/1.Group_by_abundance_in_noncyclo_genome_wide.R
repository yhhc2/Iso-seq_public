# Hank Cheng
# yhhc@uw.edu

# Usage: 
# conda activate r_env_per_isoform
# Rscript Group_by_abundance_in_noncyclo_genome_wide.R > Group_by_abundance_in_noncyclo_genome_wide_output.txt 2>&1

# Load necessary libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(stats) # For chisq.test
library(testthat)

# Step 1: Read the CSV file
dt_isoform_level <- fread("/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/4.24.24_merge_aligned_bams/3.Comparison_between_samples/Isoform/data_combined_full.csv")
dt_gene_level <- fread("/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/4.24.24_merge_aligned_bams/3.Comparison_between_samples/Gene/data_combined_full.csv")

omim_file <- "/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/3.15.24_merge_aligned_bams/4.Comparison_between_samples/Combined_OMIM_Isoform_4.3.24/genemap2_7.21.23.txt"

# Hyp1 masking proportion threshold
# Use 10 if you want no masking
# Use 1 if you want masking only if hyp1 is sig in ALL controls
hyp1_sig_proportion_masking_threshold <- 0.5

# Bin proportion for minor bin used in hypothesis 5
bin_proportion <- 0.005

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

# Step 4: Pivot the data to wide format
wide_result <- dcast(result, Sample + associated_gene ~ bin, value.var = c("Total_bin_cyclo_count", "Total_bin_noncyclo_count"))


# Replace NA values with zero. NA values appear in the Total_bin_cyclo_count_Bin2_g and Total_bin_noncyclo_count_Bin2_g
# bins if the gene has zero counts or happens to have no isoforms that make up more than 0.01 of the gene.
# Similar thing can happen to the le bins.
wide_result[is.na(Total_bin_cyclo_count_Bin2_g), Total_bin_cyclo_count_Bin2_g := 0]
wide_result[is.na(Total_bin_noncyclo_count_Bin2_g), Total_bin_noncyclo_count_Bin2_g := 0]
wide_result[is.na(Total_bin_cyclo_count_Bin1_le), Total_bin_cyclo_count_Bin1_le := 0]
wide_result[is.na(Total_bin_noncyclo_count_Bin1_le), Total_bin_noncyclo_count_Bin1_le := 0]


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

# Step 6: Calculate proportions for cyclo and noncyclo counts across bins for each gene
wide_result[, proportion_in_Bin1_cyclo := Total_bin_cyclo_count_Bin1_le / (Total_bin_cyclo_count_Bin2_g + Total_bin_cyclo_count_Bin1_le)]
wide_result[, proportion_in_Bin1_noncyclo := Total_bin_noncyclo_count_Bin1_le / (Total_bin_noncyclo_count_Bin2_g + Total_bin_noncyclo_count_Bin1_le)]


#############################################################
# Add omim data
#############################################################

# some genes are not included when I read in the data for some reason. Like PDZK1
# This is fixed by including the quote="" argument. https://www.biostars.org/p/221983/
omim_data <- read.table(omim_file, header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE, quote="")

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
wide_result <- left_join(wide_result, rename(omim_data, associated_gene = "Approved Gene Symbol"), by = "associated_gene")

# Add new column that indicates if omim phenotype is associated with gene
wide_result <- wide_result %>%
  mutate(PhenotypesNotEmpty = !is.na(Phenotypes) & Phenotypes != "")

# Add the gene-level hypothesis 1 values and NormalizedFractionDifference to the datatable
setDT(wide_result)
names(dt_gene_level)[names(dt_gene_level) == "Isoform_PBid"] <- "associated_gene"
wide_result <- merge(wide_result, dt_gene_level[, .(Sample, associated_gene, P_Value_Hyp1)], by = c("Sample", "associated_gene"), all.x = TRUE)
wide_result <- merge(wide_result, dt_gene_level[, .(Sample, associated_gene, NormalizedFractionDifference)], by = c("Sample", "associated_gene"), all.x = TRUE)

# Rename p-value to P_Value_Hyp5
setnames(wide_result, "p_value", "P_Value_Hyp5")

# Pre-filter
genes_to_keep <- wide_result[P_Value_Hyp5 < 0.1 & ((proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)) > 0, unique(associated_gene)]

# Step 2: Select rows with the genes of interest and only keep specified columns
smaller_dt <- wide_result[associated_gene %in% genes_to_keep, .(
  Sample,
  associated_gene,
  P_Value_Hyp5,
  proportion_in_Bin1_cyclo,
  proportion_in_Bin1_noncyclo,
  Phenotypes,
  P_Value_Hyp1,
  NormalizedFractionDifference
)]

# Remove rows with ((proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)) <0
smaller_dt <- smaller_dt[((proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)) > 0]

#fwrite(smaller_dt, "Hyp5_gene_data_combined_filtered_0.1_BCH_UDN_4.21.24.csv")

data <- c()

data$genes_filtered_hyp5 <- smaller_dt

saveRDS(data, "hank_data_hyp5.rds")


# Do this for plotting purposes.
smaller_dt <- wide_result[associated_gene %in% genes_to_keep, .(
  Sample,
  associated_gene,
  P_Value_Hyp5,
  proportion_in_Bin1_cyclo,
  proportion_in_Bin1_noncyclo,
  Phenotypes,
  PhenotypesNotEmpty,
  P_Value_Hyp1,
  NormalizedFractionDifference
)]

#############################################################
# Plotting
#############################################################

sample_of_interest <- "UDN212054"
gene_of_interest <- "HARS1"

# Try this
p_value_thresh <- 1
masking_thresh <- 0.001

# If the gene has a significant hypothesis 5 p-value in any of the controls, remove it.
gene_significant_in_other_samples <- subset(wide_result, Sample != sample_of_interest & P_Value_Hyp5 < masking_thresh & (((proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)) > 0))
gene_significant_in_other_samples <- unique(gene_significant_in_other_samples$associated_gene)
# Only include rows from merged_data_omim_added if the Isoform column doesn't contain a value that is present in the isoforms_significant_in_samples_other_than_UDN212054 character array
wide_result_removed_common_bin_hypothesis <- subset(wide_result, !(associated_gene %in% gene_significant_in_other_samples))

# If the gene has a significant hypothesis 1 p-value in ALL of the controls, remove it.
# This would remove CD83 from UDN633333.
# Subset to get only control samples, not including the sample of interest
control_samples <- wide_result_removed_common_bin_hypothesis[Sample != sample_of_interest]
# Check if all P_Values_Hyp1 for a gene in control samples are below the masking threshold
gene_significant_in_all_controls_hyp1 <- control_samples[, .(AllSignificant = all(P_Value_Hyp1 < masking_thresh & NormalizedFractionDifference > 0)), by = .(associated_gene)]
# Filter to keep only those where AllSignificant is TRUE
significant_genes <- gene_significant_in_all_controls_hyp1[AllSignificant == TRUE, associated_gene]
# Remove these genes from wide_result_sample
wide_result_removed_common_bin_hypothesis_and_hyp1 <- wide_result_removed_common_bin_hypothesis[!associated_gene %in% significant_genes]

# Filter only for sample of interest.
wide_result_sample <- wide_result_removed_common_bin_hypothesis_and_hyp1[Sample == sample_of_interest]
wide_result_sample <- wide_result_sample[!is.na(P_Value_Hyp5)]
wide_result_sample <- wide_result_sample[P_Value_Hyp5 < p_value_thresh]

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- as.data.frame(wide_result_sample)

# Convert PhenotypesNotEmpty to a factor with specific labels
data$PhenotypesNotEmpty <- factor(data$PhenotypesNotEmpty, levels = c("FALSE", "TRUE"), labels = c("No Phenotype", "Phenotype"))

# Create the volcano plot with appropriate legend labels
volcano_plot <- ggplot(data, aes(x = (proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo), y = -log10(P_Value_Hyp5), color = PhenotypesNotEmpty)) +
  geom_point(alpha = 0.6, size = 1) +  # Color points based on PhenotypesNotEmpty
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlight gene of interest
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Significance threshold line
  labs(x = "(proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)", 
       y = "-log10(P_Value_Hyp5)", 
       title = paste(sample_of_interest, "p-value thresh=", p_value_thresh, " masking thresh=", masking_thresh)) +
  scale_color_manual(values = c("No Phenotype" = "gray", "Phenotype" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank(),  # Optionally hide the legend title
        axis.title.x = element_text(size = 7))


# Add the label for the gene of interest using geom_text()
# The hjust and vjust arguments should be negative to position the label correctly above the point
volcano_plot_with_label <- volcano_plot +
  geom_text(
    data = subset(data, associated_gene == gene_of_interest),
    aes(label = associated_gene),
    hjust = 0, vjust = -1,
    size = 2,
    color = "Green",
    fontface = "bold"
  ) #+ ylim(0, 10)
  

# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest, "mask_with_hyp1", "p-value thresh=",p_value_thresh, " masking thresh=", masking_thresh,".pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)

#############################################################
# Plotting with smaller dt
#############################################################

#smaller_dt <- fread("Hyp5_gene_data_combined_filtered_0.1_BCH_UDN_4.21.24.csv")

sample_of_interest <- "UDN633333"
gene_of_interest <- "MFN2"

# Try this
p_value_thresh <- 1
masking_thresh <- 0.00000001

# If the gene has a significant hypothesis 5 p-value in any of the controls, remove it.
gene_significant_in_other_samples <- subset(smaller_dt, Sample != sample_of_interest & P_Value_Hyp5 < masking_thresh & (((proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)) > 0))
gene_significant_in_other_samples <- unique(gene_significant_in_other_samples$associated_gene)
wide_result_removed_common_bin_hypothesis <- subset(smaller_dt, !(associated_gene %in% gene_significant_in_other_samples))

# ###########################################################################
# # If the gene has a significant hypothesis 1 p-value in ALL of the controls, remove it.
# # This would remove CD83 from UDN633333.
# # Subset to get only control samples, not including the sample of interest
# control_samples <- wide_result_removed_common_bin_hypothesis[Sample != sample_of_interest]
# # Check if all P_Values_Hyp1 for a gene in control samples are below the masking threshold
# gene_significant_in_all_controls_hyp1 <- control_samples[, .(AllSignificant = all(P_Value_Hyp1 < masking_thresh & NormalizedFractionDifference > 0)), by = .(associated_gene)]
# # Filter to keep only those where AllSignificant is TRUE
# significant_genes <- gene_significant_in_all_controls_hyp1[AllSignificant == TRUE, associated_gene]
# # Remove these genes from wide_result_sample
# wide_result_removed_common_bin_hypothesis_and_hyp1 <- wide_result_removed_common_bin_hypothesis[!associated_gene %in% significant_genes]
# ###########################################################################

###########################################################################
# If the gene has a significant hypothesis 1 p-value in HALF of the controls, remove it.
# Subset to get only control samples, not including the sample of interest
control_samples <- wide_result_removed_common_bin_hypothesis[Sample != sample_of_interest]

# Check the proportion of controls where P_Values_Hyp1 for a gene are below the masking threshold and NormalizedFractionDifference > 0
gene_significant_in_half_controls_hyp1 <- control_samples[, .(
  ProportionSignificant = sum(P_Value_Hyp1 < masking_thresh & NormalizedFractionDifference > 0) / .N
), by = .(associated_gene)]

# Filter to keep only those genes where ProportionSignificant is greater than or equal to 0.5
significant_genes <- gene_significant_in_half_controls_hyp1[ProportionSignificant >= hyp1_sig_proportion_masking_threshold, associated_gene]

# Remove these genes from wide_result_sample
wide_result_removed_common_bin_hypothesis_and_hyp1 <- wide_result_removed_common_bin_hypothesis[!associated_gene %in% significant_genes]
###########################################################################

# Filter only for sample of interest.
wide_result_sample <- wide_result_removed_common_bin_hypothesis_and_hyp1[Sample == sample_of_interest]
wide_result_sample <- wide_result_sample[!is.na(P_Value_Hyp5)]
wide_result_sample <- wide_result_sample[P_Value_Hyp5 < p_value_thresh]

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- as.data.frame(wide_result_sample)

# Convert PhenotypesNotEmpty to a factor with specific labels
data$PhenotypesNotEmpty <- factor(data$PhenotypesNotEmpty, levels = c("FALSE", "TRUE"), labels = c("No Phenotype", "Phenotype"))

# Create the volcano plot with appropriate legend labels
volcano_plot <- ggplot(data, aes(x = (proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo), y = -log10(P_Value_Hyp5), color = PhenotypesNotEmpty)) +
  geom_point(alpha = 0.6, size = 1) +  # Color points based on PhenotypesNotEmpty
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlight gene of interest
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Significance threshold line
  labs(x = "(proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)", 
       y = "-log10(P_Value_Hyp5)", 
       title = paste(sample_of_interest, "p-value thresh=", p_value_thresh, " masking thresh=", masking_thresh)) +
  scale_color_manual(values = c("No Phenotype" = "gray", "Phenotype" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank(),  # Optionally hide the legend title
        axis.title.x = element_text(size = 7))


# Add the label for the gene of interest using geom_text()
# The hjust and vjust arguments should be negative to position the label correctly above the point
volcano_plot_with_label <- volcano_plot +
  geom_text(
    data = subset(data, associated_gene == gene_of_interest),
    aes(label = associated_gene),
    hjust = 0, vjust = -1,
    size = 2,
    color = "Green",
    fontface = "bold"
  ) #+ ylim(0, 10)


# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest, "mask_with_hyp1", "p-value thresh=",p_value_thresh, " masking thresh=", masking_thresh,".pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)