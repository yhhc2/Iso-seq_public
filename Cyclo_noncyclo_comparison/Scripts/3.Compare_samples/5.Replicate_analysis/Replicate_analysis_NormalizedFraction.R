# Hank Cheng
# yhhc@uw.edu

# Usage: 
# conda activate r_env_per_isoform
# Rscript Replicate_analysis_NormalizedFraction <file_path> <count_threshold> <plot_title> <samples> <highlight_isoforms>

#setwd("~/2022-2023/Research/Computational/Isoform counting/Per_isoform_analysis/4.28.24_Comparison/Isoseq_comparison_cyclo_vs_noncyclo_4.28.24/PCA")

# Load the data.table package
library(data.table)
library(ggplot2)

# Parsing Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if (length(args) != 5) {
  stop("Usage: Rscript Replicate_analysis_NormalizedFraction <file_path> <count_threshold> <plot_title> <samples> <highlight_isoforms>", call. = FALSE)
}

# Assign variables from arguments
file_path <- args[1]
count_threshold <- as.numeric(args[2])
plot_title <- args[3]
samples <- strsplit(args[4], ",")[[1]]  # Parse samples argument
highlight_isoforms <- strsplit(args[5], ",")[[1]]  # Parse isoforms to highlight argument

# Read the CSV file into a data table
dt <- fread(file_path)

# Print the first few rows of the data table
print(head(dt))

# Ensure that the relevant columns exist
if (!("Sample" %in% colnames(dt)) | !("noncyclo_count" %in% colnames(dt)) | !("cyclo_count" %in% colnames(dt))) {
  stop("The required 'Sample', 'noncyclo_count', or 'cyclo_count' columns are not found in the data.")
}

# Filter isoforms to keep only those that meet the threshold for either noncyclo or cyclo counts in both samples
isoforms_to_keep <- dt[Sample %in% samples & (noncyclo_count >= count_threshold | cyclo_count >= count_threshold), unique(Isoform_PBid)]

# Ensure that the filtered isoforms are present in both samples for either noncyclo or cyclo counts
for (sample in samples) {
  isoforms_to_keep <- intersect(isoforms_to_keep, dt[Sample == sample & (noncyclo_count >= count_threshold | cyclo_count >= count_threshold), unique(Isoform_PBid)])
}

# Print helpful statements
print("Number of isoforms to keep after filtering based on counts:")
print(length(isoforms_to_keep))

# Filter out the data table to keep only the desired Isoform_PBid and the selected samples
dt <- dt[Isoform_PBid %in% isoforms_to_keep & Sample %in% samples]

################################################
# Generate expression matrix
################################################

# Reshape the data table to create the expression matrix
dt_long <- melt(dt, id.vars = c("Isoform_PBid", "Sample"), measure.vars = c("NormalizedFractionDifference"))

print(head(dt_long))

# Cast the long format back to a wide format where each gene has its own row and each sample its own column
expression_matrix <- dcast(dt_long, Isoform_PBid ~ Sample, value.var = "value")

# Replace all NA values with 0
expression_matrix[is.na(expression_matrix)] <- 0

# Print the first few rows of the expression matrix to check
print(head(expression_matrix))

################################################
# Generate scatter plots for the selected samples
################################################

sample1 <- samples[1]
sample2 <- samples[2]

x_col <- grep(paste0("^", sample1, "$"), colnames(expression_matrix), value = TRUE)
y_col <- grep(paste0("^", sample2, "$"), colnames(expression_matrix), value = TRUE)

# Check if both columns are found
if (length(x_col) == 0 | length(y_col) == 0) {
  stop("Required columns for scatter plot not found in the expression matrix.")
}

# Ensure the columns are not empty and handle missing values
x_values <- expression_matrix[[x_col]]
y_values <- expression_matrix[[y_col]]

if (all(is.na(x_values)) | all(is.na(y_values))) {
  stop("One or both columns are entirely NA.")
}

# Add a column to identify highlighted isoforms
expression_matrix$highlight <- ifelse(expression_matrix$Isoform_PBid %in% highlight_isoforms, "highlight", "normal")

# Calculate Spearman and Pearson correlation coefficients
spearman_corr <- cor(x_values, y_values, use = "complete.obs", method = "spearman")
pearson_corr <- cor(x_values, y_values, use = "complete.obs", method = "pearson")

# Create scatter plot
plot <- ggplot(expression_matrix, aes_string(x = x_col, y = y_col, color = "highlight")) +
  geom_point() +
  geom_text(aes(label = ifelse(Isoform_PBid %in% highlight_isoforms, Isoform_PBid, "")), vjust = -1, hjust = 1) +
  scale_color_manual(values = c("normal" = "black", "highlight" = "red")) +
  labs(
    title = paste0(plot_title, " (Spearman: ", round(spearman_corr, 2), ", Pearson: ", round(pearson_corr, 2), ")"),
    x = paste0("NormalizedFractionDifference of ", x_col),
    y = paste0("NormalizedFractionDifference of ", y_col)
  ) +
  theme_minimal()

# Save the plot as a PNG file
ggsave(paste0(sample1, "_", sample2, "_replicate_analysis.png"), plot = plot, width = 10, height = 8)

# Save the plot as a PDF file
ggsave(paste0(sample1, "_", sample2, "_replicate_analysis.pdf"), plot = plot, width = 10, height = 8)
