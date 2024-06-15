# Hank Cheng
# yhhc@uw.edu

# Usage: 
# conda activate r_env_per_isoform
# Rscript 1.PCA.R "/path/to/data_combined_full_gene_with_Hyp5.csv" 6 "PCA_Gene_Level.pdf"

#setwd("~/2022-2023/Research/Computational/Isoform counting/Per_isoform_analysis/4.28.24_Comparison/Isoseq_comparison_cyclo_vs_noncyclo_4.28.24/PCA")

# Load the data.table package
library(data.table)
library(ggplot2)

# Parsing Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if (length(args) != 2) {
  stop("Usage: Rscript Replicate_analysis.R <file_path> <count_threshold> <plot_title>", call. = FALSE)
}

# Assign variables from arguments
file_path <- args[1]
count_threshold <- as.numeric(args[2])
plot_title <- args[3]


# Read the CSV file into a data table
dt <- fread(file_path)

# Print the first few rows of the data table
print(head(dt))

# To reduce dimensions, remove low coverage isoforms/genes
isoforms_to_keep <- dt[, .(
  keep = any(cyclo_count >= count_threshold | noncyclo_count >= count_threshold)
), by = Isoform_PBid][keep == TRUE, Isoform_PBid]

# Print helpful statements
print("Number of isoforms to keep after filtering based on counts:")
print(length(isoforms_to_keep))

# Filter out the data_with_totals to keep only the desired Isoform_PBid
dt <- dt[Isoform_PBid %in% isoforms_to_keep]

################################################
# Generate expression matrix
################################################

# Assume 'dt' is already your loaded data table as per the previous steps
# Now, reshape the data table to create the expression matrix
dt_long <- melt(dt, id.vars = c("Isoform_PBid", "Sample"), measure.vars = c("Cyclo_TPM", "Noncyclo_TPM"))

# Create a new identifier for sample type
dt_long[, variable := paste0(Sample, '_', sub("TPM", "", variable))]

# Cast the long format back to a wide format where each gene has its own row and each sample type its own column
expression_matrix <- dcast(dt_long, Isoform_PBid ~ variable, value.var = "value")

# Print the first few rows of the expression matrix to check
print(head(expression_matrix))

#write.csv(expression_matrix, "expression_matrix_ge10reads.csv", row.names = FALSE)  # Adjust row.names based on your needs

################################################
# Replicate analysis for UDN052264
################################################

# Make a scatter plot
# If the column contains "UDN052264_" and "_Noncyclo_" use this as the x-axis values
# If the column contains "UDN052264b_" and "_Noncyclo_" use this as the y-axis values

# Extract relevant columns for the scatter plot
x_col <- grep("UDN052264_Noncyclo_", colnames(expression_matrix), value = TRUE)
y_col <- grep("UDN052264b_Noncyclo_", colnames(expression_matrix), value = TRUE)

# Check if both columns are found
if (length(x_col) == 0 | length(y_col) == 0) {
  stop("Required columns for scatter plot not found in the expression matrix.")
}

# Calculate correlation coefficient
correlation <- cor(expression_matrix[[x_col]], expression_matrix[[y_col]], use = "complete.obs")

# Create scatter plot
plot <- ggplot(expression_matrix, aes_string(x = x_col, y = y_col)) +
  geom_point() +
  labs(
    title = paste0(plot_title, " (Correlation: ", round(correlation, 2), ")"),
    x = paste0("Expression of ", x_col),
    y = paste0("Expression of ", y_col)
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 5000)) +
  scale_y_continuous(limits = c(0, 5000))

# Save the plot as a PNG file
ggsave("UDN052264_replicate_analysis.png", plot = plot, width = 10, height = 8)

# Save the plot as a PDF file
ggsave("UDN052264_replicate_analysis.pdf", plot = plot, width = 10, height = 8)

################################################
# Replicate analysis for UDN212054
################################################

################################################
# Replicate analysis for UDN687128
################################################