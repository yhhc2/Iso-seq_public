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
if (length(args) != 3) {
  stop("Usage: Rscript 1.PCA.R <file_path> <count_threshold> <plot_title>", call. = FALSE)
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
# Generate PCA plot
################################################

# Assuming 'expression_matrix' is already loaded and prepared
# Remove the gene names column for PCA analysis
data_for_pca <- data.table(expression_matrix[, -1, with = FALSE])

# Transpose the data to have samples as rows and genes as columns
data_for_pca <- t(data_for_pca)

# Perform PCA
pca_results <- prcomp(data_for_pca, center = TRUE, scale. = TRUE)

# Calculate percentage of variance explained by PCs.
variance_explained <- round(pca_results$sdev^2 / sum(pca_results$sdev^2) * 100, 1)

# Create a data frame for plotting
pca_data <- data.frame(Sample = rownames(pca_results$x), PC1 = pca_results$x[,1], PC2 = pca_results$x[,2])

# Add a column to distinguish between cyclo and noncyclo samples
pca_data$Type <- ifelse(grepl("_Cyclo", pca_data$Sample), "cyclo", "noncyclo")

# Create the PCA plot with labels and variance explained in the axis titles
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Type, label = Sample)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_text(aes(label = Sample), vjust = -0.1, hjust = -0.1, size = 2, color = "black") +  # Adjust label position and size
  scale_color_manual(values = c("cyclo" = "blue", "noncyclo" = "red")) +
  labs(title = "PCA Plot",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Type") +
  theme_minimal()

# Print data for inspection
print(pca_data)

# Save the plot
ggsave(plot_title, plot = pca_plot, width = 8, height = 6)
ggsave(paste0(plot_title,".png"), plot = pca_plot, width = 8, height = 6, dpi = 300)


################################################
# Random Forest for feature selection
################################################
# Load necessary libraries
library(randomForest)

# Assuming `expression_matrix` is your data frame with samples as rows and genes as columns
# And `treatment_status` is a vector indicating which samples were treated with cycloheximide (1 for treated, 0 for not)

# Combine the data into a single data frame
data <- as.data.frame(data_for_pca)
colnames(data) <- as.character(expression_matrix[, 1])

treatment_status <- ifelse(grepl("_Cyclo", rownames(data)), 1, 0)
data$treatment_status <- factor(treatment_status)

print(head(data))

# Split the data into features and target
features <- data[, -ncol(data)]  # All columns except the last one (treatment_status)
target <- data$treatment_status  # The last column (treatment_status)

# Train the random forest model
set.seed(123)  # For reproducibility
rf_model <- randomForest(x = features, y = target, importance = TRUE, ntree = 500)


# Print the model to see the error rates
print(rf_model)

# Extract and print the OOB error rate
oob_error_rate <- rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
cat("\nOOB Error Rate: ", oob_error_rate, "\n")

# Evaluate feature importance
importance_scores <- importance(rf_model)
importance_df <- data.frame(Gene = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])

# Sort by importance
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

# Display the top 10 most important genes
print(head(importance_df, 30))

# Create the plot and assign it to a variable
plot <- ggplot(importance_df[1:30, ], aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("Gene") +
  ylab("Importance") +
  ggtitle("Top 30 Most Important Genes for Predicting Cycloheximide Treatment")

# Save the plot as a PNG file
ggsave("gene_importance_plot.png", plot = plot, width = 10, height = 8)

# Save the plot as a PDF file
ggsave("gene_importance_plot.pdf", plot = plot, width = 10, height = 8)
