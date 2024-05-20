# Usage: 
# conda activate r_env_per_isoform
# Rscript 1.PCA.R > 1.PCA_output.txt 2>&1

#setwd("~/2022-2023/Research/Computational/Isoform counting/Per_isoform_analysis/4.28.24_Comparison/Isoseq_comparison_cyclo_vs_noncyclo_4.28.24/PCA")

# Load the data.table package
library(data.table)
library(ggplot2)

count_threshold <- 6

# Specify the path to the CSV file
file_path <- "/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/3.Novel_iso_abundance_in_gene/data_combined_full_gene_with_Hyp5.csv"

# Read the CSV file into a data table
dt <- fread(file_path)

# Remove BCH_1199-01 because it's an outlier
# dt <- dt[Sample != "BCH_1199-01"]

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
  labs(title = "PCA Plot Gene-Level",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Type") +
  theme_minimal()

# Print data for inspection
print(pca_data)

# Save the plot
ggsave("PCA_Gene_level_with_Hyp5_data_combined.pdf", plot = pca_plot, width = 8, height = 6)