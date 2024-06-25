
library(data.table)
library(ggplot2)

################################################
# Generate PCA plot
################################################

gene_pca_data_from_csv <- fread("PCA_gene_outlier_analysis.csv")

# Create the PCA plot with labels and variance explained in the axis titles
gene_pca_plot_cyclo_noncyclo <- ggplot(gene_pca_data_from_csv, aes(x = PC1, y = PC2, color = Type, label = Sample)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_text(aes(label = Sample), vjust = -0.1, hjust = -0.1, size = 2, color = "black") +  # Adjust label position and size
  scale_color_manual(values = c("cyclo" = "red", "noncyclo" = "blue")) +
  labs(title = "PCA Plot",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Type") +
  theme_minimal()

gene_pca_plot_Sex <- ggplot(gene_pca_data_from_csv, aes(x = PC1, y = PC2, color = Sex, label = Sample)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_text(aes(label = Sample), vjust = -0.1, hjust = -0.1, size = 2, color = "black") +  # Adjust label position and size
  scale_color_manual(values = c("F" = "red", "M" = "blue")) +
  labs(title = "PCA Plot",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Sex") +
  theme_minimal()

gene_pca_plot_Project <- ggplot(gene_pca_data_from_csv, aes(x = PC1, y = PC2, color = Project, label = Sample)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_text(aes(label = Sample), vjust = -0.1, hjust = -0.1, size = 2, color = "black") +  # Adjust label position and size
  scale_color_manual(values = c("UDN" = "red", "BCH" = "blue")) +
  labs(title = "PCA Plot",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Project") +
  theme_minimal()


isoform_pca_data_from_csv <- fread("PCA_isoform_outlier_analysis.csv")

# Create the PCA plot with labels and variance explained in the axis titles
isoform_pca_plot_cyclo_noncyclo <- ggplot(isoform_pca_data_from_csv, aes(x = PC1, y = PC2, color = Type, label = Sample)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_text(aes(label = Sample), vjust = -0.1, hjust = -0.1, size = 2, color = "black") +  # Adjust label position and size
  scale_color_manual(values = c("cyclo" = "red", "noncyclo" = "blue")) +
  labs(title = "PCA Plot",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Type") +
  theme_minimal()

isoform_pca_plot_Sex <- ggplot(isoform_pca_data_from_csv, aes(x = PC1, y = PC2, color = Sex, label = Sample)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_text(aes(label = Sample), vjust = -0.1, hjust = -0.1, size = 2, color = "black") +  # Adjust label position and size
  scale_color_manual(values = c("F" = "red", "M" = "blue")) +
  labs(title = "PCA Plot",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Sex") +
  theme_minimal()

isoform_pca_plot_Project <- ggplot(isoform_pca_data_from_csv, aes(x = PC1, y = PC2, color = Project, label = Sample)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_text(aes(label = Sample), vjust = -0.1, hjust = -0.1, size = 2, color = "black") +  # Adjust label position and size
  scale_color_manual(values = c("UDN" = "red", "BCH" = "blue")) +
  labs(title = "PCA Plot",
       x = paste("PC1 (", variance_explained[1], "% variance)", sep = ""),
       y = paste("PC2 (", variance_explained[2], "% variance)", sep = ""),
       color = "Sample Project") +
  theme_minimal()


# Save the plots
ggsave("gene_pca_plot_cyclo_noncyclo.pdf", plot = gene_pca_plot_cyclo_noncyclo, width = 8, height = 6)
ggsave("gene_pca_plot_Sex.pdf", plot = gene_pca_plot_Sex, width = 8, height = 6)
ggsave("gene_pca_plot_Project.pdf", plot = gene_pca_plot_Project, width = 8, height = 6)

ggsave("isoform_pca_plot_cyclo_noncyclo.pdf", plot = isoform_pca_plot_cyclo_noncyclo, width = 8, height = 6)
ggsave("isoform_pca_plot_Sex.pdf", plot = isoform_pca_plot_Sex, width = 8, height = 6)
ggsave("isoform_pca_plot_Project.pdf", plot = isoform_pca_plot_Project, width = 8, height = 6)
