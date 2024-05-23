
##################################################
# Volcano plot to visualize p-values
##################################################

##################################################
# Hypothesis 1
##################################################

##################################################
# Make volcano plot. Hypothesis 1.
##################################################

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = `NormalizedFractionDifference`, y = -log10(Hypothesis_1_P_values))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "CycloFraction-NoncycloFraction", y = "-log10(hypothesis 1 p-value)", title = paste(sample_of_interest)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest,".pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)

##################################################
# Make volcano plot. Hypothesis 1. _removed_if_any_controls_sig_for_hyp1.pdf
##################################################

# If the isoform has a significant hypothesis 1 p-value in any of the controls, 
# remove it.

isoforms_significant_in_other_samples <- subset(results_hypothesis_1, Sample != sample_of_interest & P_Value < p_value_thresh)

isoforms_significant_in_other_samples <- unique(isoforms_significant_in_other_samples$Isoform)

# Only include rows from merged_data_omim_added if the Isoform column doesn't contain a value that is present in the isoforms_significant_in_samples_other_than_UDN212054 character array
merged_data_omim_added_removed_common_hypothesis_1_isoforms <- subset(merged_data_omim_added, !(Isoform %in% isoforms_significant_in_other_samples))

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- merged_data_omim_added_removed_common_hypothesis_1_isoforms

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = `NormalizedFractionDifference`, y = -log10(Hypothesis_1_P_values))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "CycloFraction-NoncycloFraction", y = "-log10(hypothesis 1 p-value)", title = paste(sample_of_interest, " removed isoforms if any controls significant for hypothesis 1")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
#print(volcano_plot_with_label)  
title <- paste(sample_of_interest,"_removed_if_any_controls_sig_for_hyp1.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)

##################################################
# Hypothesis 2
##################################################


##################################################
# Hypothesis 2. _largestP_hypothesis_2_TPM_Norm.pdf
##################################################

data <- merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 Worst P-Value")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_TPM_Norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 2. _largestP_hypothesis_2_Normalized_Frac.pdf
##################################################

data <- merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedCycloFrac - NormalizedNoncycloFrac", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 worst P-Value")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_Normalized_Frac.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)



##################################################
# Hypothesis 2. _largestP_hypothesis_2_keepHyp1OnlySig_TPM_Norm.pdf
##################################################

# Keep isoforms that have significant hyp1 p-values ONLY in the sample of interest.

data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 worst, isoforms only significant for hyp1 in sample of interest")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_keepHyp1OnlySig_TPM_Norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 2. _largestP_hypothesis_2_keepHyp1OnlySig_Normalized_Frac.pdf
##################################################

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedCycloFrac - NormalizedNoncycloFrac", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 Fisher, removed isoforms sig in controls")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_keepHyp1OnlySig_Normalized_Frac.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 3
##################################################


##################################################
# Hypothesis 3. _largestP_hypothesis_3_TPM_Norm.pdf
##################################################

data <- merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM, y = -log10(Combined_Hypothesis_3_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Cyclo_TPM - Cyclo_TPM)/Cyclo_TPM", y = "-log10(Combined_Hypothesis_3_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 3 Worst P-Value")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_3_TPM_Norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 3. _largestP_hypothesis_3_Normalized_Frac.pdf
##################################################

data <- merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = -log10(Combined_Hypothesis_3_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedCycloFrac - NormalizedNoncycloFrac", y = "-log10(Combined_Hypothesis_3_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 3 worst P-Value")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
#print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_3_Normalized_Frac.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)



##################################################
# Hypothesis 3. _largestP_hypothesis_3_keepHyp1OnlySig_TPM_Norm.pdf
##################################################

# Keep isoforms that have significant hyp1 p-values ONLY in the sample of interest.

data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM, y = -log10(Combined_Hypothesis_3_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM", y = "-log10(Combined_Hypothesis_3_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 3 worst, isoforms only significant for hyp1 in sample of interest")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_3_keepHyp1OnlySig_TPM_Norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 3. _largestP_hypothesis_3_keepHyp1OnlySig_Normalized_Frac.pdf
##################################################

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = -log10(Combined_Hypothesis_3_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedCycloFrac - NormalizedNoncycloFrac", y = "-log10(Combined_Hypothesis_3_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 3 Fisher, removed isoforms sig in controls")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_3_keepHyp1OnlySig_Normalized_Frac.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)



##################################################
# Hypothesis 3. _Cyclo_TPM_Zscore_TPM_Norm.pdf
##################################################

data <- merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM, y = Cyclo_Z_Score)) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM", y = "Cyclo_Z_Score", title = paste(sample_of_interest, "Cyclo Z-score")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_Cyclo_TPM_Zscore_TPM_Norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 3. _Cyclo_TPM_Zscore_Normalized_Frac.pdf
##################################################

data <- merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = Cyclo_Z_Score)) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedFractionDifference", y = "Cyclo_Z_Score", title = paste(sample_of_interest, "Cyclo Z-score")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_Cyclo_TPM_Zscore_Normalized_Frac.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 3. _Cyclo_TPM_Zscore_keepHyp1OnlySig__TPM_Norm.pdf
##################################################

data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM, y = Cyclo_Z_Score)) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM", y = "Cyclo_Z_Score", title = paste(sample_of_interest, "Cyclo Z-score")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_Cyclo_TPM_Zscore_keepHyp1OnlySig__TPM_Norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 3. _Cyclo_TPM_Zscore__keepHyp1OnlySig_Normalized_Frac.pdf
##################################################

data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = Cyclo_Z_Score)) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedFractionDifference", y = "Cyclo_Z_Score", title = paste(sample_of_interest, "Cyclo Z-score")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_Cyclo_TPM_Zscore__keepHyp1OnlySig_Normalized_Frac.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


# LEFT OFF HERE need to add hypothesis 2 gene level


##################################################
# Gene level plots
##################################################

##################################################
# Hypothesis 2. _largestP_hypothesis_2_TPM_Norm_perGene.pdf
##################################################

data <- collapsed_by_gene_merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 worst p-value")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_TPM_Norm_perGene.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 2. _largestP_hypothesis_2_Normalized_Frac_perGene.pdf
##################################################

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- collapsed_by_gene_merged_data_omim_added

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedCycloFrac - NormalizedNoncycloFrac", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 worst p-value")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_Normalized_Frac_perGene.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 2. _largestP_hypothesis_2_keepHyp1OnlySig_TPM_Norm_perGene.pdf
##################################################

# LEFT OFF HERE. Will need to do hypothesis 1 p-values per gene. 

# If the isoform has a significant hypothesis 1 p-value in any of the controls, 
# remove it.
isoforms_significant_in_other_samples <- subset(collapsed_by_gene_merged_data_omim_added, Sample != sample_of_interest & P_Value < p_value_thresh)
isoforms_significant_in_other_samples <- unique(isoforms_significant_in_other_samples$Isoform)
# Only include rows from merged_data_omim_added if the Isoform column doesn't contain a value that is present in the isoforms_significant_in_samples_other_than_UDN212054 character array
merged_data_omim_added_removed_common_hypothesis_1_isoforms <- subset(merged_data_omim_added, !(Isoform %in% isoforms_significant_in_other_samples))



# Keep isoforms that have significant hyp1 p-values ONLY in the sample of interest.

data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 worst, isoforms only significant for hyp1 in sample of interest")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_keepHyp1OnlySig_TPM_Norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Hypothesis 2. _largestP_hypothesis_2_keepHyp1OnlySig_Normalized_Frac.pdf
##################################################

# LEFT OFF HERE. Will need to do hypothesis 1 p-values per gene. 

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = NormalizedFractionDifference, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedCycloFrac - NormalizedNoncycloFrac", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 Fisher, removed isoforms sig in controls")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_keepHyp1OnlySig_Normalized_Frac.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


##################################################
# Plot if isoforms are significant across all comparison with sample of interest. 
##################################################

##################################################
# Hypothesis 2 per gene
##################################################

##################################################
# Hypothesis 2. _largestP_hypothesis_2_TPM_Norm_perGene_hyp2_sig_for_all_comparisons.pdf
##################################################


# Subsetting the dataframe
filtered_df_for_sample_of_interest <- results_hypothesis_2_per_gene %>%
  filter(Sample1 == sample_of_interest | Sample2 == sample_of_interest)

filtered_df_for_sample_of_interest_all_sig_comparisons <- filtered_df_for_sample_of_interest %>%
  group_by(Isoform) %>%
  filter(all(P_Value < p_value_thresh)) %>%
  ungroup()

isoforms_sig_across_all_comparisons <- unique(filtered_df_for_sample_of_interest_all_sig_comparisons$Isoform)

data <- subset(collapsed_by_gene_merged_data_omim_added, (Isoform %in% isoforms_sig_across_all_comparisons))

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM, y = -log10(Combined_Hypothesis_2_Worst_P_Value))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(Noncyclo_TPM - Avg_Noncyclo_TPM)/Avg_Noncyclo_TPM", y = "-log10(Combined_Hypothesis_2_Worst_P_Value)", title = paste(sample_of_interest, "Combined hypothesis 2 worst p-value")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_largestP_hypothesis_2_TPM_Norm_perGene_hyp2_sig_for_all_comparisons.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)

##################################################
# Plotting isoform vs gene
##################################################

##################################################
# _isoform_vs_gene_Cyclo_TPM_norm.pdf
##################################################


data <- subset(merged_data_omim_added_removed_common_hypothesis_1_isoforms, Hypothesis_1_P_values < p_value_thresh)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM, y = Gene_Cyclo_TPM_norm)) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene == gene_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "Isoform (Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM", y = "Gene (Cyclo_TPM - Avg_Cyclo_TPM)/Avg_Cyclo_TPM", title = paste(sample_of_interest, "Cyclo Z-score")) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

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
  )

# Print the plot
# print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_isoform_vs_gene_Cyclo_TPM_norm.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)


