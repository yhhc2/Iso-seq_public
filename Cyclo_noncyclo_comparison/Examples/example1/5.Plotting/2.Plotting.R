library(ggplot2)
library(data.table)

data_for_plotting <- readRDS("~/Research/Projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/5.Plotting/data_for_plotting.rds")

to_save_into_rds <- c()

##################################################
# HARS1 Hypothesis 1
##################################################

significance_thresh <- 0.05

masking_thresh <- 0.01

sample_of_interest <- "UDN212054"
genes_of_interest <- c("HARS1", "EDN1")

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- data_for_plotting$isoforms_filtered_hyp1

isoforms_significant_in_other_samples <- subset(data, Sample != sample_of_interest & P_Value_Hyp1 < masking_thresh)

isoforms_significant_in_other_samples <- unique(isoforms_significant_in_other_samples$Isoform_PBid)

# Only include rows from merged_data_omim_added if the Isoform column doesn't contain a value that is present in the isoforms_significant_in_samples_other_than_UDN212054 character array
data <- subset(data, !(Isoform_PBid %in% isoforms_significant_in_other_samples))

data <- subset(data, Sample == sample_of_interest)

data <- subset(data, P_Value_Hyp1 < significance_thresh)

data$PhenotypesNotEmpty <- data$PhenotypesNotEmpty <- data$Phenotypes != "" & !is.na(data$Phenotypes)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = `NormalizedFractionDifference`, y = -log10(P_Value_Hyp1))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene %in% genes_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "NormalizedFractionDifference", y = "-log10(P_Value_Hyp1)", title = paste(sample_of_interest)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

# Add the label for the gene of interest using geom_text()
# The hjust and vjust arguments should be negative to position the label correctly above the point
volcano_plot_with_label <- volcano_plot +
  geom_text(
    data = subset(data, associated_gene %in% genes_of_interest),
    aes(label = associated_gene),
    hjust = 0.5, vjust = 1.5,
    size = 2,
    color = "Green",
    fontface = "bold"
  )

# Print the plot
#print(volcano_plot_with_label)
title <- paste(sample_of_interest,".pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)
ggsave(filename = "Hyp1.png", plot = volcano_plot_with_label, width = 8, height = 6, dpi = 300)


data <- subset(data, associated_gene %in% genes_of_interest)
# Save the data as a dataframe.
to_save_into_rds$hyp1 <- data

##################################################
# VTA1 Hypothesis 2 gene below median
##################################################

sample_of_interest <- "UDN687128"
genes_of_interest <- c("VTA1")

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- data_for_plotting$genes_filtered_hyp2_below_median

data <- subset(data, Sample == sample_of_interest)

data <- subset(data, Max_P_Value_Hyp2_below_median < significance_thresh)

data$PhenotypesNotEmpty <- data$PhenotypesNotEmpty <- data$Phenotypes != "" & !is.na(data$Phenotypes)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = `Noncyclo_Z_Score`, y = -log10(Max_P_Value_Hyp2_below_median))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene %in% genes_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "Noncyclo_Z_Score", y = "-log10(Max_P_Value_Hyp2_below_median)", title = paste(sample_of_interest)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

# Add the label for the gene of interest using geom_text()
# The hjust and vjust arguments should be negative to position the label correctly above the point
volcano_plot_with_label <- volcano_plot +
  geom_text(
    data = subset(data, associated_gene %in% genes_of_interest),
    aes(label = associated_gene),
    hjust = 0.5, vjust = 1.5,
    size = 2,
    color = "Green",
    fontface = "bold"
  )

# Print the plot
#print(volcano_plot_with_label)
title <- paste(sample_of_interest,".pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)
ggsave(filename = "Hyp2_gene_below.png", plot = volcano_plot_with_label, width = 8, height = 6, dpi = 300)


data <- subset(data, associated_gene %in% genes_of_interest)
# Save the data as a dataframe.
to_save_into_rds$hyp2_gene_below <- data


##################################################
# SET Hypothesis 3 isoform above median
##################################################


sample_of_interest <- "UDN215640"
genes_of_interest <- c("SET")

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- data_for_plotting$isoforms_filtered_hyp3_above_median

data <- subset(data, Sample == sample_of_interest)

data <- subset(data, Max_P_Value_Hyp3_above_median < significance_thresh)

data$PhenotypesNotEmpty <- data$PhenotypesNotEmpty <- data$Phenotypes != "" & !is.na(data$Phenotypes)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = `Cyclo_Z_Score`, y = -log10(Max_P_Value_Hyp3_above_median))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene %in% genes_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "Cyclo_Z_Score", y = "-log10(Max_P_Value_Hyp3_above_median)", title = paste(sample_of_interest)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

# Add the label for the gene of interest using geom_text()
# The hjust and vjust arguments should be negative to position the label correctly above the point
volcano_plot_with_label <- volcano_plot +
  geom_text(
    data = subset(data, associated_gene %in% genes_of_interest),
    aes(label = associated_gene),
    hjust = 0.5, vjust = 1.5,
    size = 2,
    color = "Green",
    fontface = "bold"
  )

# Print the plot
#print(volcano_plot_with_label)
title <- paste(sample_of_interest,".pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)
ggsave(filename = "Hyp3_isoform_above.png", plot = volcano_plot_with_label, width = 8, height = 6, dpi = 300)


data <- subset(data, associated_gene %in% genes_of_interest)
# Save the data as a dataframe.
to_save_into_rds$hyp3_isoform_above <- data


##################################################
# MFN2 Hypothesis 5
##################################################

masking_thresh <- 0.0000001
hyp1_sig_proportion_masking_threshold <- 0.5

sample_of_interest <- "UDN633333"
genes_of_interest <- c("MFN2")

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- data_for_plotting$genes_filtered_hyp5

isoforms_significant_in_other_samples <- subset(data, Sample != sample_of_interest & P_Value_Hyp5 < masking_thresh)
isoforms_significant_in_other_samples <- unique(isoforms_significant_in_other_samples$associated_gene)
# Only include rows from merged_data_omim_added if the Isoform column doesn't contain a value that is present in the isoforms_significant_in_samples_other_than_UDN212054 character array
data <- subset(data, !(associated_gene %in% isoforms_significant_in_other_samples))

data_data_table <- setDT(data_for_plotting$genes_filtered_hyp5)
control_data <- data_data_table[Sample != sample_of_interest]
# Check the proportion of controls where P_Values_Hyp1 for a gene are below the masking threshold and NormalizedFractionDifference > 0
gene_significant_in_controls_hyp1 <- control_data[, .(
  ProportionSignificant = sum(P_Value_Hyp1 < masking_thresh & NormalizedFractionDifference > 0) / .N
), by = .(associated_gene)]
# Filter to keep only those genes where ProportionSignificant is greater than or equal to 0.5
significant_genes <- gene_significant_in_controls_hyp1[ProportionSignificant >= hyp1_sig_proportion_masking_threshold, associated_gene]
# Remove these genes from data
data <- subset(data, !(associated_gene %in% significant_genes))

data <- subset(data, Sample == sample_of_interest)

data <- subset(data, P_Value_Hyp5 < significance_thresh)

data$PhenotypesNotEmpty <- data$PhenotypesNotEmpty <- data$Phenotypes != "" & !is.na(data$Phenotypes)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = (proportion_in_Bin1_cyclo - proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo), y = -log10(P_Value_Hyp5))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene %in% genes_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "(proportion_in_Bin1_cyclo − proportion_in_Bin1_noncyclo) / (proportion_in_Bin1_cyclo + proportion_in_Bin1_noncyclo)", y = "-log10(P_Value_Hyp5)", title = paste(sample_of_interest)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size = 7))  

# Add the label for the gene of interest using geom_text()
# The hjust and vjust arguments should be negative to position the label correctly above the point
volcano_plot_with_label <- volcano_plot +
  geom_text(
    data = subset(data, associated_gene %in% genes_of_interest),
    aes(label = associated_gene),
    hjust = 0.5, vjust = 1.5,
    size = 2,
    color = "Green",
    fontface = "bold"
  )

# Print the plot
#print(volcano_plot_with_label)
title <- paste(sample_of_interest,".pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)
ggsave(filename = "Hyp5.png", plot = volcano_plot_with_label, width = 8, height = 6, dpi = 300)


data <- subset(data, associated_gene %in% genes_of_interest)
# Save the data as a dataframe.
to_save_into_rds$hyp5 <- data


##################################################
# SET Hypothesis 2 isoform above median
##################################################


sample_of_interest <- "UDN215640"
genes_of_interest <- c("SET")

# Replace 'gene_of_interest' with the name of the gene you want to label in the plot
data <- data_for_plotting$isoforms_filtered_hyp2_above_median

data <- subset(data, Sample == sample_of_interest)

data <- subset(data, Max_P_Value_Hyp2_above_median < significance_thresh)

data$PhenotypesNotEmpty <- data$PhenotypesNotEmpty <- data$Phenotypes != "" & !is.na(data$Phenotypes)

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(data, aes(x = `Noncyclo_Z_Score`, y = -log10(Max_P_Value_Hyp2_above_median))) +
  geom_point(aes(color = PhenotypesNotEmpty), alpha = 0.6, size = 1) +  # Gray points with colors based on phenotypes
  geom_point(data = subset(data, associated_gene %in% genes_of_interest), color = "blue", size = 3) +  # Highlighted gene point
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Add significance threshold line
  labs(x = "Noncyclo_Z_Score", y = "-log10(Max_P_Value_Hyp2_above_median)", title = paste(sample_of_interest)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), labels = c("No Phenotype", "Phenotype")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Hide the legend title (optional)

# Add the label for the gene of interest using geom_text()
# The hjust and vjust arguments should be negative to position the label correctly above the point
volcano_plot_with_label <- volcano_plot +
  geom_text(
    data = subset(data, associated_gene %in% genes_of_interest),
    aes(label = associated_gene),
    hjust = 0.5, vjust = 1.5,
    size = 2,
    color = "Green",
    fontface = "bold"
  )

# Print the plot
#print(volcano_plot_with_label)
title <- paste(sample_of_interest,"_Hyp2.pdf")
ggsave(title, plot = volcano_plot_with_label, width = 8, height = 6)
ggsave(filename = "SET_Hyp2_isoform_above.png", plot = volcano_plot_with_label, width = 8, height = 6, dpi = 300)


data <- subset(data, associated_gene %in% genes_of_interest)
# Save the data as a dataframe.
to_save_into_rds$SET_Hyp2_isoform_above <- data

saveRDS(to_save_into_rds, "markdown_tables.rds")


