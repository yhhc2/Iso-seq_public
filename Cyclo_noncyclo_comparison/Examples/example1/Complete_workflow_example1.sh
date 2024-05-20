#!/bin/bash
# Hank Cheng
# 2024/05

#########################################################
# Define inputs
#########################################################

# Make sure the output files are written in the correct locations. 

#########################################################
# 1.Split_read_stats
#########################################################

#bash 1.Split_read_stats_awk.sh \
#    /path/to/read_stat.txt \
#    /path/to/sample_file.tsv \
#    > output.txt \
#    2>&1

#########################################################
# 2.Collapse_counts
#########################################################

#bash 1.Cyclo_noncyclo_comparison.sh \
#    /path/to/sample_file.tsv \
#    /path/to/split_read_stats \
#    /path/to/output_base_dir \
#    /path/to/classification.txt \
#    /path/to/helper_script.sh \
#    > output.txt \
#    2>&1

#########################################################
# 3.Compare_samples
#########################################################

# 1.Isoform
#Rscript 1.Isoform_analysis_p_value_generation.R \
#  <sample_info_path> \
#  <omim_file_path> \
#  <gene_level> \
#  <significance_threshold> \
#  <masking_threshold> \
#  <count_threshold>

# 2.Gene
#Rscript 1.Isoform_analysis_p_value_generation.R \
#  <sample_info_path> \
#  <omim_file_path> \
#  <gene_level> \
#  <significance_threshold> \
#  <masking_threshold> \
#  <count_threshold>

# 3.Novel_iso_abundance_in_gene
#Rscript 1.Group_by_abundance_in_noncyclo_genome_wide.R \
#  <dt_isoform_level_file> \
#  <dt_gene_level_file> \
#  <hyp1_sig_proportion_masking_threshold> \
#  <bin_proportion> \
#  <significance_threshold> \
#  <masking_threshold>

# 4.PCA
#Rscript 1.PCA.R \
#  <file_path> \
#  <count_threshold> \
#  <plot_title>


