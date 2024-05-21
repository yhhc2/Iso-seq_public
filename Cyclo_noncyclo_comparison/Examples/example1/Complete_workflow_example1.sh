#!/bin/bash
# Hank Cheng
# 2024/05

# Usage bash Complete_workflow_example1.sh > Complete_workflow_example1.output.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate

conda activate r_env_per_isoform


#########################################################
# Define inputs
#########################################################

# Scripts:
Split_read_stats_awk="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/1.Split_read_stats/Split_read_stats_awk.sh"
Cyclo_noncyclo_comparison="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/2.Collapse_counts/Cyclo_noncyclo_comparison.sh"
Isoform_analysis_p_value_generation="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/1.Isoform/Isoform_analysis_p_value_generation.R"
Group_by_abundance_in_noncyclo_genome_wide="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/3.Novel_iso_abundance_in_gene"
PCA="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/4.PCA/PCA.R"

# Directories:
Split_read_stats="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/1.Split_read_stats"
Collapse_counts="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/2.Collapse_counts"
Compare_samples_Isoform="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/3.Compare_samples/1.Isoform"
Compare_samples_Gene="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/3.Compare_samples/2.Gene"
Compare_samples_Novel_iso_abundance_in_gene="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/3.Compare_samples/3.Novel_iso_abundance_in_gene"
Compare_samples_PCA="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/3.Compare_samples/4.PCA"

# Input files:
read_stat_file="/mmfs1/gscratch/stergachislab/asedeno/data/Hank_MASseq/6-collapsed/collapsed.read_stat.txt"
sample_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Examples/example1/Inputs/Sample_names_and_bam_locations.tsv"
classification_file="/mmfs1/gscratch/stergachislab/asedeno/data/Hank_MASseq/7-pigeon/pigeon_classification.txt"
helper_script="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/2.Collapse_counts/Cyclo_noncyclo_comparison_general_collapse_isoform.sh"
omim_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/3.15.24_merge_aligned_bams/4.Comparison_between_samples/Combined_OMIM_Isoform_4.3.24/genemap2_7.21.23.txt"

#########################################################
# 1.Split_read_stats
#########################################################
mkdir -p "$Split_read_stats"
cd "$Split_read_stats"

#bash "$Split_read_stats_awk" \
#    /path/to/read_stat.txt \
#    /path/to/sample_file.tsv \
#    > output.txt \
#    2>&1

bash "$Split_read_stats_awk" \
    "$read_stat_file" \
    "$sample_file" \
    > output.txt \
    2>&1


#########################################################
# 2.Collapse_counts
#########################################################
mkdir -p "$Collapse_counts"
cd "$Collapse_counts"


#bash "$Cyclo_noncyclo_comparison" \
#    /path/to/sample_file.tsv \
#    /path/to/split_read_stats \
#    /path/to/output_base_dir \
#    /path/to/classification.txt \
#    /path/to/helper_script.sh \
#    > output.txt \
#    2>&1

bash "$Cyclo_noncyclo_comparison" \
    "$sample_file" \
    "$Split_read_stats" \
    "${Collapse_counts}/Results" \
    "$classification_file" \
    "$helper_script" \
    > output.txt \
    2>&1


#########################################################
# 3.Compare_samples
#########################################################

# 1.Isoform
mkdir -p "$Compare_samples_Isoform"
cd "$Compare_samples_Isoform"

#Rscript "$Isoform_analysis_p_value_generation" \
#  <sample_info_path> \
#  <omim_file_path> \
#  <gene_level> \
#  <significance_threshold> \
#  <masking_threshold> \
#  <count_threshold>

Rscript "$Isoform_analysis_p_value_generation" \
  "$sample_file" \
  "$omim_file" \
  FALSE \
  0.01 \
  0.01 \
  6 \
  > output.txt \
  2>&1


# 2.Gene
mkdir -p "$Compare_samples_Gene"
cd "$Compare_samples_Gene"

#Rscript "$Isoform_analysis_p_value_generation" \
#  <sample_info_path> \
#  <omim_file_path> \
#  <gene_level> \
#  <significance_threshold> \
#  <masking_threshold> \
#  <count_threshold>

Rscript "$Isoform_analysis_p_value_generation" \
  "$sample_file" \
  "$omim_file" \
  TRUE \
  0.01 \
  0.01 \
  6 \
  > output.txt \
  2>&1



# 3.Novel_iso_abundance_in_gene
mkdir -p "$Compare_samples_Novel_iso_abundance_in_gene"
cd "$Compare_samples_Novel_iso_abundance_in_gene"

#Rscript "$Group_by_abundance_in_noncyclo_genome_wide" \
#  <dt_isoform_level_file> \
#  <dt_gene_level_file> \
#  <hyp1_sig_proportion_masking_threshold> \
#  <bin_proportion> \
#  <significance_threshold> \
#  <masking_threshold>

Rscript "$Group_by_abundance_in_noncyclo_genome_wide" \
  "${Compare_samples_Isoform}/data_combined_full.csv" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  0.5 \
  0.005 \
  1 \
  0.00000001 \
  > output.txt \
  2>&1



# 4.PCA
mkdir -p "$Compare_samples_PCA"
cd "$Compare_samples_PCA"

#Rscript "$PCA" \
#  <file_path> \
#  <count_threshold> \
#  <plot_title>


# Isoform-level
Rscript "$PCA" \
  "${Compare_samples_Isoform}/data_combined_full.csv" \
  6 \
  "PCA_isoform_level.pdf"

# Gene-level
Rscript "$PCA" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  6 \
  "PCA_gene_level.pdf"




