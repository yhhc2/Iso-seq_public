#!/bin/bash
# Hank Cheng
# 2024/05

# Usage bash Complete_workflow_analysis_6.14.24_added_QC.sh > Complete_workflow_analysis_6.14.24_added_QC.output.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate

#conda activate r_env_per_isoform


#########################################################
# Define inputs
#########################################################

# Scripts:
Split_read_stats_awk="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/1.Split_read_stats/Split_read_stats_awk.sh"
Cyclo_noncyclo_comparison="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/2.Collapse_counts/Cyclo_noncyclo_comparison.sh"
Isoform_analysis_p_value_generation="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/1.Isoform/Isoform_analysis_p_value_generation.R"
Group_by_abundance_in_noncyclo_genome_wide="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/3.Novel_iso_abundance_in_gene/Group_by_abundance_in_noncyclo_genome_wide.R"
PCA="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/4.PCA/PCA.R"
helper_script="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/2.Collapse_counts/Cyclo_noncyclo_comparison_general_collapse_isoform.sh"
Filter_for_plotting_script="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/4.Filter_for_plotting/Filter_results_for_plotting.R"

# Directories:
Split_read_stats="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/1.Split_read_stats"
Collapse_counts="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/2.Collapse_counts"
Compare_samples_Isoform="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/3.Compare_samples/1.Isoform"
Compare_samples_Gene="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/3.Compare_samples/2.Gene"
Compare_samples_Novel_iso_abundance_in_gene="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/3.Compare_samples/3.Novel_iso_abundance_in_gene"
Compare_samples_PCA="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/3.Compare_samples/4.PCA"
Filter_for_plotting_dir="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/4.Filter_for_plotting"

# Input files:
read_stat_file="/mmfs1/gscratch/stergachislab/asedeno/data/Hank_MASseq/6-collapsed/collapsed.read_stat.txt"
classification_file="/mmfs1/gscratch/stergachislab/asedeno/data/Hank_MASseq/7-pigeon/pigeon_classification.txt"
sample_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/Inputs/Sample_names_and_bam_locations.tsv"
omim_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/Inputs/genemap2_6.10.2024.txt"


# Additional QC stuff
QC_SRSF6="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/QC_Values/SRSF6_pileup.sh"
QC_SRSF6_dir="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/QC_Values"

#########################################################
# QC with SRSF6 cassette exon
#########################################################

conda activate samtools

mkdir -p "$QC_SRSF6_dir"
cd "$QC_SRSF6_dir"

bash "$QC_SRSF6" \
  "$sample_file" \
  "/gscratch/stergachislab/assemblies/simple-names/hg38.fa" \
  "chr20" \
  43459200 \
  "QC_SRSF6_results.tsv" \
  > output.txt \
  2>&1
