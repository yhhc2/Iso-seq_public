#!/bin/bash
# Hank Cheng
# 2024/05

# Runs faster if using more cores: salloc -A stergachislab --mem=500G --time="1-0" -p compute-ultramem -c 20

# Usage bash Complete_workflow_analysis_10.12.24.sh > Complete_workflow_analysis_10.12.24.output.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate

conda activate r_env_per_isoform


#########################################################
# Define inputs
#########################################################

# Scripts:
Split_read_stats_awk="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/1.Split_read_stats/Split_read_stats_awk.sh"
Cyclo_noncyclo_comparison="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/2.Collapse_counts/Cyclo_noncyclo_comparison.sh"
Isoform_analysis_p_value_generation="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/1.Isoform/Isoform_analysis_test_statistic_generation.R"
Group_by_abundance_in_noncyclo_genome_wide="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/3.Novel_iso_abundance_in_gene/Group_by_abundance_in_noncyclo_genome_wide_No_pvalue_calc.R"
PCA="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/4.PCA/PCA.R"
helper_script="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/2.Collapse_counts/Cyclo_noncyclo_comparison_general_collapse_isoform.sh"
#Filter_for_plotting_script="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/4.Filter_for_plotting/Filter_results_for_plotting.R"

# Directories:
Split_read_stats="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/1.Split_read_stats"
Collapse_counts="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/2.Collapse_counts"
Compare_samples_Isoform="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/3.Compare_samples/1.Isoform"
Compare_samples_Gene="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/3.Compare_samples/2.Gene"
Compare_samples_Novel_iso_abundance_in_gene="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/3.Compare_samples/3.Novel_iso_abundance_in_gene"
Compare_samples_PCA="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/3.Compare_samples/4.PCA"
#Filter_for_plotting_dir="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/4.Filter_for_plotting"

# Input files:
read_stat_file="/mmfs1/gscratch/stergachislab/asedeno/data/Hank_MASseq/6-collapsed/collapsed.read_stat.txt"
classification_file="/mmfs1/gscratch/stergachislab/asedeno/data/Hank_MASseq/7-pigeon/pigeon_classification.txt"
sample_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/Inputs/Sample_names_and_bam_locations.tsv"
omim_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/Inputs/genemap2.txt"

# Additional QC stuff
QC_SRSF6="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/QC_Values/SRSF6_pileup.sh"
QC_diversity="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/QC_Values/Isoform_gene_diversity.sh"
QC_dir="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/4.QC_Values"


# #########################################################
# # 1.Split_read_stats
# #########################################################
# mkdir -p "$Split_read_stats"
# cd "$Split_read_stats"

# #bash "$Split_read_stats_awk" \
# #    /path/to/read_stat.txt \
# #    /path/to/sample_file.tsv \
# #    > output.txt \
# #    2>&1

# bash "$Split_read_stats_awk" \
#     "$read_stat_file" \
#     "$sample_file" \
#     > output.txt \
#     2>&1


# #########################################################
# # 2.Collapse_counts
# #########################################################
# mkdir -p "$Collapse_counts"
# cd "$Collapse_counts"


# #bash "$Cyclo_noncyclo_comparison" \
# #    /path/to/sample_file.tsv \
# #    /path/to/split_read_stats \
# #    /path/to/output_base_dir \
# #    /path/to/classification.txt \
# #    /path/to/helper_script.sh \
# #    > output.txt \
# #    2>&1

# bash "$Cyclo_noncyclo_comparison" \
#     "$sample_file" \
#     "$Split_read_stats" \
#     "${Collapse_counts}/Results" \
#     "$classification_file" \
#     "$helper_script" \
#     > output.txt \
#     2>&1


# #########################################################
# # 3.Compare_samples
# #########################################################

# # 1.Isoform
# mkdir -p "$Compare_samples_Isoform"
# cd "$Compare_samples_Isoform"

# #Rscript "$Isoform_analysis_p_value_generation" \
# #  <sample_info_path> \
# #  <omim_file_path> \
# #  <gene_level> \
# #  <significance_threshold> \
# #  <masking_threshold> \
# #  <count_threshold>

# Rscript "$Isoform_analysis_p_value_generation" \
#   "$sample_file" \
#   "$omim_file" \
#   FALSE \
#   0.01 \
#   0.01 \
#   6 \
#   > output.txt \
#   2>&1


# # 2.Gene
# mkdir -p "$Compare_samples_Gene"
# cd "$Compare_samples_Gene"

# #Rscript "$Isoform_analysis_p_value_generation" \
# #  <sample_info_path> \
# #  <omim_file_path> \
# #  <gene_level> \
# #  <significance_threshold> \
# #  <masking_threshold> \
# #  <count_threshold>

# Rscript "$Isoform_analysis_p_value_generation" \
#   "$sample_file" \
#   "$omim_file" \
#   TRUE \
#   0.01 \
#   0.01 \
#   6 \
#   > output.txt \
#   2>&1


# # 3.Novel_iso_abundance_in_gene
# mkdir -p "$Compare_samples_Novel_iso_abundance_in_gene"
# cd "$Compare_samples_Novel_iso_abundance_in_gene"

# #Rscript "$Group_by_abundance_in_noncyclo_genome_wide" \
# #  <dt_isoform_level_file> \
# #  <dt_gene_level_file> \
# #  <hyp1_sig_proportion_masking_threshold> \
# #  <bin_proportion> \
# #  <significance_threshold> \
# #  <masking_threshold>

# Rscript "$Group_by_abundance_in_noncyclo_genome_wide" \
#   "${Compare_samples_Isoform}/data_combined_full.csv" \
#   "${Compare_samples_Gene}/data_combined_full.csv" \
#   0.5 \
#   0.005 \
#   1 \
#   0.00000001 \
#   > output.txt \
#   2>&1



# #########################################################
# # 4.QC
# #########################################################


# # 4.PCA
# mkdir -p "$Compare_samples_PCA"
# cd "$Compare_samples_PCA"

# #Rscript "$PCA" \
# #  <file_path> \
# #  <count_threshold> \
# #  <plot_title>


# # Isoform-level
# Rscript "$PCA" \
#   "${Compare_samples_Isoform}/data_combined_full.csv" \
#   6 \
#   "PCA_isoform_level.pdf" \
#   > output_isoform.txt \
#   2>&1


# # Gene-level
# Rscript "$PCA" \
#   "${Compare_samples_Gene}/data_combined_full.csv" \
#   6 \
#   "PCA_gene_level.pdf" \
#   > output_gene.txt \
#   2>&1



# # QC with SRSF6 cassette exon
# conda deactivate
# conda activate samtools
# mkdir -p "$QC_dir"
# cd "$QC_dir"
# bash "$QC_SRSF6" \
#   "$sample_file" \
#   "/gscratch/stergachislab/assemblies/simple-names/hg38.fa" \
#   "chr20" \
#   43459200 \
#   "QC_SRSF6_results.tsv" \
#   > output_QC_SRSF6.txt \
#   2>&1
# # QC with isoform and gene diversity
# bash "$QC_diversity" \
#   "${Compare_samples_Isoform}/data_combined_full.csv" \
#   unique_isoforms_counts.tsv \
#   > output_QC_diversity_isoform.txt \
#   2>&1
# bash "$QC_diversity" \
#   "${Compare_samples_Gene}/data_combined_full.csv" \
#   unique_genes_counts.tsv \
#   > output_QC_diversity_gene.txt \
#   2>&1


#########################################################
# 5.Test_statistics
#########################################################

Test_statistics="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/5.Test_statistics"
# conda deactivate
# conda activate python_env_isoseq

# mkdir -p "$Test_statistics"
# cd "$Test_statistics"

# # Define variables for input scripts
# HYP1_SCRIPT="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/5.Test_statistics/Test_statistic_Hyp1.ipynb"
# HYP2_LOE_SCRIPT="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/5.Test_statistics/Test_statistic_Hyp2_LOE.ipynb"
# HYP2_GOE_SCRIPT="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/5.Test_statistics/Test_statistic_Hyp2_GOE.ipynb"
# HYP3_GOE_SCRIPT="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/5.Test_statistics/Test_statistic_Hyp3_GOE.ipynb"
# HYP5_SCRIPT="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/5.Test_statistics/Test_statistic_Hyp5.ipynb"

# # Define variables for input files
# GENE_FILE="${Compare_samples_Novel_iso_abundance_in_gene}/data_combined_full_gene_with_Hyp5.csv"
# ISOFORM_FILE="${Compare_samples_Isoform}/data_combined_full.csv"

# # Hypothesis 1
# echo "Starting hypothesis 1"
# mkdir -p output/Hypothesis1/Gene
# mkdir -p output/Hypothesis1/Isoform
# cd output/Hypothesis1/Gene
# papermill $HYP1_SCRIPT output_Test_statistic_Hyp1_Gene.ipynb -p filename $GENE_FILE &
# cd ../Isoform
# papermill $HYP1_SCRIPT output_Test_statistic_Hyp1_Isoform.ipynb -p filename $ISOFORM_FILE &
# cd ../../..
# wait # Ensure all background processes for Hypothesis 1 finish
# echo "Completed hypothesis 1"

# # Hypothesis 2 LOE
# echo "Starting hypothesis 2 LOE"
# mkdir -p output/Hypothesis2LOE/Gene
# mkdir -p output/Hypothesis2LOE/Isoform
# cd output/Hypothesis2LOE/Gene
# papermill $HYP2_LOE_SCRIPT output_Test_statistic_Hyp2_LOE_Gene.ipynb -p filename $GENE_FILE &
# cd ../Isoform
# papermill $HYP2_LOE_SCRIPT output_Test_statistic_Hyp2_LOE_Isoform.ipynb -p filename $ISOFORM_FILE &
# cd ../../..
# wait # Ensure all background processes for Hypothesis 2 LOE finish
# echo "Completed hypothesis 2 LOE"

# # Hypothesis 2 GOE
# echo "Starting hypothesis 2 GOE"
# mkdir -p output/Hypothesis2GOE/Gene
# mkdir -p output/Hypothesis2GOE/Isoform
# cd output/Hypothesis2GOE/Gene
# papermill $HYP2_GOE_SCRIPT output_Test_statistic_Hyp2_GOE_Gene.ipynb -p filename $GENE_FILE &
# cd ../Isoform
# papermill $HYP2_GOE_SCRIPT output_Test_statistic_Hyp2_GOE_Isoform.ipynb -p filename $ISOFORM_FILE &
# cd ../../..
# wait # Ensure all background processes for Hypothesis 2 GOE finish
# echo "Completed hypothesis 2 GOE"

# # Hypothesis 3 GOE
# echo "Starting hypothesis 3 GOE"
# mkdir -p output/Hypothesis3GOE/Gene
# mkdir -p output/Hypothesis3GOE/Isoform
# cd output/Hypothesis3GOE/Gene
# papermill $HYP3_GOE_SCRIPT output_Test_statistic_Hyp3_GOE_Gene.ipynb -p filename $GENE_FILE &
# cd ../Isoform
# papermill $HYP3_GOE_SCRIPT output_Test_statistic_Hyp3_GOE_Isoform.ipynb -p filename $ISOFORM_FILE &
# cd ../../..
# wait # Ensure all background processes for Hypothesis 3 GOE finish
# echo "Completed hypothesis 3 GOE"

# # Hypothesis 5
# echo "Starting hypothesis 5"
# mkdir -p output/Hypothesis5/Gene
# cd output/Hypothesis5/Gene
# papermill $HYP5_SCRIPT output_Test_statistic_Hyp5_Gene.ipynb -p filename $GENE_FILE &
# cd ../../..
# wait # Ensure all background processes for Hypothesis 5 finish
# echo "Completed hypothesis 5"

# echo "Completed all hypotheses"


# #########################################################
# # 5.Test_statistics filtering and HPO
# #########################################################

# cd "$Test_statistics"

# # Define the input notebook as a variable
# SHINY_INPUT_NOTEBOOK='/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/5.Test_statistics/Prep_for_shiny.ipynb'

# # Define variables for paths
# PROBANDS_FILE='/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/Inputs/7.12.24_manualListOfSamples.tsv'
# HPO_FILE='/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/Inputs/phenotype.hpoa'
# SEQR_FILES_DIR='/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/Inputs/Seqr'

# papermill "$SHINY_INPUT_NOTEBOOK" Prep_for_shiny_output_notebook.ipynb \
#     -p input_file_hyp1_gene 'output/Hypothesis1/Gene/Test_statistic_hyp1_notebook_variables.pkl' \
#     -p input_file_hyp1_isoform 'output/Hypothesis1/Isoform/Test_statistic_hyp1_notebook_variables.pkl' \
#     -p input_file_hyp2loe_gene 'output/Hypothesis2LOE/Gene/Test_statistic_hyp2LOE_notebook_variables.pkl' \
#     -p input_file_hyp2loe_isoform 'output/Hypothesis2LOE/Isoform/Test_statistic_hyp2LOE_notebook_variables.pkl' \
#     -p input_file_hyp2goe_gene 'output/Hypothesis2GOE/Gene/Test_statistic_hyp2GOE_notebook_variables.pkl' \
#     -p input_file_hyp2goe_isoform 'output/Hypothesis2GOE/Isoform/Test_statistic_hyp2GOE_notebook_variables.pkl' \
#     -p input_file_hyp3goe_gene 'output/Hypothesis3GOE/Gene/Test_statistic_hyp3GOE_notebook_variables.pkl' \
#     -p input_file_hyp3goe_isoform 'output/Hypothesis3GOE/Isoform/Test_statistic_hyp3GOE_notebook_variables.pkl' \
#     -p input_file_hyp5_gene 'output/Hypothesis5/Gene/Test_statistic_hyp5_notebook_variables.pkl' \
#     -p probands_file "$PROBANDS_FILE" \
#     -p hpo_file "$HPO_FILE" \
#     -p genemap_file "$omim_file" \
#     -p seqr_files_dir "$SEQR_FILES_DIR"

#########################################################
# 6.Intersect_genetics
#########################################################

conda deactivate
conda activate vep_environment

VCF_dir="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/9.19.24/6.Intersect_genetics/VCF"
mkdir -p "$VCF_dir"
cd "$VCF_dir"

vep_script="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/6.Intersect_genetics/vep_annotate_and_filter_conda_gnomad_varinfo.sh"

bash "$vep_script" "$sample_file" . > vep_annotate_and_filter_conda_gnomad_varinfo_out.txt 2>&1

conda deactivate
conda activate python_env_isoseq

Concatenated_and_intersected_dir="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/6.Intersect_genetics/Concatenated_and_intersected"
mkdir -p "$Concatenated_and_intersected_dir"
cd "$Concatenated_and_intersected_dir"

concat_and_intersect_notebook="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/6.Intersect_genetics/Combine_csvs_and_filter_with_var_info.ipynb"

papermill "$concat_and_intersect_notebook" concat_and_intersect_output_notebook.ipynb \
    -p csv_files_gene '["added_HPO_and_seqr_Hyp1_Gene.csv", "added_HPO_and_seqr_Hyp2LOE_Gene.csv", "added_HPO_and_seqr_Hyp2GOE_Gene.csv", "added_HPO_and_seqr_Hyp3GOE_Gene.csv", "added_HPO_and_seqr_Hyp5_Gene.csv"]' \
    -p csv_files_isoform '["added_HPO_and_seqr_Hyp1_Isoform.csv", "added_HPO_and_seqr_Hyp2LOE_Isoform.csv", "added_HPO_and_seqr_Hyp2GOE_Isoform.csv", "added_HPO_and_seqr_Hyp3GOE_Isoform.csv"]' \
    -p tsv_file "$sample_file" \
    -p vcf_directory "$VCF_dir" \
    -p base_dir "$Test_statistics"



# Todo
# Add HPO for new patients. Get SpliceAI to work in VEP. 