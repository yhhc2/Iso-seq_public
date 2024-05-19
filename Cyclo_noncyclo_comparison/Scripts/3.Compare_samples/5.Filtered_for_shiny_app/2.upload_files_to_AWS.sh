# !/bin/bash
# Hank Cheng

# sh upload_files_to_AWS.sh > upload_files_to_AWS_out.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate
conda activate ft-pipeline

S3URL=""


##########################################################################
# Upload each file to AWS
##########################################################################


aws s3 cp "/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/4.24.24_merge_aligned_bams/3.Comparison_between_samples/Filter_results_for_shiny/hank_data_hyp123.rds" $S3URL
aws s3 cp "/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/4.24.24_merge_aligned_bams/3.Comparison_between_samples/Gene/hank_data_hyp5.rds" $S3URL
