# !/bin/bash
# Hank Cheng

# sh upload_files_to_AWS.sh > upload_files_to_AWS_out.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate
conda activate ft-pipeline

#S3URL="s3://stergachis-public1/yhhc/Isoseq_results_for_shiny/"
S3URL="s3://stergachis-public1/yhhc/Isoseq_results_for_shiny/"

##########################################################################
# Upload each file to AWS
##########################################################################

aws s3 cp "/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/10.12.24/5.Test_statistics/hank_test_statistic_results_for_shiny.rds" $S3URL
