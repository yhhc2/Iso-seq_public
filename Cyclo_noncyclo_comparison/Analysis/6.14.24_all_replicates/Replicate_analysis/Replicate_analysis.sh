#!/bin/bash
# Hank Cheng
# 2024/05

# Usage bash Replicate_analysis.sh > Replicate_analysis_analysis_output.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate

conda activate r_env_per_isoform


#########################################################
# Define inputs
#########################################################

# Scripts:
Replicate_analysis="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/5.Replicate_analysis/Replicate_analysis.R"
Replicate_analysis_NormalizedFraction="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Scripts/3.Compare_samples/5.Replicate_analysis/Replicate_analysis_NormalizedFraction.R"

# Directories:
Compare_samples_Isoform="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/3.Compare_samples/1.Isoform"
Compare_samples_Gene="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq_public/Cyclo_noncyclo_comparison/Analysis/6.14.24_all_replicates/3.Compare_samples/2.Gene"

#########################################################
# Replicate analysis
#########################################################

# IMPORTANT: Remember to play around with the coverage threshold. 

mkdir Gene_TPM_0_1000
cd Gene_TPM_0_1000

# Gene-level TPM
Rscript "$Replicate_analysis" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  6 \
  "UDN212054 Gene Level TPM" \
  "UDN212054,UDN212054b" \
  "VTA1"
  > output.txt \
  2>&1

# Gene-level TPM
Rscript "$Replicate_analysis" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  6 \
  "UDN687128 Gene Level TPM" \
  "UDN687128,UDN687128b" \
  "VTA1"
  > output.txt \
  2>&1

cd .. 
mkdir Isoform_TPM_0_1000
cd Isoform_TPM_0_1000

# Isoform-level TPM
Rscript "$Replicate_analysis" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  6 \
  "UDN212054 Isoform Level TPM" \
  "UDN212054,UDN212054b" \
  "N/A"
  > output.txt \
  2>&1

Rscript "$Replicate_analysis" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  6 \
  "UDN687128 Isoform Level TPM" \
  "UDN687128,UDN687128b" \
  "N/A"
  > output.txt \
  2>&1

Rscript "$Replicate_analysis" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  6 \
  "UDN212054 and UDN687128 Isoform Level TPM" \
  "UDN212054,UDN687128" \
  "N/A"
  > output.txt \
  2>&1


cd ..
mkdir Gene_NormalizedFractionDifference
cd Gene_NormalizedFractionDifference

# Gene-level NormalizedFractionDifference
Rscript "$Replicate_analysis_NormalizedFraction" \
  "${Compare_samples_Gene}/data_combined_full.csv" \
  101 \
  "UDN212054 Gene Level NormalizedFractionDifference" \
  "UDN212054,UDN212054b" \
  "N/A"
  > output.txt \
  2>&1


cd ..
mkdir Isoform_NormalizedFractionDifference
cd Isoform_NormalizedFractionDifference

# Isoform-level NormalizedFractionDifference
Rscript "$Replicate_analysis_NormalizedFraction" \
  "${Compare_samples_Isoform}/data_combined_full.csv" \
  70 \
  "UDN212054 Isoform Level NormalizedFractionDifference" \
  "UDN212054,UDN212054b" \
  "PB.60131.293"
  > output.txt \
  2>&1

# Isoform-level NormalizedFractionDifference
Rscript "$Replicate_analysis_NormalizedFraction" \
  "${Compare_samples_Isoform}/data_combined_full.csv" \
  70 \
  "UDN212054 and UDN687128 Isoform Level NormalizedFractionDifference" \
  "UDN212054,UDN687128" \
  "PB.60131.293"
  > output.txt \
  2>&1


