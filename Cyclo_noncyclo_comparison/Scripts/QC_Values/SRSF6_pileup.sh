#!/bin/bash

# Usage bash SRSF6_pileup.sh > SRSF6_pileup.output.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate

conda activate samtools

# Ensure the script stops if any command fails
set -e

# Variables for the reference genome, BAM file, and positions to compare
REFERENCE="/gscratch/stergachislab/assemblies/simple-names/hg38.fa"
BAM_FILE="/mmfs1/gscratch/stergachislab/asedeno/data/Hank_MASseq/8-tag/PS00429.tagged.bam"

# Cassette exon position
CHR1="chr20"
POS1=43459200

FINAL_OUTPUT_FILE="adjusted_pileup_results.txt"

# Step 1: Generate the pileup file for specific locations
echo "Generating pileup for specific locations..."
samtools mpileup -f $REFERENCE $BAM_FILE -r $CHR1:$POS1-$POS1 > temp1.pileup

# Step 2: Compare the pileup information
echo "Comparing pileup information..."

# Using awk to adjust read depth by subtracting reference skips and print comparison in a readable format.
# Reference skips are indicated by < and > characters in the pileup results. These are intronic regions.
awk '
BEGIN {
    printf "%-10s %-10s %-15s %-10s %-20s %-20s %-50s %-20s\n", "Chromosome", "Position", "Reference_Base", "Original_Read_Depth", "Adjusted_Read_Depth", "Exonic_Proportion", "Read_Bases", "Base_Qualities";
}
{
    skip_count = gsub(/[<>]/, "", $5);
    adjusted_read_depth = $4 - skip_count;
    exonic_proportion = ($4 > 0) ? adjusted_read_depth / $4 : 0;
    printf "%-10s %-10d %-15s %-10d %-20d %-20.2f %-50s %-20s\n", $1, $2, $3, $4, adjusted_read_depth, exonic_proportion, $5, $6;
}' temp1.pileup > $FINAL_OUTPUT_FILE

# Step 3: Cleanup (Optional)
# Remove intermediate files
rm temp1.pileup

echo "Comparison complete. Output saved to $FINAL_OUTPUT_FILE"