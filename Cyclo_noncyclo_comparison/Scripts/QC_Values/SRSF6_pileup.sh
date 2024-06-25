#!/bin/bash

# Script that uses pileup at cassette exon in SRSF6 as a quality check to see that
# cycloheximide is working an intended.

# Usage: bash SRSF6_pileup.sh <INPUT_FILE> <REFERENCE_GENOME> <CHROMOSOME> <POSITION> <OUTPUT_FILE> > SRSF6_pileup.output.txt 2>&1

# Ensure the script stops if any command fails
set -e

# Check if the correct number of arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <INPUT_FILE> <REFERENCE_GENOME> <CHROMOSOME> <POSITION> <OUTPUT_FILE>"
    exit 1
fi

# Variables for input file, reference genome, chromosome, position, and final output file
INPUT_FILE=$1
REFERENCE=$2
CHR=$3
POS=$4
FINAL_OUTPUT_FILE=$5

# Temporary file to store intermediate results
TEMP_OUTPUT_FILE="temp_combined_results.txt"

# Initialize the temporary output file and print the header
echo -e "Source\tChromosome\tPosition\tReference_Base\tOriginal_Read_Depth\tAdjusted_Read_Depth\tExonic_Proportion\tRead_Bases\tBase_Qualities" > $TEMP_OUTPUT_FILE

# Read unique BAM files and associated metadata from the input file and process each one
awk -F'\t' 'NR > 1 {print $1, $10, $11, $14}' $INPUT_FILE | sort -k4 | uniq -f3 | while read SAMPLE ID CYCLO_NONCYCLO BAM_FILE; do
    SOURCE="${SAMPLE}_${ID}_${CYCLO_NONCYCLO}"
    echo "Processing BAM file: $BAM_FILE from source: $SOURCE"
    
    if [ -f "$BAM_FILE" ]; then
        # Generate the pileup file for the specific location
        samtools mpileup -f $REFERENCE $BAM_FILE -r $CHR:$POS-$POS > temp.pileup

        # Compare the pileup information and append to the temporary output file
        awk -v source="$SOURCE" '
        {
            skip_count = gsub(/[<>]/, "", $5);
            adjusted_read_depth = $4 - skip_count;
            exonic_proportion = ($4 > 0) ? adjusted_read_depth / $4 : 0;
            printf "%s\t%s\t%d\t%s\t%d\t%d\t%.2f\t%s\t%s\n", source, $1, $2, $3, $4, adjusted_read_depth, exonic_proportion, $5, $6;
        }' temp.pileup >> $TEMP_OUTPUT_FILE

        # Remove the temporary pileup file
        rm temp.pileup
    else
        # Append a line with N/A for the missing BAM file
        echo -e "${SOURCE}\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A" >> $TEMP_OUTPUT_FILE
    fi
done

# Move the temporary output file to the final output file
mv $TEMP_OUTPUT_FILE $FINAL_OUTPUT_FILE

echo "Processing complete. Results saved to $FINAL_OUTPUT_FILE"