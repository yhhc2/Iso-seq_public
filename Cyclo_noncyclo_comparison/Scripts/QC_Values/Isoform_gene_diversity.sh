#!/bin/bash

# I have a file called data_combined_full.csv with these columns Isoform_PBid,Sample,cyclo_count,noncyclo_count,associated_gene
# For each Sample in data_combined_full.csv, look individually at the cyclo_count column. Only look at rows with >0 cyclo_count. Get all the unique values
# for isoform_pbid. This is the amount of cyclo unique isoforms. Do the same thing but with noncyclo_count to get noncyclo unique isoforms. 
# I want an output file with these columns: Sample, cyclo unique isoforms, cyclo total reads, noncyclo unique isoforms, noncyclo total reads

# Usage: bash Isoform_diversity.sh data_combined_full.csv unique_isoforms_counts.csv 

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <INPUT_FILE> <OUTPUT_FILE>"
    exit 1
fi

# Input and output file variables
INPUT_FILE=$1
OUTPUT_FILE=$2

# Initialize the output file with headers
echo -e "Sample\tCyclo Unique Isoforms\tCyclo Total Reads\tCyclo Unique Isoforms >1\tCyclo Unique Isoforms >10\tCyclo Unique Isoforms >100\tCyclo Unique Isoforms >1000\tNoncyclo Unique Isoforms\tNoncyclo Total Reads\tNoncyclo Unique Isoforms >1\tNoncyclo Unique Isoforms >10\tNoncyclo Unique Isoforms >100\tNoncyclo Unique Isoforms >1000" > $OUTPUT_FILE

# Extract unique sample names from the input file
samples=$(awk -F',' 'NR > 1 {print $2}' "$INPUT_FILE" | sort | uniq)

# Process each sample
for sample in $samples; do

    echo "$sample"

    # Cyclo analysis: Filter rows with cyclo_count > 0 for the current sample
    cyclo_data=$(awk -F',' -v sample="$sample" 'NR > 1 && $2 == sample && $3 > 0 {print $1 "\t" $3}' "$INPUT_FILE")
    
    # Count unique isoforms for cyclo
    cyclo_unique_isoforms=$(echo "$cyclo_data" | awk '{print $1}' | sort | uniq | wc -l)
    
    # Sum total reads for cyclo
    cyclo_total_reads=$(echo "$cyclo_data" | awk '{sum += $2} END {print sum}')
    
    # Count unique isoforms with >1 reads for cyclo
    cyclo_isoforms_gt1=$(echo "$cyclo_data" | awk '$2 > 1 {print $1}' | sort | uniq | wc -l)
    
    # Count unique isoforms with >10 reads for cyclo
    cyclo_isoforms_gt10=$(echo "$cyclo_data" | awk '$2 > 10 {print $1}' | sort | uniq | wc -l)
    
    # Count unique isoforms with >100 reads for cyclo
    cyclo_isoforms_gt100=$(echo "$cyclo_data" | awk '$2 > 100 {print $1}' | sort | uniq | wc -l)
    
    # Count unique isoforms with >1000 reads for cyclo
    cyclo_isoforms_gt1000=$(echo "$cyclo_data" | awk '$2 > 1000 {print $1}' | sort | uniq | wc -l)

    # Noncyclo analysis: Filter rows with noncyclo_count > 0 for the current sample
    noncyclo_data=$(awk -F',' -v sample="$sample" 'NR > 1 && $2 == sample && $4 > 0 {print $1 "\t" $4}' "$INPUT_FILE")
    
    # Count unique isoforms for noncyclo
    noncyclo_unique_isoforms=$(echo "$noncyclo_data" | awk '{print $1}' | sort | uniq | wc -l)
    
    # Sum total reads for noncyclo
    noncyclo_total_reads=$(echo "$noncyclo_data" | awk '{sum += $2} END {print sum}')
    
    # Count unique isoforms with >1 reads for noncyclo
    noncyclo_isoforms_gt1=$(echo "$noncyclo_data" | awk '$2 > 1 {print $1}' | sort | uniq | wc -l)
    
    # Count unique isoforms with >10 reads for noncyclo
    noncyclo_isoforms_gt10=$(echo "$noncyclo_data" | awk '$2 > 10 {print $1}' | sort | uniq | wc -l)
    
    # Count unique isoforms with >100 reads for noncyclo
    noncyclo_isoforms_gt100=$(echo "$noncyclo_data" | awk '$2 > 100 {print $1}' | sort | uniq | wc -l)
    
    # Count unique isoforms with >1000 reads for noncyclo
    noncyclo_isoforms_gt1000=$(echo "$noncyclo_data" | awk '$2 > 1000 {print $1}' | sort | uniq | wc -l)

    # Append results to the output file for the current sample
    echo -e "$sample\t$cyclo_unique_isoforms\t$cyclo_total_reads\t$cyclo_isoforms_gt1\t$cyclo_isoforms_gt10\t$cyclo_isoforms_gt100\t$cyclo_isoforms_gt1000\t$noncyclo_unique_isoforms\t$noncyclo_total_reads\t$noncyclo_isoforms_gt1\t$noncyclo_isoforms_gt10\t$noncyclo_isoforms_gt100\t$noncyclo_isoforms_gt1000" >> $OUTPUT_FILE
done

echo "Analysis complete. Results saved to '$OUTPUT_FILE'."
