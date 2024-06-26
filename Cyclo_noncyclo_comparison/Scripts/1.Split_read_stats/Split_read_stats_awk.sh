#!/bin/bash
# Hank Cheng
# 2024/03

# Usage: bash 1.Split_read_stats_awk.sh /path/to/read_stat.txt /path/to/sample_file.tsv > output.txt 2>&1

# Check if proper arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <read_stat_file> <sample_file>"
    exit 1
fi

# Assign input files from command line arguments
read_stat_file="$1"
sample_file="$2"


# Pre-process the sample_file to create a map of PS# to id_cyclo_noncyclo
awk -F '\t' 'NR > 1 {print $1, $10 "_" $11}' "$sample_file" > /tmp/ps_id_cyclo_map.txt

# Process the read_stat.txt file with awk, using the pre-created map
awk -v mapfile=/tmp/ps_id_cyclo_map.txt '
BEGIN {
    # Load the mapping file into the info array
    FS=" ";
    while ((getline < mapfile) > 0) {
        info[$1] = $2;
    }
    close(mapfile);

    # Set FS to tab for processing read_stat.txt
    FS="\t";
}

{
    # Extract PS# from the read name
    split($1, arr, "_");
    ps = arr[1];
    
    filename = ps "_" info[ps] "_read_stat.txt";

    # Output the current line to the corresponding file
    print $0 > filename;
    
}' "$read_stat_file"


# Count the total number of reads in the original read_stat.txt
total_reads=$(wc -l < "$read_stat_file")

# Count the total number of reads across all smaller files
total_reads_in_smaller_files=$(cat *_read_stat.txt | wc -l)

# Verification step
if [ "$total_reads" -eq "$total_reads_in_smaller_files" ]; then
    echo "Verification successful: The number of reads matches."
else
    echo "Verification failed: The number of reads does not match."
    echo "Reads in original file: $total_reads"
    echo "Reads in smaller files: $total_reads_in_smaller_files"
fi

# Calculate and report the number of lines in each output file
echo "Line counts for each output file:"
for file in *_read_stat.txt; do
    count=$(wc -l < "$file")
    echo "$file: $count lines"
done