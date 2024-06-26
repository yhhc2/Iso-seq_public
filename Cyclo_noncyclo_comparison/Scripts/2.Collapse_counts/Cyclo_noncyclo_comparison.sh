# !/bin/bash
# Hank Cheng
# 7/1/2023
# Usage: bash 1.Cyclo_noncyclo_comparison.sh /path/to/sample_file.tsv /path/to/split_read_stats /path/to/output_base_dir /path/to/classification.txt /path/to/helper_script.sh > output.txt 2>&1


# Run Cyclo_noncyclo_comparison_general_collapse_isoform.sh for each patient. 

# For each unique id in the Sample_names_and_bam_locations.tsv file, do this:

# Look in the 2.Split_read_stats folder for two files:
# - A file that contains _${id}_noncyclo in its name. Assign the path of this file into a variable called noncyclo_readstat.
# - A file that contains _${id}_cyclo in its name. Assign the path of this file into a variable called cylo_readstat.

# Concat fileA and fileB together and store the path of this file into a variable called combined_readstat. 

# Then run this script:
#sh Cyclo_noncyclo_comparison_general_collapse_isoform.sh $combined_readstat $classification $noncyclo_readstat $cyclo_readstat

#"PS100_UDN123_cyclo_read_stat.txt"
#"PS101_UDN123_noncyclo_read_stat.txt"

#############################################################################
# Script
#############################################################################

# Check for correct number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <sample_file> <split_read_stats_dir> <output_base_dir> <classification> <helper_script>"
    exit 1
fi

# Assign command line arguments to variables
sample_file="$1"
split_read_stats_dir="$2"
output_base_dir="$3"
classification="$4"
helper_script="$5"

# Extract unique IDs from the TSV file (assuming 'id' is in the 10th column and skipping the header)
unique_ids=$(cut -f10 $sample_file | tail -n +2 | sort | uniq)

# Iterate through each unique id
for id in $unique_ids; do
    # Define the files based on id
    noncyclo_readstat=$(find $split_read_stats_dir -type f -name "*_${id}_noncyclo*")
    cyclo_readstat=$(find $split_read_stats_dir -type f -name "*_${id}_cyclo*")

    # Check if both files exist
    if [[ -f "$noncyclo_readstat" && -f "$cyclo_readstat" ]]; then
        # Concatenate the files and save to a new file in a new directory for the id
        mkdir -p "${output_base_dir}/${id}"
        combined_readstat="${output_base_dir}/${id}/${id}_combined_readstat.txt"
        cat "$noncyclo_readstat" "$cyclo_readstat" > "$combined_readstat"
        
        # Navigate to the id's directory
        pushd "${output_base_dir}/${id}"

        # Run the provided script with the paths as arguments and direct the output to the id's directory
        sh "$helper_script" "$combined_readstat" "$classification" "$noncyclo_readstat" "$cyclo_readstat"

        # Return to the previous directory
        popd

	# Print new line so that the pushd and popd directory output are not all squished
	echo ""

    else
        echo "Error: Missing file(s) for ID $id"
    fi
done