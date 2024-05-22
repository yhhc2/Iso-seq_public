# !/bin/bash
# Hank Cheng
# 7/1/2023
# Usage: sh Cyclo_noncyclo_comparison_general_collapse_isoform.sh > Cyclo_noncyclo_comparison_general_collapse_isoform_out.txt 2>&1


##########################################################################
# Define inputs
##########################################################################

#read_stat_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/11.22.23_merge_aligned_bams/2.Split_read_stats/UDN052264_read_stat.txt"
#classification_file="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/Merge_more_than_two_bams/11.22.23_merge_aligned_bams/1.Isoseq/08_pigeon/sorted_classification.filtered_lite_classification.txt"

#file_with_noncyclo_readIDs="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/ReadIDs/UDN052264_noncyclo_readIDs.txt"
#file_with_cyclo_readIDs="/mmfs1/gscratch/stergachislab/yhhc/projects/Iso-seq/Cyclo_noncyclo_comparison/ReadIDs/UDN052264_cyclo_readIDs.txt"

read_stat_file=$1
classification_file=$2

file_with_noncyclo_readIDs=$3
file_with_cyclo_readIDs=$4


# Select a single gene for positive control
selectedGene="MFN2"

##########################################################################
# Organize inputs
##########################################################################

# For some reason, the UDN633333 isoseq files had an additional second column that contains the read length of reads.
# So when I wrote my script, I relied on that. However, it looks like the isoseq files that I generated with my
# Isoseq_processing.sh script doesn't have this additional second column. So to make this script work, I'll
# just add a dummy second column. 

awk '{print $1, "000", $2}' $read_stat_file > read_stat_file_added_secondCol.txt

read_stat_file="read_stat_file_added_secondCol.txt"


##########################################################################
# 1. Add gene name to each row of this file collapse_isoforms.read_stat.txt, 
# using the information in this file isoseq_classification.txt
##########################################################################

# Although these files are txt files, it appears that the delimiter is a space instead of tab. 
# So all awk commands I used use the default space as the delimiter. 

# CHECK if read_stat_file is space delimited by seeing if the first row is space delimited. 
awk 'NR==1 { if (NF<=1) print "Error: First row of read_stat_file is not space delimited" }' $read_stat_file

# CHECK if classification_file is space delimited by seeing if the first row is space delimited.
awk 'NR==1 { if (NF<=1) print "Error: First row of classification_file is not space delimited" }' $classification_file


# CHECK if column 1 of read_stat_file contains smrt/zmw/ccs read ID. The value in the column should contain "ccs"
awk 'NR>1 {if ($1 !~ /ccs/) {error=1; print "Error: Line " NR ": column 1 of read_stat_file does not contain the smrt/zmw/ccs read ID"}} END{if(error==1) exit 1}' $read_stat_file

# CHECK if column 3 of read_stat_file contains the isoform ID (has PB identifier). The value in the column should contain "PB"
awk 'NR>1 {if ($3 !~ /PB/) {error=1; print "Error: Line " NR ": column 3 of read_stat_file does not contain the isoform ID"}} END{if(error==1) exit 1}' $read_stat_file

# CHECK if column 1 of classification_file contains the isoform ID (has PB identifier). The value in the column should contain "PB"
awk 'NR>1 {if ($1 !~ /PB/) {error=1; print "Error: Line " NR ": column 1 of classification_file does not contain the isoform ID"}} END{if(error==1) exit 1}' $classification_file

# CHECK if column 7 of classification_file contains the gene name. The column title should be "associated_gene"
awk 'NR==1 && $7 == "associated_gene" {found=1; exit} END {if(found!=1) print "Error: column 7 of classification_file does not contain the associated gene"}' $classification_file


# Need to first sort based on isoform ID.

sort -k3 $read_stat_file > read_stat_sorted.txt

sort -k1 $classification_file > classification_sorted.txt

# Then use join command to add the gene name to the collapse_isoforms.read_stat.txt file. Also add if the read corresponds to novel or canonical isoform. 
# Also add the subcategory label. 

#join -1 3 -2 1 -a 1 -o 1.1 1.2 1.3 2.7 read_stat_sorted.txt classification_sorted.txt > read_stat_sorted_with_gene.txt

# If the gene name is not found, then put N/A as the gene name. Do the same if the novel/non-novel label is not found. This will add N/A as gene name to reads that are not found in the classification file. 
join -1 3 -2 1 -a 1 -o 1.1 1.2 1.3 2.7 2.8 2.15 read_stat_sorted.txt classification_sorted.txt | awk '{ if ($4 == "") $4 = "N/A"; print }' | awk '{ if ($5 == "") $5 = "N/A"; print }' | awk '{ if ($6 == "") $6 = "N/A"; print }' > read_stat_sorted_with_gene.txt


# CHECK if read_stat_sorted_with_gene.txt has 6 columns.
awk 'BEGIN {found=0} NF==6 {found=1; exit} END {if(found!=1) print "Error: read_stat_sorted_with_gene.txt does not have the expected number of columns"}' read_stat_sorted_with_gene.txt
# 1st column smrt/zmw/ccs ID.
awk 'BEGIN {found=0} $1 ~ /ccs/ {found=1; exit} END {if (found==0) print "Error: The first column of read_stat_sorted_with_gene does not appear to contain smrt/zmw/ccs IDs"}' read_stat_sorted_with_gene.txt
# 2nd column is read length.
# 3rd column is isoform PB ID.
awk 'BEGIN {found=0} $3 ~ /PB/ {found=1; exit} END {if (found==0) print "Error: The third column of read_stat_sorted_with_gene does not appear to contain the isoform PB ID"}' read_stat_sorted_with_gene.txt
# 4th column is gene name. Assume that the gene name contains alphabet letters. 
awk 'BEGIN {found=0} $4 ~ /[^0-9]+/ {found=1; exit} END {if (found==0) print "Error: The fourth column of read_stat_sorted_with_gene does not appear to contain gene names"}' read_stat_sorted_with_gene.txt
# 5th column is whether the read corresponds to novel isoform or not.
# 6th column is what subcategory the isoform belogs in. 
# 7th column will be added below and will indicate if it came from cyclo or noncyclo. 

##########################################################################
# 2.For each row in collapse_isoforms.read_stat.txt, determine if the read is cyclo or non-cyclo. 
##########################################################################

# # This file contains read IDs for non-cycloheximide sample in the first column. 
awk 'NR==FNR{a[$1];next} $1 in a{$7=$7"noncyclo"}1' $file_with_noncyclo_readIDs read_stat_sorted_with_gene.txt > read_stat_sorted_with_gene_noncycloLabeled.txt

# # This file contains read IDs for cycloheximide sample in the first column. 
awk 'NR==FNR{a[$1];next} $1 in a{$7=$7"cyclo"}1' $file_with_cyclo_readIDs read_stat_sorted_with_gene_noncycloLabeled.txt > read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled.txt

# CHECK if all reads have a cyclo or noncyclo label at the 7th column.
awk '{if ($7 != "cyclo" && $7 != "noncyclo") {print "Error: In read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled.txt, row " NR " does not contain whether the read is from cyclo or noncyclo in the 5th column" }}' read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled.txt


##########################################################################
# 3.Get all reads that correspond to a particular gene
##########################################################################

if [[ "$selectedGene" == "GENOME_WIDE" ]]; then
    cp read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled.txt read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled_selectedGene.txt
else
    # Get all reads that correspond to a particular gene. Like MFN2
    awk -v gene="$selectedGene" '$4 == gene' read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled.txt > read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled_selectedGene.txt

fi

##########################################################################
# 4. See if there are any specific gene isoforms for the selected gene that is overwhelmingly 
# enriched in the cyclo sample
##########################################################################


if [[ "$selectedGene" == "GENOME_WIDE" ]]; then
    cp $classification_file classification_file_selectedGene.txt
else
    awk -v gene="$selectedGene" 'NR==1 || $7 == gene' $classification_file > classification_file_selectedGene.txt

fi

# Replace all ENST values with "non-novel" 
awk '{ if ($5 ~ /ENST/) $5 = "non-novel"; print }' read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled.txt > read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled_nonnovelLabeled.txt

# I have a file where I want to collapse rows that have the same 4th (gene), 5th (novel or not), and 7th (cyclo/noncyclo) column values. And I want to add a column that indicates the number of rows collapsed.
awk 'BEGIN { FS=OFS=" " } { key = $4 OFS $5 OFS $7; count[key]++; } END { for (key in count) { print key, count[key]; } }' read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled_nonnovelLabeled.txt > collapsed_by_gene_file.txt

# I have a file where I want to collapse rows that have the same 3rd (isoform) and 7th (cyclo/noncyclo) column values. And I want to add a column that indicates the number of rows collapsed.
awk 'BEGIN { FS=OFS=" " } { key = $3 OFS $7; count[key]++; } END { for (key in count) { print key, count[key]; } }' read_stat_sorted_with_gene_noncycloLabeled_cycloLabeled_nonnovelLabeled.txt > collapsed_by_isoform_file.txt


#For each unique PB identifier, I want to create two columns. One column that shows the numerical value for the cyclo category and another column that shows the numerical value for the noncyclo category

awk '
{
    id = $1;
    category = $2;
    value = $3;
    if (category == "cyclo") {
        cyclo[id] = value;
    } else if (category == "noncyclo") {
        noncyclo[id] = value;
    }
}
END {
    for (id in cyclo) {
        print id, cyclo[id], (noncyclo[id] ? noncyclo[id] : 0);
    }
    for (id in noncyclo) {
        if (!cyclo[id]) {
            print id, 0, noncyclo[id];
        }
    }
}' collapsed_by_isoform_file.txt > collapsed_by_isoform_file_cyclo_noncyclo_counts.txt



##########################################################################
# 5.Add classification info to the isoform collapsed file
##########################################################################

fileA="collapsed_by_isoform_file_cyclo_noncyclo_counts.txt"
fileB=$classification_file
outputFile="collapsed_by_isoform_file_cyclo_noncyclo_counts_classified.txt"

# Write header to outputFile
echo "Isoform_PBid cyclo_count noncyclo_count structural_category associated_gene associated_transcript subcategory" > "$outputFile"

# Use awk to process files
awk -F ' ' '
    # Load fileB into an array
    NR == FNR {
        data[$1] = $6 " " $7 " " $8 " " $15;
        next;
    }
    # Process fileA and output to outputFile. FileA should have no header, so process every line.
    FNR > 0 {
        if ($1 in data) {
            print $1, $2, $3, data[$1];
        } else {
            print $1, $2, $3;
        }
    }
' "$fileB" "$fileA" >> "$outputFile"

# NOTE: Any PBid where classification information is not written into the output are PBids where the isoform is filtered
# out of the lite classification file. 


# CHECK if the number of lines in the output is what we expect.
expected_lines=$(wc -l < "$fileA")
output_lines=$(($(wc -l < "$outputFile") - 1)) #Subtract 1 because header is added to output file.
# Output file should have header + same number of lines as input file
if [ "$((output_lines))" -ne "$expected_lines" ]; then
    echo "Error: Number of lines in $outputFile does not match $fileA. Expected $expected_lines. Actual $output_lines."
    exit 1
fi

# CHECK column names in the final output
actual_header=$(head -1 "$outputFile")
expected_header="Isoform_PBid cyclo_count noncyclo_count structural_category associated_gene associated_transcript subcategory"

if [ "$actual_header" != "$expected_header" ]; then
    echo "Error: The header in the output file does not match the expected format."
    echo "Expected: $expected_header"
    echo "Actual: $actual_header"
    exit 1
fi

##########################################################################
# Clean up directory
##########################################################################

# Remove intermediate files. 
# Comment out this section to troubleshoot.

# Define the directory and the file to keep
DIRECTORY="."
KEEP_FILE="collapsed_by_isoform_file_cyclo_noncyclo_counts_classified.txt"

# Find all files in the directory, excluding the specified file
for file in "$DIRECTORY"/*; do
    if [ "$(basename "$file")" != "$KEEP_FILE" ]; then
        rm -f "$file"
    fi
done


##########################################################################
# Further analysis
##########################################################################

# collapsed_by_isoform_file_cyclo_noncyclo_counts_classified.txt file can be inputted into R
# for further analysis. 
