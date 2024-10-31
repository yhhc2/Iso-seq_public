#!/bin/bash

# Usage: ./vep_annotate_and_filter_conda_gnomad_parallel.sh <vcf_list_file> <output_dir> > vep_annotate_and_filter_conda_gnomad_parallel_out.txt 2>&1

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <vcf_list_file> <output_dir>"
  exit 1
fi

VCF_LIST=$1
OUTPUT_DIR=$2

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate vep_environment

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Function to process each VCF file
process_vcf() {
  INPUT_VCF=$1
  OUTPUT_DIR=$2
  
  BASENAME=$(basename $INPUT_VCF .vcf.gz)
  FILTERED_VCF="${OUTPUT_DIR}/${BASENAME}_high_qual.vcf.gz"
  ANNOTATED_TXT="${OUTPUT_DIR}/${BASENAME}_annotated.txt"
  FILTERED_TXT="${OUTPUT_DIR}/${BASENAME}_filtered.txt"
  GENE_LIST="${OUTPUT_DIR}/${BASENAME}_gene_list.txt"

  echo "Processing $INPUT_VCF..."

  # Step 1: Filter for high-quality variants (QUAL > 30)
  echo "Filtering high-quality variants from $INPUT_VCF..."
  bcftools filter -i 'QUAL > 30' $INPUT_VCF -Oz -o $FILTERED_VCF
  tabix -p vcf $FILTERED_VCF

  if [ ! -f "$FILTERED_VCF" ]; then
    echo "Filtering failed for $INPUT_VCF, skipping."
    return
  fi

  # Step 2: Annotate the high-quality filtered VCF using VEP
  echo "Annotating $FILTERED_VCF with VEP..."
  vep --input_file $FILTERED_VCF --output_file $ANNOTATED_TXT \
    --offline --cache --dir_cache /mmfs1/gscratch/stergachislab/yhhc/tools/vep_cache \
    --force_overwrite \
    --everything --fork 4 --assembly GRCh38 \
    --plugin CADD,snv=/mmfs1/gscratch/stergachislab/yhhc/tools/CADD/whole_genome_SNVs.tsv.gz,indels=/mmfs1/gscratch/stergachislab/yhhc/tools/CADD/gnomad.genomes.r4.0.indel.tsv.gz \
    --quiet

  if [ ! -f "$ANNOTATED_TXT" ]; then
    echo "Annotation failed for $FILTERED_VCF, skipping."
    return
  fi

  # Step 3: Filter based on MAX_AF and CADD_PHRED values in the Extra column
  echo "Filtering variants with MAX_AF < 0.01 and CADD_PHRED > 10 from $ANNOTATED_TXT..."
  awk -F "\t" 'BEGIN {OFS="\t"} 
  NR==1 { print; next } 
  /^##/ { next }  # Skip lines starting with ##
  {
    # Extract MAX_AF and CADD_PHRED from the Extra field (assuming it is the last column)
    extra=$NF;
    
    # Match and extract MAX_AF (handling scientific notation)
    match(extra, /MAX_AF=([0-9.]+[eE]?[-+]?[0-9]*)/, arr_maxaf);
    
    # Match and extract CADD_PHRED
    match(extra, /CADD_PHRED=([0-9.]+)/, arr_cadd);

    # Apply both filters: MAX_AF < 0.01 and CADD_PHRED > 15
    if (arr_maxaf[1] < 0.01 && arr_cadd[1] > 15) print;
  }' "$ANNOTATED_TXT" > "$FILTERED_TXT"

  # Step 4: Get list of HGNC gene names with rare and high-quality variants.
  echo "Extracting gene names from $FILTERED_TXT..."
  
  # awk -F "\t" '{
  #   # Skip header lines
  #   if ($0 ~ /^##/ || NR==1) next;

  #   # Extract the Extra field (assuming it is the last field)
  #   extra = $NF;

  #   # Find and print the gene name (SYMBOL)
  #   if (match(extra, /SYMBOL=([^;]+)/, arr)) {
  #     print arr[1];  # arr[1] contains the gene name after SYMBOL=
  #   }
  # }' "$FILTERED_TXT" | sort | uniq > "$GENE_LIST"

  awk -F "\t" '{
    # Skip header lines
    if ($0 ~ /^##/ || NR==1) next;

    # Extract variant position (1st column) and variant effect (7th column)
    variant_position = $1;
    variant_effect = $7;

    # Extract the Extra field (assuming it is the last field)
    extra = $NF;

    # Find the gene name (SYMBOL)
    if (match(extra, /SYMBOL=([^;]+)/, arr_symbol)) {
      gene_name = arr_symbol[1];  # arr_symbol[1] contains the gene name after SYMBOL=

      # Match and extract CADD_PHRED
      if (match(extra, /CADD_PHRED=([0-9.]+)/, arr_cadd)) {
        cadd_phred = arr_cadd[1];
      } else {
        cadd_phred = "N/A";  # Default if CADD_PHRED not found
      }

      # Match and extract MAX_AF (handle scientific notation)
      if (match(extra, /MAX_AF=([0-9.]+[eE]?[-+]?[0-9]*)/, arr_maxaf)) {
        max_af = arr_maxaf[1];
      } else {
        max_af = "N/A";  # Default if MAX_AF not found
      }

      # Store variant information in an array for each gene
      gene_variants[gene_name] = gene_variants[gene_name] "[" variant_position ", " variant_effect ", " cadd_phred ", " max_af "], ";
    }
  }
  END {
    # Print each gene followed by its variants
    for (gene in gene_variants) {
      # Remove the trailing comma and space from the variants list
      variants = gene_variants[gene];
      sub(/, $/, "", variants);  # Remove trailing ", " at the end of the list

      print gene "\t" variants;
    }
  }' "$FILTERED_TXT" | sort | uniq > "$GENE_LIST"

  echo "Gene names have been extracted and saved to $GENE_LIST."
}

export -f process_vcf  # Export the function so GNU Parallel can use it

# Detect the number of CPUs
NUM_CPUS=$(nproc)  # Detect number of CPUs automatically (use `getconf _NPROCESSORS_ONLN` as an alternative)
echo "Detected $NUM_CPUS CPUs for parallel processing."

# Use GNU Parallel to process VCFs in parallel
cat $VCF_LIST | parallel -j $NUM_CPUS process_vcf {} $OUTPUT_DIR

echo "All VCF files processed. Results are in $OUTPUT_DIR."

# Deactivate the conda environment
#conda deactivate
