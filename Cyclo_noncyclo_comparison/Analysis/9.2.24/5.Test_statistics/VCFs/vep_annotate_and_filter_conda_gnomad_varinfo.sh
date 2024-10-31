#!/bin/bash

# Shell script to create a conda environment, install necessary tools, download gnomAD plugin and data,
# annotate VCF files using VEP with gnomAD, and filter for high-quality and rare variants.
# Outputs non-compressed VCF files.

# Usage: ./vep_annotate_and_filter_conda_gnomad.sh <vcf_list_file> <output_dir> > vep_annotate_and_filter_conda_gnomad_out.txt 2>&1
# sh vep_annotate_and_filter_conda_gnomad_varinfo.sh vcf_list_file.txt output_dir > vep_annotate_and_filter_conda_gnomad_varinfo_out.txt 2>&1

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <vcf_list_file> <output_dir>"
  exit 1
fi

# Assign input arguments to variables
VCF_LIST=$1
OUTPUT_DIR=$2

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate vep_environment

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Process each VCF in the list
while read INPUT_VCF; do
  if [ ! -f "$INPUT_VCF" ]; then
    echo "VCF file $INPUT_VCF not found, skipping."
    continue
  fi

  BASENAME=$(basename $INPUT_VCF .vcf.gz)
  FILTERED_VCF="${OUTPUT_DIR}/${BASENAME}_high_qual.vcf.gz"
  ANNOTATED_TXT="${OUTPUT_DIR}/${BASENAME}_annotated.txt"
  FILTERED_TXT="${OUTPUT_DIR}/${BASENAME}_filtered.txt"
  GENE_LIST="${OUTPUT_DIR}/${BASENAME}_gene_list.txt"  # Add missing GENE_LIST variable

  echo "Processing $INPUT_VCF..."

  # # Step 1: Filter for high-quality variants (QUAL > 30)
  # echo "Filtering high-quality variants from $INPUT_VCF..."
  # bcftools filter -i 'QUAL > 30' $INPUT_VCF -Oz -o $FILTERED_VCF
  # tabix -p vcf $FILTERED_VCF

  # if [ ! -f "$FILTERED_VCF" ]; then
  #   echo "Filtering failed for $INPUT_VCF, skipping."
  #   continue
  # fi

  # # Step 2: Annotate the high-quality filtered VCF using VEP with the default (non-VCF) output
  # echo "Annotating $FILTERED_VCF with VEP (default output format)..."
  # vep --input_file $FILTERED_VCF --output_file $ANNOTATED_TXT \
  #   --offline --cache --dir_cache /mmfs1/gscratch/stergachislab/yhhc/tools/vep_cache \
  #   --force_overwrite \
  #   --everything --fork 4 --assembly GRCh38 \
  #   --plugin CADD,snv=/mmfs1/gscratch/stergachislab/yhhc/tools/CADD/whole_genome_SNVs.tsv.gz,indels=/mmfs1/gscratch/stergachislab/yhhc/tools/CADD/gnomad.genomes.r4.0.indel.tsv.gz \
  #   --quiet

  # if [ ! -f "$ANNOTATED_TXT" ]; then
  #   echo "Annotation failed for $FILTERED_VCF, skipping."
  #   continue
  # fi

  # # Step 3: Filter based on MAX_AF and CADD_PHRED values in the Extra column
  # echo "Filtering variants with MAX_AF < 0.01 and CADD_PHRED > 10 from $ANNOTATED_TXT..."

  # awk -F "\t" 'BEGIN {OFS="\t"} 
  # NR==1 { print; next } 
  # /^##/ { next }  # Skip lines starting with ##
  # {
  #   # Extract MAX_AF and CADD_PHRED from the Extra field (assuming it is the last column)
  #   extra=$NF;
    
  #   # Match and extract MAX_AF (handling scientific notation)
  #   match(extra, /MAX_AF=([0-9.]+[eE]?[-+]?[0-9]*)/, arr_maxaf);
  #   if (!arr_maxaf[1]) arr_maxaf[1] = "N/A";  # Add safety check for MAX_AF

  #   # Match and extract CADD_PHRED
  #   match(extra, /CADD_PHRED=([0-9.]+)/, arr_cadd);
  #   if (!arr_cadd[1]) arr_cadd[1] = "N/A";  # Add safety check for CADD_PHRED

  #   # Apply both filters: MAX_AF < 0.01 and CADD_PHRED > 15
  #   if (arr_maxaf[1] < 0.01 && arr_cadd[1] > 15) print;
  # }' "$ANNOTATED_TXT" > "$FILTERED_TXT"

  # Step 4: Extract gene names with variant info
  echo "Extracting gene names and variant information from $FILTERED_TXT..."

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
      gene_name = arr_symbol[1];

      # Match and extract CADD_PHRED
      if (match(extra, /CADD_PHRED=([0-9.]+)/, arr_cadd)) {
        cadd_phred = arr_cadd[1];
      } else {
        cadd_phred = "N/A";
      }

      # Match and extract MAX_AF (handle scientific notation)
      if (match(extra, /MAX_AF=([0-9.]+[eE]?[-+]?[0-9]*)/, arr_maxaf)) {
        max_af = arr_maxaf[1];
      } else {
        max_af = "N/A";
      }

      # Create a unique key combining all variant information
      variant_info = "[" variant_position ", " variant_effect ", " cadd_phred ", " max_af "]";

      # Check if this variant for the gene is already stored
      if (!(variant_info in seen[gene_name])) {
        gene_variants[gene_name] = gene_variants[gene_name] variant_info ", ";
        seen[gene_name][variant_info] = 1;
      }
    }
  }
  END {
    # Print each gene followed by its variants
    for (gene in gene_variants) {
      variants = gene_variants[gene];
      sub(/, $/, "", variants);  # Remove trailing ", " at the end of the list
      print gene "\t" variants;
    }
  }' "$FILTERED_TXT" | sort | uniq > "$GENE_LIST"

  echo "Gene names and variant information have been extracted and saved to $GENE_LIST."

done < $VCF_LIST

echo "All VCF files processed. Results are in $OUTPUT_DIR."

# Deactivate the conda environment
conda deactivate
