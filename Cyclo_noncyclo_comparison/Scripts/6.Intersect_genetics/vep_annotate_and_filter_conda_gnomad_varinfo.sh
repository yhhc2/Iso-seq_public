#!/bin/bash

# Shell script to create a conda environment, install necessary tools, download gnomAD plugin and data,
# annotate VCF files using VEP with gnomAD, and filter for high-quality and rare variants.
# Outputs non-compressed VCF files.

# Usage: ./vep_annotate_and_filter_conda_gnomad.sh <vcf_tsv_file> <output_dir> > vep_annotate_and_filter_conda_gnomad_out.txt 2>&1

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <vcf_tsv_file> <output_dir>"
  exit 1
fi

# Assign input arguments to variables
VCF_TSV=$1
OUTPUT_DIR=$2

# Activate the conda environment
source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate vep_environment

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Find the column number for "deepvariant_vcf"
VCF_COLUMN_NUM=$(head -n 1 "$VCF_TSV" | tr -d '\r' | tr '\t' '\n' | awk '{ gsub(/[[:space:]]+$/, ""); if ($0 == "deepvariant_vcf") print NR; }')

# Check if the column was found
if [ -z "$VCF_COLUMN_NUM" ]; then
  echo "Error: deepvariant_vcf column not found in the TSV file."
  exit 1
fi

# Extract unique VCF paths from the deepvariant_vcf column using the correct column number, and trim them of newline or space characters
VCF_PATHS=$(cut -f"$VCF_COLUMN_NUM" "$VCF_TSV" | sort | uniq | sed 's/[[:space:]]*$//')

# Process each unique VCF file
for INPUT_VCF in $VCF_PATHS; do
  if [ ! -f "$INPUT_VCF" ]; then
    echo "VCF file $INPUT_VCF not found, skipping."
    continue
  fi

  BASENAME=$(basename $INPUT_VCF .vcf.gz)
  FILTERED_VCF="${OUTPUT_DIR}/${BASENAME}_high_qual.vcf.gz"
  ANNOTATED_TXT="${OUTPUT_DIR}/${BASENAME}_annotated.txt"
  FILTERED_TXT="${OUTPUT_DIR}/${BASENAME}_filtered.txt"
  GENE_LIST="${OUTPUT_DIR}/*${BASENAME}_gene_list.txt" #The wildcard is needed so that BCH_0252-01.PS00404.deepvariant_gene_list.txt and PS00404.deepvariant_gene_list.txt are recognized as the same VCF.

  ## Skip the VCF if the gene list file already exists
  #if [ -f "$GENE_LIST" ]; then
  #  echo "Gene list for $INPUT_VCF already exists, skipping."
  #  continue
  #fi

  # Skip the VCF if any gene list file matching the pattern already exists
  if find "${OUTPUT_DIR}" -type f -name "*${BASENAME}_gene_list.txt" | grep -q .; then
    echo "Gene list for $INPUT_VCF already exists, skipping."
    continue
  fi

  echo "Processing $INPUT_VCF..."

  # Step 1: Filter for high-quality variants (QUAL > 30)
  echo "Filtering high-quality variants from $INPUT_VCF..."
  bcftools filter -i 'QUAL > 30' $INPUT_VCF -Oz -o $FILTERED_VCF
  tabix -p vcf $FILTERED_VCF

  if [ ! -f "$FILTERED_VCF" ]; then
    echo "Filtering failed for $INPUT_VCF, skipping."
    continue
  fi

  # Step 2: Annotate the high-quality filtered VCF using VEP with the default (non-VCF) output
  echo "Annotating $FILTERED_VCF with VEP (default output format)..."
  vep --input_file $FILTERED_VCF --output_file $ANNOTATED_TXT \
    --offline --cache --dir_cache /mmfs1/gscratch/stergachislab/yhhc/tools/vep_cache \
    --force_overwrite \
    --everything --fork 4 --assembly GRCh38 \
    --plugin CADD,snv=/mmfs1/gscratch/stergachislab/yhhc/tools/CADD/whole_genome_SNVs.tsv.gz,indels=/mmfs1/gscratch/stergachislab/yhhc/tools/CADD/gnomad.genomes.r4.0.indel.tsv.gz \
    --quiet

  if [ ! -f "$ANNOTATED_TXT" ]; then
    echo "Annotation failed for $FILTERED_VCF, skipping."
    continue
  fi

  # Step 3: Filter based on MAX_AF and CADD_PHRED values in the Extra column
  echo "Filtering variants with MAX_AF < 0.01 and CADD_PHRED > 10 from $ANNOTATED_TXT..."

  awk -F "\t" 'BEGIN {OFS="\t"} 
  NR==1 { print; next } 
  /^##/ { next }  # Skip lines starting with ##
  {
    extra=$NF;
    
    # Match and extract MAX_AF (handling scientific notation)
    match(extra, /MAX_AF=([0-9.]+[eE]?[-+]?[0-9]*)/, arr_maxaf);
    if (!arr_maxaf[1]) arr_maxaf[1] = "N/A";  # Add safety check for MAX_AF

    # Match and extract CADD_PHRED
    match(extra, /CADD_PHRED=([0-9.]+)/, arr_cadd);
    if (!arr_cadd[1]) arr_cadd[1] = "N/A";  # Add safety check for CADD_PHRED

    # Apply both filters: MAX_AF < 0.01 and CADD_PHRED > 15
    if (arr_maxaf[1] < 0.01 && arr_cadd[1] > 15) print;
  }' "$ANNOTATED_TXT" > "$FILTERED_TXT"

  # Step 4: Extract gene names with variant info
  echo "Extracting gene names and variant information from $FILTERED_TXT..."

  awk -F "\t" '{
    if ($0 ~ /^##/ || NR==1) next;

    variant_position = $1;
    variant_effect = $7;
    extra = $NF;

    if (match(extra, /SYMBOL=([^;]+)/, arr_symbol)) {
      gene_name = arr_symbol[1];
      
      if (match(extra, /CADD_PHRED=([0-9.]+)/, arr_cadd)) {
        cadd_phred = arr_cadd[1];
      } else {
        cadd_phred = "N/A";
      }

      if (match(extra, /MAX_AF=([0-9.]+[eE]?[-+]?[0-9]*)/, arr_maxaf)) {
        max_af = arr_maxaf[1];
      } else {
        max_af = "N/A";
      }

      variant_info = "[" variant_position ", " variant_effect ", " cadd_phred ", " max_af "]";

      if (!(variant_info in seen[gene_name])) {
        gene_variants[gene_name] = gene_variants[gene_name] variant_info ", ";
        seen[gene_name][variant_info] = 1;
      }
    }
  }
  END {
    for (gene in gene_variants) {
      variants = gene_variants[gene];
      sub(/, $/, "", variants);
      print gene "\t" variants;
    }
  }' "$FILTERED_TXT" | sort | uniq > "$GENE_LIST"

  echo "Gene names and variant information have been extracted and saved to $GENE_LIST."

done

echo "All VCF files processed. Results are in $OUTPUT_DIR."

# Deactivate the conda environment
conda deactivate
