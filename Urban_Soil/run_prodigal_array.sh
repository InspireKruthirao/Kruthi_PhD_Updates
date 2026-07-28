#!/bin/bash

# Load conda and activate environment
source ~/miniconda3/bin/activate
conda activate kruthi

set -e

# Input list and output directory
ASSEMBLY_LIST="/work/microbiome/users/kruthi/intermediate_results/urban_soil/assembly_list.txt"
OUTPUT_DIR="/work/microbiome/users/kruthi/intermediate_results/urban_soil"

# Get the assembly
ASSEMBLY=$(sed -n "${PBS_ARRAY_INDEX}p" "$ASSEMBLY_LIST")

# Stop if no assembly 
if [ -z "$ASSEMBLY" ]; then
	    echo "No assembly found for array index ${PBS_ARRAY_INDEX}"
	        exit 1
fi

# Extract sample name
SAMPLE=$(basename "$ASSEMBLY" "_medaka_polypolish.fasta.PolcaCorrected.fa.gz")

echo "Processing ${SAMPLE}"

# Create output folder 
mkdir -p "${OUTPUT_DIR}/${SAMPLE}"

# Create temporary fasta
TMP_FASTA="/tmp/${SAMPLE}.fna"

gunzip -c "$ASSEMBLY" > "$TMP_FASTA"

# Prodigal
prodigal \
	    -i "$TMP_FASTA" \
	        -o "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_coords.gbk" \
		    -a "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_proteins.faa" \
		        -d "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_orfs.fna" \
			    -p meta


rm "$TMP_FASTA"

echo "Finished ${SAMPLE}"
