#!/bin/bash

export EGGNOG_DATA_DIR=/work/microbiome/users/kruthi/intermediate_results/egg_nog/databases

PROTEIN_LIST="/work/microbiome/users/kruthi/intermediate_results/egg_nog/protein_list.txt"
DATA_DIR="/work/microbiome/users/kruthi/intermediate_results/egg_nog/databases"
OUTPUT_DIR="/work/microbiome/users/kruthi/intermediate_results/egg_nog/results"

PROTEIN=$(sed -n "${PBS_ARRAY_INDEX}p" "$PROTEIN_LIST")

if [ -z "$PROTEIN" ]; then
    echo "No protein file found for array index ${PBS_ARRAY_INDEX}"
    exit 1
fi

SAMPLE=$(basename "$PROTEIN" "_proteins.faa")

mkdir -p "${OUTPUT_DIR}/${SAMPLE}"

echo "Processing ${SAMPLE}"

emapper.py \
    -m diamond \
    --itype proteins \
    --cpu 12 \
    --data_dir "$DATA_DIR" \
    -i "$PROTEIN" \
    --output "$SAMPLE" \
    --output_dir "${OUTPUT_DIR}/${SAMPLE}" \
    --override

echo "Finished ${SAMPLE}"
