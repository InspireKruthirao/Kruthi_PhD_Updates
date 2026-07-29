#!/bin/bash
PROTEIN_LIST="/work/microbiome/users/kruthi/intermediate_results/egg_nog/protein_list.txt"
DATA_DIR_SHARED="/work/microbiome/users/kruthi/intermediate_results/egg_nog/databases"
OUTPUT_DIR="/work/microbiome/users/kruthi/intermediate_results/egg_nog/results"

PROTEIN=$(sed -n "${PBS_ARRAY_INDEX}p" "$PROTEIN_LIST")
if [ -z "$PROTEIN" ]; then
    echo "No protein file found for array index ${PBS_ARRAY_INDEX}"
    exit 1
fi

SAMPLE=$(basename "$PROTEIN" "_proteins.faa")
SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "$SAMPLE_OUT"

# copy DB to local scratch (TMPDIR)
export EGGNOG_DATA_DIR="${TMPDIR}/eggnog_data"
mkdir -p "$EGGNOG_DATA_DIR"
cp "${DATA_DIR_SHARED}"/eggnog.db "$EGGNOG_DATA_DIR/"
cp "${DATA_DIR_SHARED}"/eggnog.taxa.db* "$EGGNOG_DATA_DIR/"
cp "${DATA_DIR_SHARED}"/eggnog_proteins.dmnd "$EGGNOG_DATA_DIR/"
# ^ adjust these cp lines to match exactly what files exist in databases/ — check with:
#   ls -lh /work/microbiome/users/kruthi/intermediate_results/egg_nog/databases

echo "Processing ${SAMPLE}"
emapper.py \
    -m diamond \
    --itype proteins \
    --cpu 12 \
    --data_dir "$EGGNOG_DATA_DIR" \
    -i "$PROTEIN" \
    --output "$SAMPLE" \
    --output_dir "$SAMPLE_OUT" \
    --resume
echo "Finished ${SAMPLE}"
