#!/bin/bash
set -u
set -o pipefail

PROTEIN_LIST="/work/microbiome/users/kruthi/intermediate_results/egg_nog/protein_list.txt"
DATA_DIR_SHARED="/work/microbiome/users/kruthi/intermediate_results/egg_nog/databases"
OUTPUT_DIR="/work/microbiome/users/kruthi/intermediate_results/egg_nog/6aug_results"

if [ -z "${TMPDIR:-}" ]; then
    echo "TMPDIR is not set"
    exit 1
fi

PROTEIN=$(sed -n "${PBS_ARRAY_INDEX}p" "$PROTEIN_LIST")
if [ -z "$PROTEIN" ]; then
    echo "No protein file found for array index ${PBS_ARRAY_INDEX}"
    exit 1
fi

SAMPLE=$(basename "$PROTEIN" "_proteins.faa")
SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "$SAMPLE_OUT"

LOGFILE="${SAMPLE_OUT}/${SAMPLE}.log"
exec > >(tee -a "$LOGFILE") 2>&1

log() {
    echo "[$(date '+%F %T')] $*"
}

log "Starting sample: ${SAMPLE}"
log "Protein file: ${PROTEIN}"
log "Output dir: ${SAMPLE_OUT}"
log "TMPDIR: ${TMPDIR}"

EGGNOG_DATA_DIR="${TMPDIR}/eggnog_data"
mkdir -p "$EGGNOG_DATA_DIR"

if [ ! -f "${EGGNOG_DATA_DIR}/eggnog.db" ]; then
    log "Copying eggNOG databases to local scratch"
    cp "${DATA_DIR_SHARED}/eggnog.db" "${EGGNOG_DATA_DIR}/"
    cp "${DATA_DIR_SHARED}"/eggnog.taxa.db* "${EGGNOG_DATA_DIR}/"
    cp "${DATA_DIR_SHARED}/eggnog_proteins.dmnd" "${EGGNOG_DATA_DIR}/"
else
    log "Local eggNOG database already present in TMPDIR"
fi

for f in eggnog.db eggnog.taxa.db eggnog.taxa.db.traverse.pkl eggnog_proteins.dmnd; do
    if [ ! -e "${EGGNOG_DATA_DIR}/${f}" ]; then
        log "ERROR: missing ${EGGNOG_DATA_DIR}/${f}"
        exit 1
    fi
done

log "Running emapper.py"
emapper.py \
    -m diamond \
    --itype proteins \
    --cpu 12 \
    --data_dir "$EGGNOG_DATA_DIR" \
    -i "$PROTEIN" \
    --output "$SAMPLE" \
    --output_dir "$SAMPLE_OUT" \
    --override

EMAPPER_STATUS=$?
log "emapper.py exit status: ${EMAPPER_STATUS}"

ANN_FILE="${SAMPLE_OUT}/${SAMPLE}.emapper.annotations"
HITS_FILE="${SAMPLE_OUT}/${SAMPLE}.emapper.hits"
SEED_FILE="${SAMPLE_OUT}/${SAMPLE}.emapper.seed_orthologs"

if [ -s "$HITS_FILE" ]; then
    log "OK: hits file present (${HITS_FILE})"
else
    log "ERROR: hits file missing or empty (${HITS_FILE})"
fi

if [ -s "$ANN_FILE" ]; then
    log "OK: annotations file present (${ANN_FILE})"
else
    log "ERROR: annotations file missing or empty (${ANN_FILE})"
fi

if [ -s "$SEED_FILE" ]; then
    log "OK: seed ortholog file present (${SEED_FILE})"
elif [ -e "$SEED_FILE" ]; then
    log "WARNING: seed ortholog file is empty (${SEED_FILE})"
else
    log "WARNING: seed ortholog file missing (${SEED_FILE})"
fi

if [ "$EMAPPER_STATUS" -ne 0 ]; then
    log "FAILED: emapper.py exited with status ${EMAPPER_STATUS}"
    exit "$EMAPPER_STATUS"
fi

if [ ! -s "$ANN_FILE" ] || [ ! -s "$HITS_FILE" ]; then
    log "FAILED: required output files are incomplete"
    exit 2
fi

log "Finished sample: ${SAMPLE}"

#RUN THE JOB WITH BELOW COMMAND
mqsub --mem 20 -t 6 --cpus 12 --bg --segregated-log-file --name eggnog_array --array 58 -- bash eggnog_array.sh
