#!/bin/bash -l
#PBS -N eggnog_array
#PBS -J 1-58
#PBS -l select=1:ncpus=12:mem=64gb
#PBS -l walltime=04:00:00
#PBS -q cpu_batch_exec
#PBS -j oe
#PBS -m abe

# STAGE 7a - Run eggNOG-mapper (v2.1.12, diamond mode) on every
# sample's Prodigal protein predictions.
#
# NOTE ON RECONSTRUCTION
# -----------------------
# Like 04_run_mmseqs_95nt_64bit.sh, the original per-sample driver
# (eggnog_array.sh / eggnog2_array.sh, later renamed Run_Eggnog_Mapper.sh)
# was authored with `vim` and iterated on directly on the HPC, so its
# exact original content is not in shell history. This reconstruction
# uses the exact `emapper.py` CLI flags that were tested and confirmed
# working (see history around 2026-07-27/28), driven as a PBS job
# array over the 58 Urban Soil samples, one eggNOG run per sample.
#
# Input : <BASEDIR>/<sample>/<sample>_proteins.faa   (Prodigal proteins)
# Output: egg_nog/results/<sample>/<sample>.emapper.annotations
#         (consumed downstream by 08_US95_emapper.py)
#
# eggNOG-mapper database: v2.1.12, copied locally to
#   /work/microbiome/users/kruthi/intermediate_results/egg_nog/databases
# (from /work/microbiome/db/eggnog-mapper/2.1.12)

set -euo pipefail

BASEDIR="/work/microbiome/users/kruthi/intermediate_results/urban_soil"
EGGNOG_DB="/work/microbiome/users/kruthi/intermediate_results/egg_nog/databases"
RESULTS="/work/microbiome/users/kruthi/intermediate_results/egg_nog/results"
SAMPLE_LIST="${BASEDIR}/sample_list.txt"   # one sample name per line, 58 total

mkdir -p "${RESULTS}"

SAMPLE=$(sed -n "${PBS_ARRAY_INDEX}p" "${SAMPLE_LIST}")

if [ -z "${SAMPLE}" ]; then
    echo "No sample found for array index ${PBS_ARRAY_INDEX}"
    exit 1
fi

PROTEINS="${BASEDIR}/${SAMPLE}/${SAMPLE}_proteins.faa"
OUTDIR="${RESULTS}/${SAMPLE}"

mkdir -p "${OUTDIR}"

source ~/.bashrc
conda activate eggnog212

echo "=========================================="
echo "Sample: ${SAMPLE}"
echo "Host:   $(hostname)"
echo "Start:  $(date)"
echo "=========================================="

emapper.py \
    -m diamond \
    --itype proteins \
    --cpu "${PBS_NCPUS:-12}" \
    --data_dir "${EGGNOG_DB}" \
    -i "${PROTEINS}" \
    --output "${SAMPLE}" \
    --output_dir "${OUTDIR}" \
    --override

echo "=========================================="
echo "Finished: ${SAMPLE}"
echo "End: $(date)"
echo "=========================================="
