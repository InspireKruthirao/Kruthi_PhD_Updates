#!/bin/bash -l
#PBS -N mmseqs_95nt_64bit
#PBS -l select=1:ncpus=100:mem=1200gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/mmseqs_95nt_64bit.log
#PBS -j oe

cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate kruthi
set -euo pipefail

MMSEQS_BIN="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/mmseqs_64bit/mmseqs/bin/mmseqs"

BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
INPUT_XZ="${BASE}/US.100NT_v3.fna.xz"
OUTDIR="${BASE}/mmseqs_95nt_64bit"
INPUT_FNA="${OUTDIR}/US.100NT.fna"
TMPDIR="${OUTDIR}/tmp"

mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"

echo "=========================================="
echo "MMseqs2 95% clustering - 64-bit build test"
echo "=========================================="
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "MMseqs binary: ${MMSEQS_BIN}"
"${MMSEQS_BIN}" version
echo "CPUs: 100"
echo "Memory requested: 1200 GB"
echo "Input: ${INPUT_XZ}"
echo "Output: ${OUTDIR}"
echo "=========================================="
echo ""
echo "Decompressing input..."
echo "Started: $(date)"
xzcat "${INPUT_XZ}" > "${INPUT_FNA}"
echo ""
echo "Decompression completed: $(date)"
ls -lh "${INPUT_FNA}"
echo ""
echo "Running MMseqs2 (64-bit build) easy-cluster..."
"${MMSEQS_BIN}" easy-cluster \
    "${INPUT_FNA}" \
    "${OUTDIR}/US.mmseqs.95NT" \
    "${TMPDIR}" \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 0 \
    --threads 100 \
    --split-memory-limit 1000G
echo ""
echo "=========================================="
echo "MMseqs2 finished: $(date)"
echo "=========================================="
echo ""
echo "Output files:"
ls -lh "${OUTDIR}"
echo ""
echo "DONE."
