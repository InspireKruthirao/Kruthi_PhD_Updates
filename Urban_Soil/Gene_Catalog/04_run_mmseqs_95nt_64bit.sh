#!/bin/bash -l
#PBS -N mmseqs_95nt_64bit
#PBS -l select=1:ncpus=100:mem=700gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -q cpu_batch_exec
#PBS -j oe
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/mmseqs_95nt_64bit.log

# STAGE 4 - Cluster the 256.8M 100NT representatives at 95% nucleotide
# identity / 90% coverage ("95NT" level of the catalogue).
#
# NOTE ON RECONSTRUCTION
# -----------------------
# This script was originally authored directly with `vim` on the HPC
# (25-26 Aug) and its exact byte-for-byte content was therefore not
# captured in shell history, unlike the other pipeline scripts. The
# parameters below (--min-seq-id 0.95 -c 0.9 --cov-mode 0, the use of
# a statically-linked mmseqs2 binary, and the mmseqs_95nt_64bit output
# directory name) are reconstructed from the surrounding commands and
# log messages and match the files that 05_postprocess_mmseqs95.py
# actually consumed (US.mmseqs.95NT_rep_seq.fasta /
# US.mmseqs.95NT_cluster.tsv). Re-check the resource request
# (ncpus/mem) before resubmitting.
#
# Several earlier attempts (700gb / 300gb single-node runs using the
# conda-installed mmseqs2) failed with illegal-instruction crashes on
# some compute nodes, apparently because that build required AVX-512.
# The fix was to download a portable AVX2 static build instead:
#
#   wget https://mmseqs.com/latest/mmseqs-linux-avx2-64bit.tar.gz
#   tar xvfz mmseqs-linux-avx2-64bit.tar.gz -C mmseqs_64bit
#
# Input : US.mmseqs.100NT.splitmem_rep_seq.fasta
#         (the 100NT representative sequences from Stage 2/3)
# Output: mmseqs_95nt_64bit/US.mmseqs.95NT_rep_seq.fasta
#         mmseqs_95nt_64bit/US.mmseqs.95NT_cluster.tsv

cd "$PBS_O_WORKDIR"
set -euo pipefail

BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"

# Statically-linked AVX2 mmseqs2 binary (works around illegal
# instruction crashes seen with the conda build on some nodes).
MMSEQS_BIN="${BASE}/mmseqs_64bit/mmseqs/bin/mmseqs"

INPUT="${BASE}/mmseqs_test/US.mmseqs.100NT.splitmem_rep_seq.fasta"
OUTDIR="${BASE}/mmseqs_95nt_64bit"
TMPDIR="${OUTDIR}/tmp"

mkdir -p "${OUTDIR}" "${TMPDIR}"

echo "Urban Soil 95NT clustering of 100NT representatives"
echo "mmseqs binary: ${MMSEQS_BIN}"
echo "Input:  ${INPUT}"
echo "Output: ${OUTDIR}"
echo "Started: $(date)"

"${MMSEQS_BIN}" easy-cluster "${INPUT}" "${OUTDIR}/US.mmseqs.95NT" "${TMPDIR}" \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 0 \
    --threads 100 \
    --split-memory-limit 600G

echo "Finished: $(date)"
ls -lh "${OUTDIR}"
