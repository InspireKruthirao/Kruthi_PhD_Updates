#!/bin/bash -l
#PBS -N mmseqs_linclust_urbansoil_splitmem
#PBS -l select=1:ncpus=100:mem=700gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/mmseqs_test/mmseqs_linclust_splitmem.log
#PBS -j oe
cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate kruthi
set -euo pipefail
BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
INPUT="${BASE}/UrbanSoil.ORF.fna"
OUTDIR="${BASE}/mmseqs_test"
TMPDIR="${OUTDIR}/tmp_splitmem"
mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"
echo "[1/1] Running mmseqs easy-linclust on full ORF dataset (with --split-memory-limit)..."
echo "Input: ${INPUT}"
echo "Output dir: ${OUTDIR}"
echo "Started: $(date)"
mmseqs easy-linclust "${INPUT}" "${OUTDIR}/US.mmseqs.100NT.splitmem" "${TMPDIR}" \
    --min-seq-id 1.0 \
    -c 1.0 \
    --cov-mode 0 \
    --threads 100 \
    --split-memory-limit 600G
echo "Finished: $(date)"
echo "DONE."
echo "Outputs in: ${OUTDIR}"
ls -lh "${OUTDIR}"
