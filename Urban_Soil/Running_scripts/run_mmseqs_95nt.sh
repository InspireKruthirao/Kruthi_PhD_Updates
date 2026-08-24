cat > /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/run_mmseqs_95nt.sh << 'EOF'
#!/bin/bash -l
#PBS -N mmseqs_95nt_urbansoil
#PBS -l select=1:ncpus=100:mem=700gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/mmseqs_95nt/mmseqs_95nt.log
#PBS -j oe

cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate kruthi
set -euo pipefail

BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
INPUT_XZ="${BASE}/US.100NT.fna.xz"
INPUT_FNA="${BASE}/mmseqs_95nt/US.100NT.fna"
OUTDIR="${BASE}/mmseqs_95nt"
TMPDIR="${OUTDIR}/tmp"

mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"

echo "Decompressing input (${INPUT_XZ}) for mmseqs..."
echo "Started: $(date)"
xzcat "${INPUT_XZ}" > "${INPUT_FNA}"
echo "Decompression done: $(ls -lh ${INPUT_FNA})"

echo "Running mmseqs easy-cluster at 95% identity, 90% coverage..."
mmseqs easy-cluster "${INPUT_FNA}" "${OUTDIR}/US.mmseqs.95NT" "${TMPDIR}" \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 0 \
    --threads 100 \
    --split-memory-limit 600G

echo "Finished: $(date)"
echo "DONE."
echo "Outputs in: ${OUTDIR}"
ls -lh "${OUTDIR}"
EOF
