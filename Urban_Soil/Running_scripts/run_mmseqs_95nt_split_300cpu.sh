cat > run_mmseqs_95nt_split_300cpu.sh << 'EOF'
#!/bin/bash -l
#PBS -N mmseqs_95nt_split_300cpu
#PBS -l select=1:ncpus=300:mem=400gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/mmseqs_95nt_split_300cpu.log
#PBS -j oe

cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate kruthi
set -euo pipefail

BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
INPUT_XZ="${BASE}/US.100NT_v3.fna.xz"
OUTDIR="${BASE}/mmseqs_95nt_split_300cpu"
INPUT_FNA="${OUTDIR}/US.100NT.fna"
N_CHUNKS=4

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "=========================================="
echo "MMseqs2 95% clustering - split approach, 300 CPUs"
echo "=========================================="
echo "Started: $(date)"
echo "Host: $(hostname)"

echo ""
echo "Decompressing input..."
xzcat "${INPUT_XZ}" > "${INPUT_FNA}"
echo "Decompression done: $(date)"
ls -lh "${INPUT_FNA}"

echo ""
echo "Splitting into ${N_CHUNKS} chunks..."
python3 "${BASE}/split_fasta.py" "${INPUT_FNA}" "${N_CHUNKS}" chunk_300cpu
echo "Split done: $(date)"
ls -lh chunk_300cpu_*.fna

for i in $(seq 0 $((N_CHUNKS - 1))); do
    echo ""
    echo "=========================================="
    echo "Clustering chunk ${i}..."
    echo "Started: $(date)"
    echo "=========================================="
    mkdir -p "tmp_chunk_300cpu_${i}"
    mmseqs easy-cluster \
        "chunk_300cpu_${i}.fna" \
        "chunk_300cpu_${i}_95NT" \
        "tmp_chunk_300cpu_${i}" \
        --min-seq-id 0.95 \
        -c 0.9 \
        --cov-mode 0 \
        --threads 300 \
        --split-memory-limit 350G
    echo "Chunk ${i} done: $(date)"
done

echo ""
echo "All chunks clustered. Merging representative sequences..."
cat chunk_300cpu_*_95NT_rep_seq.fasta > merged_reps_300cpu.fna
echo "Merged: $(wc -l < merged_reps_300cpu.fna) lines"

echo ""
echo "Running final clustering pass on merged representatives..."
mkdir -p tmp_final_300cpu
mmseqs easy-cluster \
    merged_reps_300cpu.fna \
    US.mmseqs.95NT.final_300cpu \
    tmp_final_300cpu \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 0 \
    --threads 300 \
    --split-memory-limit 350G

echo ""
echo "=========================================="
echo "DONE: $(date)"
echo "=========================================="
echo "Final outputs:"
ls -lh US.mmseqs.95NT.final_300cpu*
EOF
qsub run_mmseqs_95nt_split_300cpu.sh
