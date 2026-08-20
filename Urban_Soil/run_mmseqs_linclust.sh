cat > /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/run_mmseqs_linclust.sh << 'EOF'
#!/bin/bash -l
#PBS -N mmseqs_linclust_urbansoil
#PBS -l select=1:ncpus=30:mem=500gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/mmseqs_test/mmseqs_linclust.log
#PBS -j oe
cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate kruthi
BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
INPUT="${BASE}/UrbanSoil.ORF.fna"
OUTDIR="${BASE}/mmseqs_test"
TMPDIR="${OUTDIR}/tmp"
mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"
echo "[1/1] Running mmseqs easy-linclust on full ORF dataset..."
echo "Input: ${INPUT}"
echo "Output dir: ${OUTDIR}"
echo "Started: $(date)"
mmseqs easy-linclust "${INPUT}" "${OUTDIR}/US.mmseqs.100NT" "${TMPDIR}" \
    --min-seq-id 1.0 \
    -c 1.0 \
    --cov-mode 0 \
    --threads 30
echo "Finished: $(date)"
echo "DONE."
echo "Outputs in: ${OUTDIR}"
ls -lh "${OUTDIR}"
EOF
