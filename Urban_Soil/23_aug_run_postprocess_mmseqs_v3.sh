cat > /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/run_postprocess_mmseqs_v3.sh << 'EOF'
#!/bin/bash -l
#PBS -N postprocess_mmseqs_v3
#PBS -l select=1:ncpus=80:mem=300gb
#PBS -l walltime=48:00:00
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/postprocess_mmseqs_v3.log
#PBS -j oe

cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate kruthi

python3 postprocess_mmseqs_clusters_v3.py
EOF
