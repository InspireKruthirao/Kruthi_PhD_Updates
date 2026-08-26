cat > /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/run_postprocess_95nt.sh << 'EOF'
#!/bin/bash -l
#PBS -N postprocess_95nt
#PBS -l select=1:ncpus=16:mem=128gb
#PBS -l walltime=24:00:00
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/postprocess_95nt.log
#PBS -j oe

cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate kruthi

python3 postprocess_mmseqs_95nt.py
EOF
qsub /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/run_postprocess_95nt.sh
