cat > /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/check_integrity.sh << 'EOF'
#!/bin/bash -l
#PBS -N check_integrity
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=01:00:00
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/check_integrity.log
#PBS -j oe

cd "$PBS_O_WORKDIR"

echo "Testing integrity of US.100NT.fna.xz..."
xz -t US.100NT.fna.xz
echo "Exit code: $?"
echo "If exit code is 0, file is NOT corrupted."
EOF
qsub check_integrity.sh
