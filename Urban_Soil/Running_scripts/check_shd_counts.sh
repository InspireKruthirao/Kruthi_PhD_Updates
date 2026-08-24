cat > check_shd_counts.sh << 'EOF'
#!/bin/bash -l
#PBS -N check_shd_counts
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=01:00:00
#PBS -q cpu_batch_exec
#PBS -o /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/check_shd_counts.log
#PBS -j oe

SHD_DIR="/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog"

echo "ORF count:"
xzcat "${SHD_DIR}/SHD.ORF.fna.xz" | grep -c '^>'

echo "100NT count:"
xzcat "${SHD_DIR}/SHD.100NT.fna.xz" | grep -c '^>'

echo "95NT count:"
xzcat "${SHD_DIR}/SHD.95NT.fna.xz" | grep -c '^>'
EOF
qsub check_shd_counts.sh
