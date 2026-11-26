#!/bin/bash -l
#PBS -N urban_orf_catalog
#PBS -l select=1:ncpus=4:mem=20gb
#PBS -l walltime=06:00:00
#PBS -q cpu_batch_exec
#PBS -j oe
#PBS -m abe

set -euo pipefail

cd /work/microbiome/users/kruthi/intermediate_results

echo "Urban Soil ORF catalogue (Stage 1)"
echo "Started: $(date)"

python3 /work/microbiome/users/kruthi/intermediate_results/urban_soil/Gene_Catalog/01_Build_urban_soil_orf_catalog.py

echo "Finished: $(date)"
