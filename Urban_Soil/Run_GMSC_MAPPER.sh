#!/bin/bash -l
#PBS -N gmsc_mapper_test
#PBS -l select=1:ncpus=6:mem=80gb
#PBS -l walltime=15:00:00
#PBS -m abe
#PBS -q cpu_batch

cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate gmscmapper

# paths
input="/work/microbiome/urban_soil/data/UrbanSoilAssemblies/snj16_medaka_polypolish.fasta.PolcaCorrected.fa.gz"
output="/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil"
GMSC="/home/n12228516/tools/GMSC-mapper_git_version"
logfile="gmsc_errors.log"
summaryfile="gmsc_summary_report.log"

echo "GMSC Mapper error log - $(date)" > "$logfile"
echo "GMSC Mapper summary report - $(date)" > "$summaryfile"

filename=$(basename "$input")
filename_no_ext="${filename%.fasta.PolcaCorrected.fa.gz}"
outdir="$output/$filename_no_ext"
mkdir -p "$outdir"

echo "🚀 Running gmsc-mapper on $filename"
if ! gmsc-mapper -i "$input" -o "$outdir" --dbdir "$GMSC/db" -t 6; then
    echo "GMSC-mapper execution failed for $filename at $(date)" >> "$logfile"
    echo "$filename: FAILED during execution" >> "$summaryfile"
else
    echo "$filename: Processed successfully" >> "$summaryfile"
fi
