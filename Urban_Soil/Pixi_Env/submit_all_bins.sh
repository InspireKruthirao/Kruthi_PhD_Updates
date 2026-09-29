#!/bin/bash
# Submits one PBS job per SmORF_neighbourhoods_* bin, all at once, so they
# run concurrently instead of one script working through every bin serially.
#
# Usage: bash submit_all_bins.sh
# Run it from anywhere -- it discovers bins from PARENT_DIR below.

set -eo pipefail

PARENT_DIR="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL_ge6"
PBS_SCRIPT="$PARENT_DIR/run_cog_vote_bin.pbs"
LOG_DIR="$PARENT_DIR/logs"

mkdir -p "$LOG_DIR"

shopt -s nullglob
bins=("$PARENT_DIR"/SmORF_neighbourhoods_*/)
shopt -u nullglob

if [ ${#bins[@]} -eq 0 ]; then
    echo "No SmORF_neighbourhoods_* folders found under $PARENT_DIR"
    exit 1
fi

echo "Found ${#bins[@]} bins. Submitting one job per bin..."

for bin_path in "${bins[@]}"; do
    bin_path="${bin_path%/}"          # strip trailing slash
    bin_name="$(basename "$bin_path")"
    job_id=$(qsub -N "cog_${bin_name}" \
                   -o "$LOG_DIR/${bin_name}.log" \
                   -v BIN_DIR="$bin_path" \
                   "$PBS_SCRIPT")
    echo "  $bin_name -> $job_id"
done

echo ""
echo "All jobs submitted. Watch progress live with:"
echo "  tail -f $LOG_DIR/*.log"
echo ""
echo "Check job status with:"
echo "  qstat -u \$USER"
