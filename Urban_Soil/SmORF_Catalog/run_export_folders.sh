#!/bin/bash
# Rewrites the old neighbourhood folders (Neighbourhood_Analysis_SQL_ge6) with
# the corrected neighbors.tsv, EggNOG_annotations.tsv and <contig>.gff from
# smorf_neighbourhoods.sqlite, as N pieces running at the same time.
#     bash run_export_folders.sh        (N = 100)
# Safe to run again: finished pieces are skipped.

set -eo pipefail
N=${1:-100}
PY=/home/n12228516/.conda/envs/kruthi/bin/python3
CODE=$(cd "$(dirname "$0")" && pwd)
OUT=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL
LOGS=$OUT/logs
mkdir -p "$LOGS"

JOB=$(qsub -N nb_folders -q cpu_batch_exec -J 0-$((N - 1)) \
    -l select=1:ncpus=1:mem=4gb -l walltime=24:00:00 -j oe -o "$LOGS/" << SCRIPT
#!/bin/bash
set -eo pipefail
export PYTHONNOUSERSITE=1
cd "$CODE"
i=\$PBS_ARRAY_INDEX
if [ -e "$LOGS/folders_\${i}_of_$N.done" ]; then echo "piece \$i already done"; exit 0; fi
$PY -s -u build_smorf_neighbourhoods.py export-folders \$i $N > "$LOGS/folders_\${i}_of_$N.log" 2>&1
touch "$LOGS/folders_\${i}_of_$N.done"
SCRIPT
)
echo "folders: $JOB ($N pieces)"
echo "Check:   grep -h DONE $LOGS/folders_*_of_$N.log | wc -l     (should reach $N)"
