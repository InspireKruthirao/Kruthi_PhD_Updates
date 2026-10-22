#!/bin/bash
# Runs the 188 rules (Pixi_Env/sweep_db.py) on a neighbourhood database, then
# the held-out check and cascade (Pixi_Env/cascade.py). Run inside tmux:
#     bash run_sweep.sh                                     smorf_neighbourhoods.sqlite
#     bash run_sweep.sh smorf_neighbourhoods_q50_s50_p30_e1e-10.sqlite   a filtered copy
# Results: Neighbourhood_Analysis_SQL/sweep_results/ for the full database,
# sweep_results_<filter>/ for a filtered copy: combinations_results.csv,
# results_by_bin.csv, predictions_all_smorfs.tsv, cascade_results.tsv,
# cascade_predictions.tsv and sweep.log. About 15-20 min.

set -eo pipefail
PIXI_ENV=$(cd "$(dirname "$0")/../Pixi_Env" && pwd)
NB=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL
NAME=${1:-smorf_neighbourhoods.sqlite}
DB=$NB/$(basename "$NAME")
[ -f "$DB" ] || { echo "No such database: $DB"; exit 1; }
TAG=$(basename "$DB" .sqlite)
RESULTS=$NB/sweep_results${TAG#smorf_neighbourhoods}

mkdir -p "$RESULTS"
cd "$RESULTS"
echo "Database: $DB"
echo "Results:  $RESULTS"
pixi run --manifest-path "$PIXI_ENV/pixi.toml" python "$PIXI_ENV/sweep_db.py" "$DB" 2>&1 | tee sweep.log
pixi run --manifest-path "$PIXI_ENV/pixi.toml" python "$PIXI_ENV/cascade.py" . 2>&1 | tee -a sweep.log
