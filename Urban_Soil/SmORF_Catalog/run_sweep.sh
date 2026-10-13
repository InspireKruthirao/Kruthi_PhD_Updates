#!/bin/bash
# Runs the 150 rules (Pixi_Env/sweep_db.py) on smorf_neighbourhoods.sqlite.
# Results go to Neighbourhood_Analysis_SQL/sweep_results/: combinations_results.csv,
# results_by_bin.csv, predictions_all_smorfs.tsv and sweep.log. About 10 min;
# run it inside tmux:   tmux new -s sweep   then   bash run_sweep.sh

set -eo pipefail
PIXI_ENV=$(cd "$(dirname "$0")/../Pixi_Env" && pwd)
NB=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL
DB=$NB/smorf_neighbourhoods.sqlite

mkdir -p "$NB/sweep_results"
cd "$NB/sweep_results"
pixi run --manifest-path "$PIXI_ENV/pixi.toml" python "$PIXI_ENV/sweep_db.py" "$DB" 2>&1 | tee sweep.log
