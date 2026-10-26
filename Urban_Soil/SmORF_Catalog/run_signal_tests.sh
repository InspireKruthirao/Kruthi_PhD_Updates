#!/bin/bash
# Is the neighbourhood signal real, and does it fade with distance?
#   Pixi_Env/permutation_test.py  3 rules on the real neighbourhoods against N randomised ones
#   Pixi_Env/distance_decay.py    hit rate by rank, strand and bp distance; context size
# Run inside tmux:
#     bash run_signal_tests.sh [db] [N=100]
# Results: Neighbourhood_Analysis_SQL/signal_tests/ (permutation_test.tsv,
# permutation_null.tsv, distance_decay.tsv, distance_decay.png/.pdf,
# signal_tests.log). Needs ~8 GB of memory; ~15-30 min.

set -eo pipefail
PIXI_ENV=$(cd "$(dirname "$0")/../Pixi_Env" && pwd)
NB=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL
DB=$NB/$(basename "${1:-smorf_neighbourhoods.sqlite}")
[ -f "$DB" ] || { echo "No such database: $DB"; exit 1; }
TAG=$(basename "$DB" .sqlite)
RESULTS=$NB/signal_tests${TAG#smorf_neighbourhoods}

mkdir -p "$RESULTS"
cd "$RESULTS"
echo "Database: $DB"
echo "Results:  $RESULTS"
pixi run --manifest-path "$PIXI_ENV/pixi.toml" python "$PIXI_ENV/distance_decay.py" "$DB" 2>&1 | tee signal_tests.log
pixi run --manifest-path "$PIXI_ENV/pixi.toml" python "$PIXI_ENV/permutation_test.py" "$DB" "${2:-100}" 2>&1 | tee -a signal_tests.log
