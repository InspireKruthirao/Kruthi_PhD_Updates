#!/bin/bash
# Permutation test of the neighbourhood signal (Pixi_Env/permutation_test.py):
# 3 rules on the real neighbourhoods against N randomised ones. Run inside tmux:
#     bash run_permutation_test.sh [db] [N=100]
# Results: Neighbourhood_Analysis_SQL/permutation_test/ (permutation_test.tsv,
# permutation_null.tsv, permutation.log). Needs ~8 GB of memory; ~10-20 min.

set -eo pipefail
PIXI_ENV=$(cd "$(dirname "$0")/../Pixi_Env" && pwd)
NB=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL
DB=$NB/$(basename "${1:-smorf_neighbourhoods.sqlite}")
[ -f "$DB" ] || { echo "No such database: $DB"; exit 1; }
TAG=$(basename "$DB" .sqlite)
RESULTS=$NB/permutation_test${TAG#smorf_neighbourhoods}

mkdir -p "$RESULTS"
cd "$RESULTS"
echo "Database: $DB"
echo "Results:  $RESULTS"
pixi run --manifest-path "$PIXI_ENV/pixi.toml" python "$PIXI_ENV/permutation_test.py" "$DB" "${2:-100}" 2>&1 | tee permutation.log
