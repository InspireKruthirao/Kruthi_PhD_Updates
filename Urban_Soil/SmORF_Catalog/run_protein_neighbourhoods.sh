#!/bin/bash
# Normal proteins instead of smORFs: builds protein_neighbourhoods.sqlite in a
# NEW folder (the smORF databases are not touched), then runs the 236 rules
# on it (sweep_db.py) and the protein -> smORF tests (protein_tests.py).
# Three chained PBS jobs, each starting only if the one before finished:
#   1. chunks   N pieces as one array job (protein_neighbourhoods_sqlite.py)
#   2. merge    the N pieces -> protein_neighbourhoods.sqlite
#   3. tests    sweep_db.py on the proteins, then protein_tests.py
#     bash run_protein_neighbourhoods.sh        (N = 50)
# Results: Neighbourhood_Analysis_SQL/Protein_Neighbourhoods/

set -eo pipefail
N=${1:-50}
CODE=$(cd "$(dirname "$0")" && pwd)
PIXI_ENV=$(cd "$CODE/../Pixi_Env" && pwd)
PY="pixi run --manifest-path $PIXI_ENV/pixi.toml python -u"
NB=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL
OUT=$NB/Protein_Neighbourhoods
SMORF_RULES=$NB/sweep_results_236rules/combinations_results.csv
QUEUE=cpu_batch_exec
RPID=DFAZCB7230
LOGS=$OUT/logs
mkdir -p "$LOGS"

CHUNKS=$(qsub -N prot_chunk -P $RPID -q $QUEUE -J 0-$((N - 1)) \
    -l select=1:ncpus=1:mem=8gb -l walltime=12:00:00 -j oe -o "$LOGS/" << JOB
#!/bin/bash
set -eo pipefail
cd "$CODE"
i=\$PBS_ARRAY_INDEX
if [ -e "$OUT/parts/chunk_\${i}_of_$N.sqlite" ]; then echo "piece \$i already done"; exit 0; fi
$PY protein_neighbourhoods_sqlite.py chunk \$i $N > "$LOGS/chunk_\${i}_of_$N.log" 2>&1
JOB
)
echo "chunks: $CHUNKS ($N pieces)"

MERGE=$(qsub -N prot_merge -P $RPID -q $QUEUE -W depend=afterok:$CHUNKS \
    -l select=1:ncpus=1:mem=16gb -l walltime=6:00:00 -j oe -o "$LOGS/merge.log" << JOB
#!/bin/bash
set -eo pipefail
cd "$CODE"
grep -h DONE "$LOGS"/chunk_*_of_$N.log | sort -V
$PY protein_neighbourhoods_sqlite.py merge $N 2>&1 | grep -v "^ WARN"
JOB
)
echo "merge:  $MERGE"

TESTS=$(qsub -N prot_tests -P $RPID -q $QUEUE -W depend=afterok:$MERGE \
    -l select=1:ncpus=16:mem=180gb -l walltime=8:00:00 -j oe -o "$LOGS/tests.log" << JOB
#!/bin/bash
set -eo pipefail
mkdir -p "$OUT/sweep_results" "$OUT/protein_tests"
cd "$OUT/sweep_results"
$PY "$PIXI_ENV/sweep_db.py" "$OUT/protein_neighbourhoods.sqlite" --save 2>&1 | grep -v "^ WARN" > sweep.log
$PY "$PIXI_ENV/protein_tests.py" "$OUT/protein_neighbourhoods.sqlite" "$NB/smorf_neighbourhoods.sqlite" \
    "$OUT/protein_tests" --protein-rules "$OUT/sweep_results/combinations_results.csv" \
    --smorf-rules "$SMORF_RULES" 2>&1 | grep -v "^ WARN" > "$OUT/protein_tests/protein_tests.log"
JOB
)
echo "tests:  $TESTS"
echo
echo "Watch:   qstat -u \$USER"
echo "Results: $OUT/protein_tests/protein_tests_report.txt"
