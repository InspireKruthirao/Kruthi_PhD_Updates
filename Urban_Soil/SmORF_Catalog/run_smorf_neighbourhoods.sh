#!/bin/bash
# Submits the whole smorf_neighbourhoods.sqlite run as three chained PBS jobs:
#   1. build   urban_soil.db (finishes in seconds if it is already complete)
#   2. chunks  N pieces as one array job, all running at the same time
#   3. merge   the N pieces -> smorf_neighbourhoods.sqlite
# Each job starts only if the one before finished without errors.
#     bash run_smorf_neighbourhoods.sh        (N = 100)
#     bash run_smorf_neighbourhoods.sh 50     (N = 50)
# Safe to run again: finished pieces are skipped.

set -eo pipefail
N=${1:-100}
PY=/home/n12228516/.conda/envs/kruthi/bin/python3
CODE=$(cd "$(dirname "$0")" && pwd)
OUT=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL
QUEUE=cpu_batch_exec
LOGS=$OUT/logs
mkdir -p "$LOGS"

BUILD=$(qsub -N nb_build -q $QUEUE -l select=1:ncpus=1:mem=32gb -l walltime=24:00:00 \
    -j oe -o "$LOGS/build.log" << JOB
#!/bin/bash
set -eo pipefail
export PYTHONNOUSERSITE=1
cd "$CODE"
$PY -s -u smorf_neighbourhoods.py build
JOB
)
echo "build:  $BUILD"

CHUNKS=$(qsub -N nb_chunk -q $QUEUE -J 0-$((N - 1)) -W depend=afterok:$BUILD \
    -l select=1:ncpus=1:mem=8gb -l walltime=24:00:00 -j oe -o "$LOGS/" << JOB
#!/bin/bash
set -eo pipefail
export PYTHONNOUSERSITE=1
cd "$CODE"
i=\$PBS_ARRAY_INDEX
if [ -e "$OUT/parts/chunk_\${i}_of_$N.sqlite" ]; then echo "piece \$i already done"; exit 0; fi
$PY -s -u smorf_neighbourhoods.py chunk \$i $N > "$LOGS/chunk_\${i}_of_$N.log" 2>&1
JOB
)
echo "chunks: $CHUNKS ($N pieces)"

MERGE=$(qsub -N nb_merge -q $QUEUE -W depend=afterok:$CHUNKS \
    -l select=1:ncpus=1:mem=16gb -l walltime=24:00:00 -j oe -o "$LOGS/merge.log" << JOB
#!/bin/bash
set -eo pipefail
export PYTHONNOUSERSITE=1
cd "$CODE"
grep -h DONE "$LOGS"/chunk_*_of_$N.log | sort -V
$PY -s -u smorf_neighbourhoods.py merge $N
JOB
)
echo "merge:  $MERGE"
echo
echo "Watch:   qstat -u \$USER      and      qstat -t $CHUNKS | tail"
echo "Result:  $LOGS/merge.log"
