#!/bin/bash
# Run the whole Supacow gene catalogue in ONE queue job: Snakemake runs every
# rule inside that job (2,444 small Prodigal runs would flood the queue as
# separate jobs). Safe to resubmit: finished steps are skipped.
#
#   bash run_pipeline.sh            submit (from inside mqyolo, via mqsub)
#   bash run_pipeline.sh --dry-run  only show what would run
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=$(grep '^out_dir:' "$HERE/config.yaml" | awk '{print $2}')
CPUS=32
mkdir -p "$OUT"

CMD="snakemake -s $HERE/Snakefile --configfile $HERE/config.yaml --directory $OUT --cores $CPUS --resources mem_mb=380000 --rerun-incomplete --keep-going --printshellcmds"

if [ "${1:-}" = "--dry-run" ]; then
    $CMD --dry-run --quiet
    exit
fi

mqsub --no-email --bg -t $CPUS -m 400 --days 4 -- bash -c "cd $OUT && $CMD > $OUT/pipeline_\$(date +%Y%m%d_%H%M).log 2>&1"
echo "submitted; follow progress in $OUT/pipeline_*.log and $OUT/logs/"
