#!/bin/bash -l
#PBS -N gmsc_mapper_all
#PBS -l select=1:ncpus=60:mem=700gb
#PBS -l walltime=26:00:00
#PBS -m abe
#PBS -q cpu_batch

cd "$PBS_O_WORKDIR"
source ~/.bashrc
conda activate gmscmapper

INPUT_DIR="/work/microbiome/urban_soil/data/UrbanSoilAssemblies"
OUTPUT_DIR="/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil"
GMSC="/home/n12228516/tools/GMSC-mapper_git_version"
LOGDIR="$OUTPUT_DIR/logs"
LOGFILE="$LOGDIR/gmsc_master_$(date +%Y%m%d_%H%M%S).log"
SUMMARYFILE="$LOGDIR/gmsc_summary_$(date +%Y%m%d_%H%M%S).log"
MAX_PARALLEL=10

mkdir -p "$LOGDIR"

echo "========================================" | tee -a "$LOGFILE"
echo "GMSC Mapper batch run - $(date)"         | tee -a "$LOGFILE"
echo "Input dir : $INPUT_DIR"                  | tee -a "$LOGFILE"
echo "Output dir: $OUTPUT_DIR"                 | tee -a "$LOGFILE"
echo "Max parallel: $MAX_PARALLEL"             | tee -a "$LOGFILE"
echo "========================================" | tee -a "$LOGFILE"

run_sample() {
    local input="$1"
    local filename
    filename=$(basename "$input")
    local sample="${filename%.fasta.PolcaCorrected.fa.gz}"
    local outdir="$OUTPUT_DIR/$sample"
    local samplelog="$LOGDIR/${sample}.log"
    local marker="$outdir/done.txt"

    if [[ -f "$marker" ]]; then
        echo "[SKIP] $sample — already done" | tee -a "$LOGFILE"
        echo "$sample: SKIPPED (already done)" >> "$SUMMARYFILE"
        return 0
    fi

    mkdir -p "$outdir"
    echo "[START] $sample — $(date)" | tee -a "$LOGFILE"

    if gmsc-mapper \
        -i "$input" \
        -o "$outdir" \
        --dbdir "$GMSC/db" \
        -t 6 \
        --tool diamond \
        > "$samplelog" 2>&1; then
        touch "$marker"
        echo "[DONE]  $sample — $(date)" | tee -a "$LOGFILE"
        echo "$sample: SUCCESS" >> "$SUMMARYFILE"
    else
        echo "[FAIL]  $sample — $(date) — see $samplelog" | tee -a "$LOGFILE"
        echo "$sample: FAILED — check $samplelog" >> "$SUMMARYFILE"
    fi
}

export -f run_sample
export OUTPUT_DIR GMSC LOGDIR LOGFILE SUMMARYFILE

running=0
pids=()

for input in "$INPUT_DIR"/*.fasta.PolcaCorrected.fa.gz; do
    run_sample "$input" &
    pids+=($!)
    (( running++ ))

    if (( running >= MAX_PARALLEL )); then
        wait "${pids[0]}"
        pids=("${pids[@]:1}")
        (( running-- ))
    fi
done

wait

echo "" | tee -a "$LOGFILE"
echo "========================================" | tee -a "$LOGFILE"
echo "All samples processed — $(date)" | tee -a "$LOGFILE"
echo "Summary:" | tee -a "$LOGFILE"
grep "SUCCESS" "$SUMMARYFILE" | wc -l | xargs -I{} echo "  Succeeded: {}" | tee -a "$LOGFILE"
grep "FAILED"  "$SUMMARYFILE" | wc -l | xargs -I{} echo "  Failed   : {}" | tee -a "$LOGFILE"
grep "SKIPPED" "$SUMMARYFILE" | wc -l | xargs -I{} echo "  Skipped  : {}" | tee -a "$LOGFILE"
echo "Full log : $LOGFILE" | tee -a "$LOGFILE"
echo "Summary  : $SUMMARYFILE" | tee -a "$LOGFILE"
echo "========================================" | tee -a "$LOGFILE"
