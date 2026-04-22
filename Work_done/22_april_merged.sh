#!/usr/bin/env bash
set -euo pipefail

BASE="/work/microbiome/shanghai_dogs/resource_generation/MAGs_Onehealth/External_cohorts/Skani_lists"
QUERY="$BASE/SHD_All_MAGs_list.txt"

QUAL="$BASE/Quality_MAGs"

RESULTS="$BASE/Skani_results"
LOGS="$RESULTS/Logs"

THREADS=40

mkdir -p "$RESULTS" "$LOGS" "$QUAL/Skani_Quality_Results"

COHORTS=(
  Coelho_2018_dog
  Wang_2019_dogs
  Yarlagadda_2022_global_dog
  Allaway_2020_dogs
  Liu_2021_Canidae
  Xu_2019_dogs
  Worsley-Tonks_2020_dog
)

echo "===== STEP 3: ALL MAGs vs SHD ====="
for cohort in "${COHORTS[@]}"; do

    REF="$BASE/${cohort}_MAGs_list.txt"
    OUT="$RESULTS/${cohort}_vs_SHD_ani.tsv"
    LOG="$LOGS/${cohort}.log"

    if [ ! -s "$REF" ]; then
        echo "Skipping $REF"
        continue
    fi

    echo "Running skani: $cohort (ALL MAGs)"

    skani dist \
        --ql "$QUERY" \
        --rl "$REF" \
        -t "$THREADS" \
        -o "$OUT" \
        > "$LOG" 2>&1

done


echo "===== STEP 5: HQ & MQ ====="
for cohort in "${COHORTS[@]}"; do
  for quality in HQ MQ; do

    REF="$QUAL/${cohort}_${quality}_list.txt"
    OUT="$QUAL/Skani_Quality_Results/${cohort}_${quality}_ani.tsv"

    if [ ! -s "$REF" ]; then
      echo "Skipping $REF"
      continue
    fi

    echo "Running skani: $cohort $quality"

    skani dist \
      --ql "$QUERY" \
      --rl "$REF" \
      -t "$THREADS" \
      -o "$OUT"

  done
done


echo "===== STEP 6: ALL (quality-filtered lists) ====="
for cohort in "${COHORTS[@]}"; do

    REF="$QUAL/${cohort}_ALL_list.txt"
    OUT="$QUAL/Skani_Quality_Results/${cohort}_ALL_ani.tsv"

    if [ ! -s "$REF" ]; then
        echo "Skipping $REF"
        continue
    fi

    echo "Running skani: $cohort ALL"

    skani dist \
        --ql "$QUERY" \
        --rl "$REF" \
        -t "$THREADS" \
        -o "$OUT"

done

echo "===== ALL SKANI RUNS COMPLETED ====="
