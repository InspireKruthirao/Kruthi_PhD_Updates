#!/bin/bash
set -uo pipefail

source ~/.bashrc
conda activate kruthi

BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
TESTDIR="${BASE}/test_1M"
mkdir -p "${TESTDIR}"

ORF_FULL="${BASE}/UrbanSoil.ORF.fna.xz"
ORF_1M="${TESTDIR}/UrbanSoil.ORF.1M.fna"

echo "[1/4] Extracting first 1M sequences..."
xzcat "${ORF_FULL}" | awk '/^>/{n++} n>1000000{exit} {print}' > "${ORF_1M}" || true
echo "Extracted: $(grep -c '^>' ${ORF_1M}) sequences"

echo "[2/4] Running original Python dedup approach..."
python3 /work/microbiome/users/kruthi/intermediate_results/redundant100_test.py \
    "${ORF_1M}" "${TESTDIR}/test_python.100NT.fna" "${TESTDIR}/test_python.100NT.matches.tsv"

echo "[3/4] Running cd-hit-est approach..."
cd-hit-est -c 1.0 -aS 1.0 -G 0 -d 0 -g 1 -T 12 -M 18000 \
    -i "${ORF_1M}" -o "${TESTDIR}/test_cdhit.100NT"

echo "[4/4] Comparing outputs..."
python3 /work/microbiome/users/kruthi/intermediate_results/compare_dedup_outputs.py \
    "${TESTDIR}/test_python.100NT.fna" \
    "${TESTDIR}/test_cdhit.100NT" \
    "${TESTDIR}/comparison_report.txt"

echo "DONE. See ${TESTDIR}/comparison_report.txt"
