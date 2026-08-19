#!/bin/bash
set -uo pipefail

source ~/.bashrc
conda activate kruthi

BASE="/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
ORF_FULL="${BASE}/UrbanSoil.ORF.fna.xz"
ORF_DECOMP="${BASE}/UrbanSoil.ORF.fna"
CDHIT_OUT="${BASE}/US.100NT_dedup"

echo "[1/3] Decompressing full ORF file (this will take a while, ~186GB of sequence data)..."
if [ ! -f "${ORF_DECOMP}" ]; then
    xzcat "${ORF_FULL}" > "${ORF_DECOMP}"
else
    echo "Decompressed file already exists, skipping."
fi
echo "Decompression done: $(ls -lh ${ORF_DECOMP})"

echo "[2/3] Running cd-hit-est on full dataset..."
cd-hit-est -c 1.0 -aS 1.0 -G 0 -d 0 -g 1 -T 16 -M 0 \
    -i "${ORF_DECOMP}" -o "${CDHIT_OUT}"

echo "[3/3] Post-processing into pipeline-standard output files..."
python3 /work/microbiome/users/kruthi/intermediate_results/postprocess_cdhit_clusters.py \
    "${CDHIT_OUT}" "${CDHIT_OUT}.clstr" \
    "${BASE}/US.100NT.fna.xz" "${BASE}/US.100NT.matches.xz"

echo "DONE."
echo "Final outputs:"
echo "  ${BASE}/US.100NT.fna.xz"
echo "  ${BASE}/US.100NT.matches.xz"
