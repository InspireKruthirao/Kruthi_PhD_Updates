#!/bin/bash
# Compare old vs new Urban Soil gene catalogue. Reads only; writes report to $OUT.
#   OLD=... NEW=... OUT=... bash run_compare.sh
set -euo pipefail
U=/work/microbiome/users/kruthi/intermediate_results/urban_soil
OLD=${OLD:-$U/Gene_Catalog_Urban_Soil}
OLD_ALT=${OLD_ALT:-$U/Prodigal}          # old US.95NT.matches.tsv.xz lives here
NEW=${NEW:-$U/Gene_Catalog_v2}
OUT=${OUT:-/mnt/hpccs01/home/n12228516/gc_compare}
PY=/home/n12228516/.conda/envs/kruthi/bin/python
SCRIPT=$(cd "$(dirname "$0")" && pwd)/compare_catalogs.py
export PATH=/mnt/weka/pkg/cmr/sw/conda/bin:$PATH   # xz, pigz
: "${TMPDIR:?TMPDIR not set}"
W=$TMPDIR/gc_compare_work; mkdir -p "$W" "$OUT"
echo "start $(date)  work=$W"; df -h "$TMPDIR" | tail -1

oldpath() { if [ -e "$OLD/$1" ]; then echo "$OLD/$1"; else echo "$OLD_ALT/$1"; fi; }

md5job() {  # file side path
  xz -T4 -dc "$3" | tee >(wc -l > "$W/lines.$1.$2") | md5sum | cut -d' ' -f1 > "$W/md5only.$1.$2"
  sleep 2; echo "$(cat "$W/md5only.$1.$2") $(cat "$W/lines.$1.$2")" > "$W/md5.$1.$2"
  echo "md5 $1 $2 done $(date +%T)"
}

pids=()
for side in old new; do
  for f in UrbanSoil.ORF.fna.xz UrbanSoil.ORF.orig.tsv.xz; do
    if [ $side = old ]; then p=$(oldpath $f); else p=$NEW/$f; fi
    md5job $f $side "$p" & pids+=($!)
  done
  for kf in fasta100:US.100NT.fna.xz fasta95:US.95NT.fna.xz clusters:US.clusters.tsv.xz \
            matches100:US.100NT.matches.xz matches95:US.95NT.matches.tsv.xz \
            emapper:US.95NT.emapper.annotations.gz; do
    k=${kf%%:*}; f=${kf#*:}
    if [ $side = old ]; then p=$(oldpath $f); else p=$NEW/$f; fi
    $PY "$SCRIPT" parse $k $side "$p" "$W" & pids+=($!)
  done
done
fail=0; for p in "${pids[@]}"; do wait $p || fail=1; done
[ $fail = 0 ] || { echo "a parse step failed"; exit 1; }
echo "all parsing done $(date)"
$PY "$SCRIPT" report "$W" "$OUT"
echo "finished $(date)"
