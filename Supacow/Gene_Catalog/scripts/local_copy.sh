#!/bin/bash
# local_copy.sh SRC NAME -> prints the path of a copy of directory SRC on the
# node's local disk ($TMPDIR, else /tmp), made once per job: parallel rules
# share it (flock), later ones reuse it. eResearch rule: databases are not
# read from /work during jobs.
set -euo pipefail
SRC="$1"; NAME="$2"
BASE="${TMPDIR:-/tmp}"
DEST="$BASE/$NAME"
exec 9>"$BASE/$NAME.lock"
flock 9
if [ ! -f "$DEST/.copied" ]; then
    rm -rf "$DEST"
    need=$(find "$SRC" -type f -printf "%s\n" | awk '{s+=$1} END {print int(s/1024)}')  # du is wrong on home
    free=$(df -Pk "$BASE" | awk 'NR==2 {print $4}')
    echo "copying $SRC ($need kB) to $DEST ($free kB free): $(date)" >&2
    [ "$free" -gt $((need + 20000000)) ] || { echo "not enough local disk in $BASE" >&2; exit 1; }
    mkdir -p "$DEST"
    cp -r "$SRC"/. "$DEST"/
    touch "$DEST/.copied"
    echo "copied: $(date)" >&2
fi
echo "$DEST"
