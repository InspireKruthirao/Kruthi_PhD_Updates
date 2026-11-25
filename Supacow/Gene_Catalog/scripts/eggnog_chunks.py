#!/usr/bin/env python3
"""
Split the 95NT representative proteins into chunks for eggNOG-mapper, or
merge the chunks' annotations back into one file.

    eggnog_chunks.py split FAA OUT_DIR N_CHUNKS     -> OUT_DIR/chunkNN.faa
    eggnog_chunks.py merge OUT_GZ ANNOTATIONS...    -> one .emapper.annotations.gz
"""

import gzip
import os
import sys

from common import AtomicWriter, check, fasta_iter


def split(faa, out_dir, n):
    n = int(n)
    os.makedirs(out_dir, exist_ok=True)
    outs = [AtomicWriter(os.path.join(out_dir, f"chunk{i:02d}.faa")) for i in range(n)]
    k = 0
    for header, seq in fasta_iter(faa):
        outs[k % n].write(f">{header.split()[0]}\n{seq}\n")
        k += 1
    for o in outs:
        o.__exit__(None, None, None)
    check(k > 0, f"{k:,} proteins in {n} chunks")


def merge(out_gz, files):
    rows = 0
    with AtomicWriter(out_gz) as out:
        for i, path in enumerate(files):
            with open(path) as f:
                for line in f:
                    if line.startswith("#"):
                        # keep the header lines of the first chunk only
                        if i == 0 and (line.startswith("#query") or rows == 0):
                            out.write(line)
                        continue
                    out.write(line)
                    rows += 1
    print(f"annotation rows: {rows:,} from {len(files)} chunks", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "split":
        split(*sys.argv[2:])
    elif sys.argv[1] == "merge":
        merge(sys.argv[2], sys.argv[3:])
    else:
        raise SystemExit(__doc__)
