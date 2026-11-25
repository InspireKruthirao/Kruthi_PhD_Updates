#!/usr/bin/env python3
"""
GMSC-mapper input: the MAGs of one chunk (every N_CHUNKS-th MAG in sorted
order, starting at CHUNK) in one FASTA, contigs renamed <MAG>|<contig> so
each smORF can be traced back to its MAG.

    gmsc_chunk_input.py MAGS_DIR N_CHUNKS CHUNK OUT_FNA
"""

import sys
from pathlib import Path

from common import AtomicWriter, check, fasta_iter


def main():
    mags_dir, n_chunks, chunk, out_fna = sys.argv[1:]
    mags = sorted(p for p in Path(mags_dir).glob("*.fasta"))[int(chunk)::int(n_chunks)]
    n = 0
    with AtomicWriter(out_fna) as out:
        for p in mags:
            mag = p.name[:-len(".fasta")]
            if "|" in mag:
                raise ValueError(f"MAG name contains '|': {mag}")
            for header, seq in fasta_iter(p):
                out.write(f">{mag}|{header.split()[0]}\n{seq}\n")
                n += 1
    check(len(mags) > 0, f"chunk {chunk}: {len(mags)} MAGs, {n:,} contigs")


if __name__ == "__main__":
    main()
