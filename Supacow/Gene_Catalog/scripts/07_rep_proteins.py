#!/usr/bin/env python3
"""
STEP 07a - Protein sequences of the 95NT representatives, named by their
SC.95NT IDs, for eggNOG-mapper. Each 95NT representative is a 100NT
representative, which is one exact ORF (its own 100NT and 95NT
representative: the ORF Urban_Soil/Gene_Catalog/08_US95_emapper.py copies
annotations from), so its protein comes straight from step 01.

    07_rep_proteins.py FNA_95NT FNA_100NT ORF_FAA OUT_FAA
"""

import sys

from common import AtomicWriter, check, fasta_iter


def main():
    fna95, fna100, orf_faa, out_faa = sys.argv[1:]
    orf_of_100 = {}
    for header, _ in fasta_iter(fna100):
        rep100, orf = header.split()[:2]
        orf_of_100[rep100] = orf
    want = {}
    for header, _ in fasta_iter(fna95):
        rep95, rep100 = header.split()[:2]
        want[orf_of_100[rep100]] = rep95
    print(f"95NT representatives: {len(want):,}", flush=True)
    n = 0
    with AtomicWriter(out_faa) as out:
        for orf, seq in fasta_iter(orf_faa):
            rep95 = want.get(orf)
            if rep95 is not None:
                out.write(f">{rep95}\n{seq}\n")
                n += 1
    check(n == len(want), f"a protein for each of the {len(want):,} 95NT representatives")


if __name__ == "__main__":
    main()
