#!/usr/bin/env python3
"""
STEP 05 - Turn the MMseqs2 95NT output (clustering of the 100NT
representatives) into the 95NT catalogue (as
Urban_Soil/Gene_Catalog/05_postprocess_mmseqs95.py). Every 100NT
representative must appear exactly once.

    05_catalog_95nt.py REP_FASTA CLUSTER_TSV FNA_100NT PREFIX THREADS OUT_FNA_XZ OUT_MATCHES

Outputs
    OUT_FNA_XZ   95NT representatives, '>SC.95NT.### <SC.100NT id>'
    OUT_MATCHES  100NT -> 95NT (.tsv.xz): 100NT id, '=' representative / 'R' reduced, SC.95NT id
"""

import sys

from common import AtomicWriter, check, fasta_iter, grouped_id


def main():
    rep_fasta, cluster_tsv, fna_100nt, prefix, threads, out_fna, out_matches = sys.argv[1:]
    threads = int(threads)

    n100 = sum(1 for _ in fasta_iter(fna_100nt))
    print(f"100NT representatives in step 03: {n100:,}", flush=True)

    new_id = {}
    with AtomicWriter(out_fna, threads) as out:
        for header, seq in fasta_iter(rep_fasta):
            rep100 = header.split()[0]
            if rep100 in new_id:
                raise RuntimeError(f"duplicate representative {rep100}")
            new_id[rep100] = grouped_id(prefix, "95NT", len(new_id))
            out.write(f">{new_id[rep100]} {rep100}\n{seq}\n")
    n95 = len(new_id)
    print(f"95NT representatives: {n95:,}", flush=True)

    seen = set()
    n_self = n_reduced = 0
    with open(cluster_tsv) as f, AtomicWriter(out_matches, threads) as out:
        for line in f:
            rep, member = (x.split()[0] for x in line.rstrip("\n").split("\t"))
            if member in seen:
                raise RuntimeError(f"{member} is in two 95NT clusters")
            seen.add(member)
            rel = "=" if member == rep else "R"
            n_self += rel == "="
            n_reduced += rel == "R"
            out.write(f"{member}\t{rel}\t{new_id[rep]}\n")

    print(f"Reduced 100NT sequences: {n_reduced:,}", flush=True)
    check(len(seen) == n100, f"every one of the {n100:,} 100NT sequences is in exactly one 95NT cluster")
    check(n_self == n95, f"{n95:,} representatives in the FASTA = {n_self:,} '=' rows")


if __name__ == "__main__":
    main()
