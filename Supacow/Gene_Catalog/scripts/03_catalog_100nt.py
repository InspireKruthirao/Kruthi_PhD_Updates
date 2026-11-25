#!/usr/bin/env python3
"""
STEP 03 - Turn the MMseqs2 100NT (exact dedup) output into the 100NT
catalogue (as Urban_Soil/Gene_Catalog/03_Run_US_100NT_catalog.py).
Instead of hard-coded expected counts, every ORF of step 01 must appear
exactly once.

    03_catalog_100nt.py REP_FASTA CLUSTER_TSV ORIG_TSV PREFIX THREADS OUT_FNA OUT_FNA_XZ OUT_MATCHES

Outputs
    OUT_FNA / OUT_FNA_XZ  100NT representatives, '>SC.100NT.### <ORF id>'
                          (plain copy is the input of the 95NT clustering)
    OUT_MATCHES           ORF -> 100NT (.xz): orf, '=' representative / 'C' collapsed, SC.100NT id
"""

import sys

from common import AtomicWriter, check, fasta_iter, grouped_id, open_text


def main():
    rep_fasta, cluster_tsv, orig_tsv, prefix, threads, out_fna, out_fna_xz, out_matches = sys.argv[1:]
    threads = int(threads)

    with open_text(orig_tsv) as f:
        n_orfs = sum(1 for _ in f) - 1
    print(f"ORFs in step 01: {n_orfs:,}", flush=True)

    new_id = {}
    with AtomicWriter(out_fna) as plain, AtomicWriter(out_fna_xz, threads) as xz:
        for header, seq in fasta_iter(rep_fasta):
            orf = header.split()[0]
            if orf in new_id:
                raise RuntimeError(f"duplicate representative {orf}")
            new_id[orf] = grouped_id(prefix, "100NT", len(new_id))
            record = f">{new_id[orf]} {orf}\n{seq}\n"
            plain.write(record)
            xz.write(record)
    n100 = len(new_id)
    print(f"100NT representatives: {n100:,}", flush=True)

    seen = set()
    n_self = n_collapsed = 0
    with open(cluster_tsv) as f, AtomicWriter(out_matches, threads) as out:
        for line in f:
            rep, member = line.rstrip("\n").split("\t")
            if member in seen:
                raise RuntimeError(f"ORF {member} is in two 100NT clusters")
            seen.add(member)
            if member == rep:
                out.write(f"{member}\t=\t{new_id[rep]}\n")
                n_self += 1
            else:
                out.write(f"{member}\tC\t{new_id[rep]}\n")
                n_collapsed += 1

    print(f"Collapsed ORFs: {n_collapsed:,}", flush=True)
    check(len(seen) == n_orfs, f"every one of the {n_orfs:,} ORFs is in exactly one 100NT cluster")
    check(n_self == n100, f"{n100:,} representatives in the FASTA = {n_self:,} '=' rows")


if __name__ == "__main__":
    main()
