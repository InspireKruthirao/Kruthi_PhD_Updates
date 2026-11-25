#!/usr/bin/env python3
"""
STEP 06 - ORF -> 100NT -> 95NT cluster-mapping table (as
Urban_Soil/Gene_Catalog/06_make_US_cluster_table.py). Small enough to
join in memory.

    06_cluster_table.py MATCHES_100NT MATCHES_95NT ORIG_TSV THREADS OUT

Output (.tsv.xz), 5 columns per ORF:
    ORF, '='/'C' to 100NT, SC.100NT id, '='/'R' to 95NT, SC.95NT id
"""

import sys

from common import AtomicWriter, check, open_text


def main():
    m100, m95, orig_tsv, threads, out_path = sys.argv[1:]

    to95 = {}
    with open_text(m95) as f:
        for line in f:
            rep100, rel, rep95 = line.split()
            to95[rep100] = (rel, rep95)
    print(f"100NT -> 95NT mappings: {len(to95):,}", flush=True)

    with open_text(orig_tsv) as f:
        n_orfs = sum(1 for _ in f) - 1

    n = 0
    with open_text(m100) as f, AtomicWriter(out_path, int(threads)) as out:
        for line in f:
            orf, rel100, rep100 = line.split()
            rel95, rep95 = to95[rep100]
            out.write(f"{orf}\t{rel100}\t{rep100}\t{rel95}\t{rep95}\n")
            n += 1
    check(n == n_orfs, f"{n:,} rows = {n_orfs:,} ORFs")


if __name__ == "__main__":
    main()
