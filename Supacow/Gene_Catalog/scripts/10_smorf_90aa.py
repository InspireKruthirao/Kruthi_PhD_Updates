#!/usr/bin/env python3
"""
STEP 10 - 90AA smORF clusters from CD-HIT (as
Urban_Soil/SmORF_Catalog/Process_90AA_clusters.py). Clusters are ordered
by how many smORF origins they hold (most first; ties: representative ID)
and get pad9 IDs COW_SM.90AA.000_000_000.

    10_smorf_90aa.py CLSTR FAA_100AA ORIGINS PREFIX OUT_CLUSTERS OUT_FAA

Outputs
    OUT_CLUSTERS  100AA SmORF ID, 90AA SmORF ID, Similarity ('*' = representative)
    OUT_FAA       >COW_SM.90AA.### <representative COW_SM.100AA id>
"""

import re
import sys
from collections import Counter

from common import AtomicWriter, check, fasta_iter, grouped_id, open_text

LINE = re.compile(r"\d+\s+\d+aa,\s+>(.*?)\.\.\.\s*(.*)")


def main():
    clstr, faa100, origins, prefix, out_clusters, out_faa = sys.argv[1:]

    clusters, current = [], None
    with open(clstr) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current = []
                clusters.append(current)
                continue
            m = LINE.match(line)
            if not m:
                continue
            sim = m.group(2).strip()
            if sim == "*":
                current.append((m.group(1), "*"))
            else:
                p = re.search(r"at\s+(\d+\.\d+%)", sim)
                if p:
                    current.append((m.group(1), p.group(1)))

    n_origins = Counter()
    with open_text(origins) as f:
        next(f)
        for line in f:
            n_origins[line.split("\t", 1)[0]] += 1

    def key(members):
        [rep] = [s for s, sim in members if sim == "*"]
        return (-sum(n_origins[s] for s, _ in members), rep)

    clusters.sort(key=key)
    rows, rep_of = [], {}
    for i, members in enumerate(clusters):
        cid = grouped_id(prefix, "90AA", i)
        for s, sim in members:
            rows.append((s, cid, sim))
            if sim == "*":
                rep_of[s] = cid
    rows.sort()

    n100 = sum(1 for _ in fasta_iter(faa100))
    check(len(rows) == n100, f"each of the {n100:,} 100AA smORFs is in a 90AA cluster")
    check(len({r[0] for r in rows}) == len(rows), "no 100AA smORF is in two clusters")
    check(len(rep_of) == len(clusters), f"each of the {len(clusters):,} clusters has one representative")

    with AtomicWriter(out_clusters) as out:
        out.write("100AA SmORF ID\t90AA SmORF ID\tSimilarity\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    seqs = sorted((rep_of[h.split()[0]], h.split()[0], s) for h, s in fasta_iter(faa100) if h.split()[0] in rep_of)
    with AtomicWriter(out_faa) as out:
        for cid, h100, s in seqs:
            out.write(f">{cid} {h100}\n{s}\n")
    print(f"90AA clusters: {len(clusters):,} from {n100:,} 100AA smORFs", flush=True)


if __name__ == "__main__":
    main()
