#!/usr/bin/env python3
"""
STEP 08 - One summary of the catalogue: genes per stage, partial genes,
eggNOG coverage and COG categories of the 95NT representatives.

    08_summary.py ORIG_TSV CLUSTERS EMAPPER_ANNOTATIONS OUT_TSV
"""

import sys
from collections import Counter

from common import open_text


def main():
    orig, clusters, annotations, out_path = sys.argv[1:]
    mags, partial = set(), Counter()
    n_orfs = 0
    with open_text(orig) as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            mags.add(x[1])
            partial[x[6]] += 1
            n_orfs += 1
    rep100, rep95 = set(), set()
    with open_text(clusters) as f:
        for line in f:
            x = line.split()
            rep100.add(x[2])
            rep95.add(x[4])
    n_ann = n_cog = 0
    cog = Counter()
    col = None
    with open_text(annotations) as f:
        for line in f:
            if line.startswith("#query"):
                col = line.lstrip("#").rstrip("\n").split("\t").index("COG_category")
            if line.startswith("#"):
                continue
            n_ann += 1
            c = line.rstrip("\n").split("\t")[col]
            if c not in ("-", ""):
                n_cog += 1
                for letter in c:
                    cog[letter] += 1 / len(c)

    pct = lambda a, b: f"{100 * a / b:.1f}" if b else "0"
    rows = [
        ("MAGs", len(mags), ""),
        ("ORFs (Prodigal genes)", n_orfs, ""),
        ("ORFs complete (partial=00)", partial.get("00", 0), pct(partial.get("00", 0), n_orfs)),
        ("100NT representatives", len(rep100), pct(len(rep100), n_orfs)),
        ("95NT representatives", len(rep95), pct(len(rep95), n_orfs)),
        ("95NT with an eggNOG hit", n_ann, pct(n_ann, len(rep95))),
        ("95NT with a COG category", n_cog, pct(n_cog, len(rep95))),
    ]
    rows += [(f"COG category {k}", round(v), pct(v, n_cog)) for k, v in sorted(cog.items(), key=lambda kv: -kv[1])]
    with open(out_path, "w") as out:
        out.write("measure\tcount\tpercent\n")
        for name, count, p in rows:
            out.write(f"{name}\t{count}\t{p}\n")
            print(f"{name:32s} {count:>12,} {p:>6}", flush=True)


if __name__ == "__main__":
    main()
