#!/usr/bin/env python3
"""
How much KEGG annotation is there in the Urban Soil data, compared with COG?

eggNOG-mapper already gives every annotated gene its KEGG annotations
(KEGG_ko, KEGG_Pathway, KEGG_Module, BRITE), stored in urban_soil.db with
the full eggNOG lines. Before using KEGG as labels instead of COG
categories, this reports:

  1. genes with eggNOG (the first SAMPLE lines; all would take hours): the %
     with a COG category, KO, KEGG pathway, KEGG module and BRITE
  2. the neighbour genes and the target genes in smorf_neighbourhoods.sqlite:
     the same %
  3. the smORFs: how many have a label (the target gene of at least one
     occurrence has one), i.e. how many could be evaluated with each
     label, and how many different labels there are

Pathways are counted as map IDs (eggNOG lists each pathway twice, as koNNNNN
and mapNNNNN). Only reads the two databases; writes nothing.

    python3 kegg_analysis.py [urban_soil.db] [smorf_neighbourhoods.sqlite]
"""

import sqlite3
import sys
import time
from collections import Counter, defaultdict

NB_DIR = "/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL"
LABELS = ("COG_category", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "BRITE")
SAMPLE = 2_000_000


def label_values(fields, cols):
    """{label: set of its values} for one eggNOG line (empty sets for '-')"""
    out = {}
    for name in LABELS:
        i = cols.get(name)
        v = fields[i] if i is not None and i < len(fields) else "-"
        if v in ("-", ""):
            out[name] = set()
        elif name == "COG_category":
            out[name] = set(v)  # 'EG' is E and G
        elif name == "KEGG_Pathway":
            out[name] = {p for p in v.split(",") if p.startswith("map")}
        else:
            out[name] = set(v.split(","))
    return out


def main(eggnog_db, nb_db):
    t0 = time.time()
    db = sqlite3.connect(f"file:{eggnog_db}?mode=ro", uri=True)
    header = db.execute("SELECT annotation_header FROM eggnog_metadata").fetchone()[0].lstrip("#").split("\t")
    cols = {c: header.index(c) for c in LABELS if c in header}
    missing = [c for c in LABELS if c not in cols]
    if missing:
        print(f"(not in the eggNOG files: {', '.join(missing)})")

    def bits(line):
        """which labels a gene has, as bits in the order of LABELS"""
        if line is None:
            return 0
        v = label_values(line.split("\t"), cols)
        return sum(1 << k for k, name in enumerate(LABELS) if v[name])

    db.create_function("label_bits", 1, bits, deterministic=True)

    def show(title, counts, total):
        print(f"\n{title}: {total:,}")
        for k, name in enumerate(LABELS):
            n = sum(c for b, c in counts.items() if b >> k & 1)
            print(f"  {name:13s} {n:>12,}  {100 * n / max(total, 1):5.1f}%")

    # 1. genes with eggNOG (a sample)
    print("1. reading eggNOG genes...", flush=True)
    counts = Counter(dict(db.execute("SELECT label_bits(annotation_line), COUNT(*) FROM "
                                     "(SELECT annotation_line FROM eggnog LIMIT ?) GROUP BY 1", (SAMPLE,))))
    show(f"1. Genes with eggNOG (first {SAMPLE:,})", counts, sum(counts.values()))
    print(f"  ({time.time() - t0:.0f}s)", flush=True)

    # 2. neighbour and target genes of the smORFs
    db.execute("ATTACH DATABASE ? AS nb", (f"file:{nb_db}?mode=ro",))
    print("\n2. reading the neighbour and target genes of the smORFs...", flush=True)
    rows = db.execute(
        "SELECT n.is_target, label_bits(e.annotation_line), COUNT(*) FROM nb.neighbours n "
        "JOIN nb.occurrences o ON o.occ_id = n.occ_id "
        "LEFT JOIN eggnog e ON e.sample = o.sample AND e.query_id = n.original_id GROUP BY 1, 2").fetchall()
    for target, title in ((0, "2a. Neighbour genes (rows)"), (1, "2b. Target genes (rows)")):
        c = Counter({b: n for t, b, n in rows if t == target})
        show(title, c, sum(c.values()))
    print(f"  ({time.time() - t0:.0f}s)", flush=True)

    # 3. smORFs: labelled if the target gene of at least one occurrence has the label
    print("\n3. reading the target genes of the smORFs...", flush=True)
    n_smorfs = db.execute("SELECT COUNT(DISTINCT smorf_id) FROM nb.occurrences").fetchone()[0]
    has = defaultdict(set)  # label -> smORFs with it
    values = defaultdict(Counter)  # label -> value -> smORFs
    per_smorf = defaultdict(lambda: defaultdict(set))
    for smorf, line in db.execute(
            "SELECT o.smorf_id, e.annotation_line FROM nb.neighbours n "
            "JOIN nb.occurrences o ON o.occ_id = n.occ_id "
            "JOIN eggnog e ON e.sample = o.sample AND e.query_id = n.original_id WHERE n.is_target = 1"):
        for name, v in label_values(line.split("\t"), cols).items():
            if v:
                per_smorf[smorf][name] |= v
    for smorf, labels in per_smorf.items():
        for name, v in labels.items():
            has[name].add(smorf)
            values[name].update(v)
    print(f"\n3. smORFs with a label (target gene of at least one occurrence has it), of {n_smorfs:,}")
    print(f"  {'label':13s} {'smORFs':>9s}  {'%':>5s}  {'different labels':>16s}  most common")
    for name in LABELS:
        top = ", ".join(f"{v} ({n:,})" for v, n in values[name].most_common(3))
        print(f"  {name:13s} {len(has[name]):>9,}  {100 * len(has[name]) / max(n_smorfs, 1):5.1f}  "
              f"{len(values[name]):>16,}  {top}")
    both = has["COG_category"] & has["KEGG_Module"]
    print(f"\n  smORFs with both a COG category and a KEGG module: {len(both):,}")
    print(f"  smORFs with a KEGG module but no COG category: {len(has['KEGG_Module'] - has['COG_category']):,}")
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    if len(sys.argv) > 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1] if len(sys.argv) > 1 else f"{NB_DIR}/urban_soil.db",
         sys.argv[2] if len(sys.argv) > 2 else f"{NB_DIR}/smorf_neighbourhoods.sqlite")
