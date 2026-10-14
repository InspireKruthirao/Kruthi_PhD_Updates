#!/usr/bin/env python3
"""
How good are the eggNOG hits behind the COG categories of the neighbour and
target genes in smorf_neighbourhoods.sqlite?

For every gene with an eggNOG annotation, takes its seed ortholog hit (the
hit eggNOG took the function from) from <sample>.emapper.seed_orthologs.gz:
pident (% identity), qcov (% of our protein aligned), scov (% of the eggNOG
protein aligned), evalue and bitscore. Then reports:

  - how the neighbour genes and the target genes are spread over ranges of
    qcov, scov and pident
  - for candidate thresholds: the % of neighbour genes with a COG category
    that pass, and how many smORFs stay annotated (the target gene of at
    least one occurrence has a COG category and passes)

Writes eggnog_hit_quality.sqlite (one row per gene: sample, original_id and
the hit values) for filtering later, and eggnog_hit_quality_summary.tsv.

    python3 eggnog_hit_quality.py
"""

import gzip
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

EGGNOG_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/egg_nog/results")
OUT_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL")
NB_DB = OUT_DIR / "smorf_neighbourhoods.sqlite"
HITS_DB = OUT_DIR / "eggnog_hit_quality.sqlite"
SUMMARY = OUT_DIR / "eggnog_hit_quality_summary.tsv"

RANGES = (0, 30, 50, 70, 90, 100.01)  # % bins for qcov, scov, pident
THRESHOLDS = [  # (min qcov, min scov, min pident, max evalue)
    (0, 0, 0, 1),
    (50, 0, 0, 1),
    (70, 0, 0, 1),
    (90, 0, 0, 1),
    (50, 50, 0, 1),
    (70, 70, 0, 1),
    (0, 0, 30, 1),
    (0, 0, 50, 1),
    (50, 50, 30, 1),
    (70, 70, 30, 1e-10),
    (90, 90, 50, 1e-10),
]


def genes_needed():
    """{sample: {original_id: [neighbour rows, target rows, with category]}}
    for the genes with eggNOG, and {(sample, original_id): smORFs} for the
    target genes with a COG category"""
    db = sqlite3.connect(f"file:{NB_DB}?mode=ro", uri=True)
    need = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    target_smorfs = defaultdict(set)
    rows = db.execute(
        "SELECT o.sample, o.smorf_id, n.original_id, n.is_target, n.cog_category != '' "
        "FROM neighbours n JOIN occurrences o ON o.occ_id = n.occ_id WHERE n.has_eggnog = 1")
    for number, (sample, smorf, original, is_target, has_cat) in enumerate(rows, 1):
        counts = need[sample][original]
        counts[1 if is_target else 0] += 1
        counts[2] += has_cat
        if is_target and has_cat:
            target_smorfs[(sample, original)].add(smorf)
        if number % 2_000_000 == 0:
            print(f"  {number:,} genes read", flush=True)
    db.close()
    return need, target_smorfs


def read_hits(need):
    """{(sample, original_id): (evalue, bitscore, pident, qcov, scov)}"""
    hits = {}
    for number, (sample, ids) in enumerate(sorted(need.items()), 1):
        path = EGGNOG_DIR / sample / f"{sample}.emapper.seed_orthologs.gz"
        if not path.is_file():
            print(f"  missing {path}", flush=True)
            continue
        found = 0
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] in ids:
                    hits[(sample, f[0])] = (float(f[2]), float(f[3]), float(f[8]), float(f[9]), float(f[10]))
                    found += 1
        print(f"  [{number}/{len(need)}] {sample}: {found:,} of {len(ids):,} genes found", flush=True)
    return hits


def passes(hit, threshold):
    qcov, scov, pident, evalue = threshold
    return hit[3] >= qcov and hit[4] >= scov and hit[2] >= pident and hit[0] <= evalue


def main():
    t0 = time.time()
    print(f"Reading genes from {NB_DB}", flush=True)
    need, target_smorfs = genes_needed()
    print(f"Reading seed orthologs from {EGGNOG_DIR}", flush=True)
    hits = read_hits(need)

    db = sqlite3.connect(str(HITS_DB))
    db.execute("DROP TABLE IF EXISTS hits")
    db.execute("CREATE TABLE hits (sample TEXT, original_id TEXT, evalue REAL, bitscore REAL, "
               "pident REAL, qcov REAL, scov REAL, PRIMARY KEY (sample, original_id)) WITHOUT ROWID")
    db.executemany("INSERT INTO hits VALUES (?, ?, ?, ?, ?, ?, ?)", ((*k, *v) for k, v in hits.items()))
    db.commit()
    db.close()

    out = []
    def say(text=""):
        print(text, flush=True)
        out.append(text)

    n_nb = sum(c[0] for ids in need.values() for c in ids.values())
    n_tg = sum(c[1] for ids in need.values() for c in ids.values())
    missing = sum(c[0] + c[1] for s, ids in need.items() for g, c in ids.items() if (s, g) not in hits)
    say(f"genes with eggNOG: {n_nb:,} neighbour rows, {n_tg:,} target rows; "
        f"without a seed ortholog line: {missing:,}")

    for col, name in ((3, "qcov"), (4, "scov"), (2, "pident")):
        say(f"\n{name} (% )\tneighbour genes\t%\ttarget genes\t%")
        for lo, hi in zip(RANGES, RANGES[1:]):
            nb = tg = 0
            for s, ids in need.items():
                for g, c in ids.items():
                    h = hits.get((s, g))
                    if h and lo <= h[col] < hi:
                        nb += c[0]
                        tg += c[1]
            label = f"{lo}-{min(hi, 100):.0f}"
            say(f"{label}\t{nb:,}\t{100 * nb / max(n_nb, 1):.1f}\t{tg:,}\t{100 * tg / max(n_tg, 1):.1f}")

    all_smorfs = set().union(*target_smorfs.values()) if target_smorfs else set()
    say(f"\nthreshold (qcov, scov, pident, evalue)\tneighbour genes with category kept %\t"
        f"annotated smORFs kept\t%")
    n_cat = sum(c[2] for ids in need.values() for c in ids.values())
    for th in THRESHOLDS:
        kept = sum(c[2] for s, ids in need.items() for g, c in ids.items()
                   if (s, g) in hits and passes(hits[(s, g)], th))
        smorfs = set()
        for key, sm in target_smorfs.items():
            if key in hits and passes(hits[key], th):
                smorfs |= sm
        say(f"qcov>={th[0]} scov>={th[1]} pident>={th[2]} evalue<={th[3]:g}\t"
            f"{100 * kept / max(n_cat, 1):.1f}\t{len(smorfs):,}\t{100 * len(smorfs) / max(len(all_smorfs), 1):.1f}")

    SUMMARY.write_text("\n".join(out) + "\n")
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Saved {SUMMARY} and {HITS_DB}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(__doc__)
        sys.exit(1)
    main()
