#!/usr/bin/env python3
"""
A copy of the smorfs table of smorf_neighbourhoods.sqlite (what sweep_db.py
reads) in which genes whose eggNOG hit is weak count as having no COG
category: as neighbours (they do not vote) and as targets (the smORF loses
that known category).

A gene's hit is its seed ortholog, from eggnog_hit_quality.sqlite (written
by eggnog_hit_quality.py). It passes if qcov >= QCOV, scov >= SCOV,
pident >= PIDENT and evalue <= EVALUE. A gene without a hit does not pass.

    python3 filter_by_hit_quality.py QCOV SCOV PIDENT EVALUE
    python3 filter_by_hit_quality.py 50 50 30 1e-10
        -> smorf_neighbourhoods_q50_s50_p30_e1e-10.sqlite

With 0 0 0 1 every gene passes and the copy equals the original table.
The new file is small (only the smorfs table); the original is not changed.
"""

import json
import sqlite3
import sys
import time
from pathlib import Path

OUT_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL")
NB_DB = OUT_DIR / "smorf_neighbourhoods.sqlite"
HITS_DB = OUT_DIR / "eggnog_hit_quality.sqlite"


def main(qcov, scov, pident, evalue):
    t0 = time.time()
    tag = f"q{qcov:g}_s{scov:g}_p{pident:g}_e{evalue:g}"
    out_path = OUT_DIR / f"smorf_neighbourhoods_{tag}.sqlite"
    keep_all = qcov <= 0 and scov <= 0 and pident <= 0 and evalue >= 1

    good = set()
    if not keep_all:
        hits = sqlite3.connect(f"file:{HITS_DB}?mode=ro", uri=True)
        good = {(s, g) for s, g in hits.execute(
            "SELECT sample, original_id FROM hits WHERE qcov >= ? AND scov >= ? AND pident >= ? "
            "AND evalue <= ?", (qcov, scov, pident, evalue))}
        hits.close()
        print(f"genes passing {tag}: {len(good):,}", flush=True)

    db = sqlite3.connect(f"file:{NB_DB}?mode=ro", uri=True)
    out_path.unlink(missing_ok=True)
    out = sqlite3.connect(str(out_path))
    out.execute("PRAGMA journal_mode=OFF")
    out.execute("CREATE TABLE smorfs (bin_name TEXT NOT NULL, smorf_id TEXT NOT NULL, "
                "n_occurrence_dirs INTEGER NOT NULL, occurrences_json TEXT NOT NULL, "
                "PRIMARY KEY (bin_name, smorf_id))")

    smorfs = db.execute("SELECT bin_name, smorf_id FROM smorfs ORDER BY bin_name, smorf_id").fetchall()
    n_cat = n_dropped = n_annot = 0
    batch = []
    for number, (bin_name, smorf_id) in enumerate(smorfs, 1):
        occs = db.execute("SELECT occ_id, sample FROM occurrences WHERE smorf_id = ? "
                          "ORDER BY sample, contig, smorf_start", (smorf_id,)).fetchall()
        occ_json = []
        for occ_id, sample in occs:
            genes = db.execute("SELECT rank, side, original_id, strand, is_target, dist_bp, gap_bp, cog_id, "
                               "cog_category FROM neighbours WHERE occ_id = ? ORDER BY position",
                               (occ_id,)).fetchall()
            target = next(g for g in genes if g[4])
            cats_of = {}
            for g in genes:
                cats = list(g[8])
                if cats:
                    n_cat += 1
                    if not keep_all and (sample, g[2]) not in good:
                        cats = []
                        n_dropped += 1
                cats_of[g] = cats
            neighbours, info = [], []
            for side in ("up", "down"):
                for g in sorted((g for g in genes if g[1] == side), key=lambda g: g[0]):
                    cats = cats_of[g]
                    neighbours.append([g[0], cats])
                    info.append([side, int(g[3] == target[3]), g[5], g[6], g[7] if cats else ""])
            occ_json.append({"neighbours": neighbours, "neighbour_info": info,
                             "target_categories": cats_of[target]})
        n_annot += any(o["target_categories"] for o in occ_json)
        batch.append((bin_name, smorf_id, len(occ_json), json.dumps(occ_json)))
        if len(batch) == 5000 or number == len(smorfs):
            out.executemany("INSERT INTO smorfs VALUES (?, ?, ?, ?)", batch)
            out.commit()
            batch.clear()
            print(f"  {number:,}/{len(smorfs):,} smORFs ({time.time() - t0:.0f}s)", flush=True)
    out.close()
    print(f"\nDone: {out_path}\n  gene categories dropped: {n_dropped:,} of {n_cat:,} "
          f"({100 * n_dropped / max(n_cat, 1):.1f}%)\n  annotated smORFs left: {n_annot:,} of {len(smorfs):,}",
          flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    main(*map(float, sys.argv[1:]))
