#!/usr/bin/env python3
"""
Neighbourhoods of normal proteins, in the same format as
smorf_neighbourhoods.sqlite, so the same rules (sweep_db.py) and models can
be trained and tested on them (supervisor meeting 28 October 2026: "benchmark
the rules on all proteins"; and a model trained on proteins needs no
annotated smORFs).

A random sample of genes from urban_soil.db (read only, nothing in it is
changed) is taken; a gene is kept if it is a normal protein (> MIN_AA amino
acids) with an eggNOG COG category. Each kept protein is one entry (one
occurrence) with WINDOW genes on each side, exactly as for a smORF; its own
COG category is the known category (target_categories).

Writes a NEW database (the smORF databases are not touched):
    Neighbourhood_Analysis_SQL/Protein_Neighbourhoods/protein_neighbourhoods.sqlite
with the tables occurrences, neighbours and smorfs (one row per protein,
smorf_id = the protein's ORF ID, bin_name = PROTEINS). The files table is
left empty.

    python3 protein_neighbourhoods_sqlite.py chunk I N   piece I of N (0...N-1) -> parts/chunk_I_of_N.sqlite
    python3 protein_neighbourhoods_sqlite.py merge N     the N pieces -> protein_neighbourhoods.sqlite

run_protein_neighbourhoods.sh submits the N pieces as one PBS array job and
then the merge.
"""

import json
import os
import random
import sqlite3
import sys
import time
from pathlib import Path

import smorf_neighbourhoods_sqlite as nb  # same urban_soil.db, same neighbourhood code

OUT_DIR = Path(os.environ.get(
    "PROTEIN_OUT_DIR",
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL/Protein_Neighbourhoods"))
OUT_DB = OUT_DIR / "protein_neighbourhoods.sqlite"
N_CANDIDATES = 400_000  # genes drawn at random; about half are kept
MIN_AA = 100            # normal proteins only (smORFs are <= 100 aa)
SEED = 0
BIN = "PROTEINS"
PROGRESS_EVERY = 1000


def candidates(connection):
    """N_CANDIDATES random gene rowids, the same in every piece"""
    top = connection.execute("SELECT MAX(rowid) FROM genes").fetchone()[0]
    rng = random.Random(SEED)
    return sorted(rng.sample(range(1, top + 1), min(N_CANDIDATES, top)))


def run_chunk(i, n):
    connection = nb.connect(read_only=True)
    rowids = candidates(connection)[i::n]
    parts = OUT_DIR / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    part = parts / f"chunk_{i}_of_{n}.sqlite"
    tmp = part.with_suffix(".sqlite.part")
    tmp.unlink(missing_ok=True)
    out = sqlite3.connect(str(tmp))
    out.execute("PRAGMA journal_mode=OFF")
    out.execute("PRAGMA synchronous=OFF")
    out.executescript(nb.SCHEMA)

    label = f"chunk_{i}_of_{n}"
    print(f"[{label}] {len(rowids):,} candidate genes", flush=True)
    t0 = time.time()
    count = {"candidates": 0, "short": 0, "no_cog": 0, "kept": 0}
    occ_rows, gene_rows, smorf_rows = [], [], []

    for number, rowid in enumerate(rowids, 1):
        count["candidates"] += 1
        row = connection.execute(
            "SELECT orf, sample, base_contig, annotation_query, gene_start, gene_end FROM genes WHERE rowid = ?",
            (rowid,)).fetchone()
        if row is None:
            continue
        orf, sample, contig, query, start, end = row
        if (end - start + 1) // 3 - 1 <= MIN_AA:  # amino acids (minus the stop codon)
            count["short"] += 1
            continue
        own = connection.execute("SELECT annotation_line FROM eggnog WHERE sample = ? AND query_id = ?",
                                 (sample, query)).fetchone()
        if not nb.cog_parts(nb.cog_name(own[0] if own else None))[1]:
            count["no_cog"] += 1
            continue
        genes = connection.execute(
            "SELECT orf, gene_start, gene_end, strand, partial, annotation_query FROM genes "
            "WHERE sample = ? AND base_contig = ? ORDER BY gene_start", (sample, contig)).fetchall()
        at = next(k for k, g in enumerate(genes) if g[0] == orf)
        left = max(0, at - nb.WINDOW)
        genes, anchor = genes[left:at + nb.WINDOW + 1], at - left
        ids = [g[5] for g in genes]
        annotations = dict(connection.execute(
            f"SELECT query_id, annotation_line FROM eggnog WHERE sample = ? AND query_id IN "
            f"({','.join('?' * len(ids))})", (sample, *ids)))
        occ_id = f"{orf}|{sample}|{contig}|{start}"
        rows, occurrence = nb.occurrence_rows(occ_id, genes, anchor, annotations)
        occ_rows.append((occ_id, orf, BIN, sample, contig, start, end, "protein", len(genes)))
        gene_rows.extend(rows)
        smorf_rows.append((BIN, orf, 1, json.dumps([occurrence])))
        count["kept"] += 1
        if number % PROGRESS_EVERY == 0 or number == len(rowids):
            out.executemany("INSERT OR REPLACE INTO occurrences VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", occ_rows)
            out.executemany("INSERT OR REPLACE INTO neighbours VALUES "
                            "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", gene_rows)
            out.executemany("INSERT OR REPLACE INTO smorfs VALUES (?, ?, ?, ?)", smorf_rows)
            out.commit()
            occ_rows, gene_rows, smorf_rows = [], [], []
            print(f"[{label}] {number:,}/{len(rowids):,} genes, {count['kept']:,} kept "
                  f"({time.time() - t0:.0f}s)", flush=True)

    out.close()
    tmp.replace(part)
    print(f"[{label}] DONE: {count['kept']:,} of {count['candidates']:,} genes kept; "
          f"{count['short']:,} <= {MIN_AA} aa, {count['no_cog']:,} without a COG category", flush=True)


def merge(n):
    paths = [OUT_DIR / "parts" / f"chunk_{i}_of_{n}.sqlite" for i in range(n)]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        sys.exit(f"{len(missing)} piece(s) not finished yet: {missing[:10]}")
    OUT_DB.unlink(missing_ok=True)
    out = sqlite3.connect(str(OUT_DB))
    out.execute("PRAGMA journal_mode=OFF")
    out.executescript(nb.SCHEMA)
    for number, path in enumerate(paths, 1):
        out.execute("ATTACH DATABASE ? AS part", (str(path),))
        for table in ("occurrences", "neighbours", "smorfs"):
            out.execute(f"INSERT OR REPLACE INTO {table} SELECT * FROM part.{table}")
        out.commit()
        out.execute("DETACH DATABASE part")
    out.execute("CREATE INDEX idx_occ_smorf ON occurrences(smorf_id)")
    out.commit()
    n_prot, n_genes, n_egg = out.execute(
        "SELECT (SELECT COUNT(*) FROM smorfs), COUNT(*), SUM(has_eggnog) FROM neighbours").fetchone()
    out.close()
    print(f"Done: {OUT_DB} ({OUT_DB.stat().st_size / 1e6:.0f} MB): {n_prot:,} proteins, "
          f"{n_genes:,} genes in their neighbourhoods, {n_egg:,} with eggNOG")


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "chunk":
        run_chunk(int(sys.argv[2]), int(sys.argv[3]))
    elif len(sys.argv) == 3 and sys.argv[1] == "merge":
        merge(int(sys.argv[2]))
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
