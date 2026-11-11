#!/usr/bin/env python3
"""Map every ORF to its 100NT representative and write <prefix>.100NT.matches.xz.

Columns: ORF, '=' (is the representative) or 'C' (collapsed), 100NT ID.
"""
import argparse

import numpy as np
import polars as pl

from common import MISSING, check, flag, pad9, read_tsv, to_index, write_compressed

ap = argparse.ArgumentParser()
ap.add_argument("--prefix", required=True)
ap.add_argument("--rep-ids", required=True, help="ORF IDs of the 100NT representatives, in 100NT order")
ap.add_argument("--clusters", required=True, help="MMseqs2 cluster.tsv (rep, member)")
ap.add_argument("--orig", required=True, help="ORF origin parquet (for the ORF count)")
ap.add_argument("--matches", required=True)
ap.add_argument("--table", required=True, help="output parquet: orf, is_rep, i100")
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
P = a.prefix

n_orfs = pl.scan_parquet(a.orig).select(pl.len()).collect().item()

# ORF index -> 100NT index
reps = to_index(read_tsv(a.rep_ids, ["id"])["id"], f"{P}.ORF").to_numpy()
check(len(np.unique(reps)) == len(reps), "duplicate 100NT representatives")
check(reps.max() < n_orfs, "representative ORF outside the ORF catalogue")
lut = np.full(n_orfs, MISSING, dtype=np.uint32)
lut[reps] = np.arange(len(reps), dtype=np.uint32)

cl = read_tsv(a.clusters, ["rep", "member"])
rep = to_index(cl["rep"], f"{P}.ORF")
orf = to_index(cl["member"], f"{P}.ORF")
del cl
tab = pl.DataFrame({"orf": orf, "is_rep": orf == rep, "i100": lut[rep.to_numpy()]})

check(len(tab) == n_orfs, f"{len(tab):,} cluster rows but {n_orfs:,} ORFs")
check(tab["orf"].n_unique() == n_orfs, "some ORFs are missing or listed twice")
check((tab["i100"] != MISSING).all(), "cluster representative not in rep FASTA")
check(tab["is_rep"].sum() == len(reps), "representative rows do not match rep FASTA")

tab.write_parquet(a.table)
write_compressed(tab.select(pad9(f"{P}.ORF", pl.col("orf")), flag("is_rep", "=", "C"),
                            pad9(f"{P}.100NT", pl.col("i100"))), a.matches, a.threads)
print(f"100NT: {len(reps):,} representatives for {n_orfs:,} ORFs", flush=True)
