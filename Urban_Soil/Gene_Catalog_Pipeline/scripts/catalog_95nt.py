#!/usr/bin/env python3
"""Map 100NT sequences to 95NT clusters and build the full cluster table.

Writes:
  <prefix>.95NT.matches.tsv.xz  100NT, '=' or 'R' (reduced), 95NT
  <prefix>.clusters.tsv.xz      ORF, '='/'C', 100NT, '='/'R', 95NT  (same order as 100NT matches)
  canonical parquet             the one ORF behind each 95NT representative (for eggNOG transfer)
"""
import argparse

import numpy as np
import polars as pl

from common import MISSING, check, flag, pad9, read_tsv, to_index, write_compressed

ap = argparse.ArgumentParser()
ap.add_argument("--prefix", required=True)
ap.add_argument("--rep-ids", required=True, help="100NT IDs of the 95NT representatives, in 95NT order")
ap.add_argument("--clusters", required=True, help="MMseqs2 cluster.tsv (rep, member)")
ap.add_argument("--table100", required=True, help="parquet from catalog_100nt.py")
ap.add_argument("--matches", required=True)
ap.add_argument("--cluster-table", required=True)
ap.add_argument("--canonical", required=True)
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
P = a.prefix

t100 = pl.read_parquet(a.table100)
n100 = int(t100["is_rep"].sum())

# 100NT representative -> 95NT index
reps = to_index(read_tsv(a.rep_ids, ["id"])["id"], f"{P}.100NT").to_numpy()
check(len(np.unique(reps)) == len(reps), "duplicate 95NT representatives")
check(reps.max() < n100, "95NT representative outside the 100NT catalogue")
rep_to_95 = np.full(n100, MISSING, dtype=np.uint32)
rep_to_95[reps] = np.arange(len(reps), dtype=np.uint32)

cl = read_tsv(a.clusters, ["rep", "member"])
rep = to_index(cl["rep"], f"{P}.100NT").to_numpy()
member = to_index(cl["member"], f"{P}.100NT").to_numpy()
del cl
check(len(member) == n100, f"{len(member):,} cluster rows but {n100:,} 100NT sequences")
check(len(np.unique(member)) == n100, "some 100NT sequences are missing or listed twice")

m95 = pl.DataFrame({"i100": member, "is_rep": member == rep, "i95": rep_to_95[rep]})
check((m95["i95"] != MISSING).all(), "cluster representative not in rep FASTA")
check(m95["is_rep"].sum() == len(reps), "representative rows do not match rep FASTA")
write_compressed(m95.select(pad9(f"{P}.100NT", pl.col("i100")), flag("is_rep", "=", "R"),
                            pad9(f"{P}.95NT", pl.col("i95"))), a.matches, a.threads)

# Per 100NT sequence: its 95NT cluster and whether it is the 95NT representative
i95_of = np.empty(n100, dtype=np.uint32)
i95_of[member] = m95["i95"].to_numpy()
rep95_of = np.zeros(n100, dtype=bool)
rep95_of[member] = m95["is_rep"].to_numpy()
del m95

i100 = t100["i100"].to_numpy()
full = t100.with_columns(i95=pl.Series(i95_of[i100]), is_rep95=pl.Series(rep95_of[i100]))
write_compressed(full.select(pad9(f"{P}.ORF", pl.col("orf")), flag("is_rep", "=", "C"),
                             pad9(f"{P}.100NT", pl.col("i100")), flag("is_rep95", "=", "R"),
                             pad9(f"{P}.95NT", pl.col("i95"))), a.cluster_table, a.threads)

# Exactly one ORF is both 100NT and 95NT representative for each 95NT cluster
canon = full.filter(pl.col("is_rep") & pl.col("is_rep95")).select("orf", "i95")
check(len(canon) == len(reps) and canon["i95"].n_unique() == len(reps), "95NT clusters without exactly one ORF")
canon.write_parquet(a.canonical)
print(f"95NT: {len(reps):,} clusters from {n100:,} 100NT sequences and {len(full):,} ORFs", flush=True)
