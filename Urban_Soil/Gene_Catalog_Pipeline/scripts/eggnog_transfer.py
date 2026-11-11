#!/usr/bin/env python3
"""Give each 95NT cluster the eggNOG annotation of its representative ORF.

Output: <prefix>.95NT.emapper.annotations.gz, first column US.95NT, then the eggNOG columns,
sorted by 95NT ID. Clusters without an eggNOG hit are not listed.
"""
import argparse

import polars as pl

from common import check, pad9, write_compressed

ap = argparse.ArgumentParser()
ap.add_argument("--prefix", required=True)
ap.add_argument("--orig", required=True, help="ORF origin parquet")
ap.add_argument("--canonical", required=True, help="parquet from catalog_95nt.py")
ap.add_argument("--samples", nargs="+", required=True)
ap.add_argument("--annotations", nargs="+", required=True, help="per-sample .emapper.annotations(.gz)")
ap.add_argument("--out", required=True)
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
check(len(a.samples) == len(a.annotations), "number of samples and annotation files differ")

# Representative ORF of each 95NT -> (sample, Prodigal ID)
canon = pl.read_parquet(a.canonical)
wanted = canon.join(pl.read_parquet(a.orig, columns=["orf", "sample", "original_id"]), on="orf")
check(len(wanted) == len(canon), "representative ORFs missing from the origin table")

hits, header = [], None
for sample, path in zip(a.samples, a.annotations):
    ann = pl.read_csv(path, separator="\t", comment_prefix="##", quote_char=None, infer_schema=False)
    header = header or ann.columns
    check(ann.columns == header and header[0] == "#query", f"{sample}: unexpected eggNOG header")
    check(ann["#query"].n_unique() == len(ann), f"{sample}: duplicate eggNOG queries")
    sel = wanted.filter(pl.col("sample") == sample).select(pl.col("original_id").alias("#query"), "i95")
    hit = ann.join(sel, on="#query")
    hits.append(hit)
    print(f"{sample}: {len(hit):,} of {len(sel):,} representatives annotated", flush=True)

out = pl.concat(hits).sort("i95")
check(out["i95"].n_unique() == len(out), "a 95NT cluster got two annotations")
write_compressed(out.select(pad9(f"{a.prefix}.95NT", pl.col("i95")).alias(f"{a.prefix}.95NT"), *header),
                 a.out, a.threads, header=True)
print(f"Annotated 95NT clusters: {len(out):,} of {len(canon):,}", flush=True)
