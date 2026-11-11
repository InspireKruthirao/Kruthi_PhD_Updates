#!/usr/bin/env python3
"""Merge per-sample Prodigal ORFs into one catalogue with new IDs (<prefix>.ORF.000_000_000).

Writes: plain ORF FASTA (for MMseqs2), ORF FASTA .xz, origin table .tsv.xz and .parquet.
"""
import argparse
import subprocess

import polars as pl

from common import check, pad9, write_compressed

ap = argparse.ArgumentParser()
ap.add_argument("--prefix", required=True)
ap.add_argument("--samples", nargs="+", required=True)
ap.add_argument("--orfs", nargs="+", required=True, help="per-sample *_orfs.fna.gz, same order as --samples")
ap.add_argument("--fna", required=True)
ap.add_argument("--fna-xz", required=True)
ap.add_argument("--orig-xz", required=True)
ap.add_argument("--orig-parquet", required=True)
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
check(len(a.samples) == len(a.orfs), "number of samples and ORF files differ")

offset = 0
origins = []
with open(a.fna, "wb") as fna:
    for i, (sample, path) in enumerate(zip(a.samples, a.orfs), 1):
        # seqkit gives one row per ORF: full header <tab> sequence <tab>
        tab = subprocess.run(["seqkit", "fx2tab", "-j", "4", path], capture_output=True, check=True).stdout
        df = pl.read_csv(tab, separator="\t", has_header=False, new_columns=["header", "seq", "_"],
                         quote_char=None, infer_schema=False)
        check(len(df) > 0, f"{sample}: no ORFs in {path}")

        # Prodigal header: "<id> # start # end # strand # ID=..;partial=..;..."
        parts = pl.col("header").str.split("#").list
        df = df.with_columns(
            orf=pl.int_range(offset, offset + len(df), dtype=pl.UInt32),
            original_id=parts.get(0).str.strip_chars().str.split(" ").list.first(),
            start=parts.get(1, null_on_oob=True).str.strip_chars().cast(pl.Int64, strict=True),
            end=parts.get(2, null_on_oob=True).str.strip_chars().cast(pl.Int64, strict=True),
            strand=parts.get(3, null_on_oob=True).str.strip_chars().cast(pl.Int8, strict=True),
            partial=pl.col("header").str.extract(r"partial=([0-9]+)").fill_null("NA"),
        )
        check(not df.select(pl.any_horizontal(pl.col("start", "end", "strand", "seq").is_null()).any()).item(),
              f"{sample}: unexpected Prodigal header or empty sequence")
        check(df["original_id"].n_unique() == len(df), f"{sample}: duplicate ORF IDs")

        df.select(pl.concat_str(pl.lit(">"), pad9(f"{a.prefix}.ORF", pl.col("orf")), pl.lit("\n"), "seq")) \
          .write_csv(fna, include_header=False, quote_style="never")
        origins.append(df.select("orf", pl.lit(sample).alias("sample"), "original_id",
                                 "start", "end", "strand", "partial"))
        offset += len(df)
        print(f"{i}/{len(a.samples)} {sample}: {len(df):,} ORFs", flush=True)

print(f"Total ORFs: {offset:,}", flush=True)

# Compress the FASTA (the plain copy is kept for MMseqs2)
with open(a.fna_xz, "wb") as out:
    subprocess.run(["xz", f"-T{a.threads}", "-c", a.fna], stdout=out, check=True)

orig = pl.concat(origins)
orig.write_parquet(a.orig_parquet)

# Origin table in the original column layout
write_compressed(
    orig.select(pad9(f"{a.prefix}.ORF", pl.col("orf")).alias("ORF"), pl.col("sample").alias("Sample"),
                pl.col("original_id").alias("Original_ID"), pl.col("start").alias("Start"),
                pl.col("end").alias("End"), pl.col("strand").alias("Strand"), pl.col("partial").alias("Partial")),
    a.orig_xz, a.threads, header=True)
