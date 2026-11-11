"""Shared helpers for the gene catalogue scripts."""
import os
import subprocess

import polars as pl

MISSING = 0xFFFFFFFF  # "no value" marker in uint32 lookup arrays


def check(ok, msg):
    """Stop the job with a clear message if a sanity check fails."""
    if not ok:
        raise SystemExit(f"CHECK FAILED: {msg}")


def pad9(prefix, idx):
    """Polars expression: 1234 -> '<prefix>.000_001_234' (same as Shanghai Dogs lib.pad9)."""
    s = idx.cast(pl.String).str.zfill(9)
    return pl.concat_str(
        pl.lit(f"{prefix}."), s.str.slice(0, 3), pl.lit("_"),
        s.str.slice(3, 3), pl.lit("_"), s.str.slice(6, 3),
    ).alias(prefix)


def to_index(ids, prefix):
    """'<prefix>.000_001_234' -> 1234. Fails on any ID with a different prefix."""
    check(ids.str.starts_with(f"{prefix}.").all(), f"IDs not starting with {prefix}.")
    return ids.str.strip_prefix(f"{prefix}.").str.replace_all("_", "").cast(pl.UInt32, strict=True)


def read_tsv(path, names):
    """Read a header-less TSV with all columns as strings."""
    return pl.read_csv(path, separator="\t", has_header=False, new_columns=names,
                       quote_char=None, infer_schema=False)


def write_compressed(df, dst, threads, header=False):
    """Write df as TSV, then compress to dst (.xz with xz, .gz with pigz)."""
    plain = f"{dst}.plain"
    df.write_csv(plain, separator="\t", include_header=header, quote_style="never")
    tool = ["xz", f"-T{threads}"] if dst.endswith(".xz") else ["pigz", "-p", str(threads)]
    with open(dst, "wb") as out:
        subprocess.run([*tool, "-c", plain], stdout=out, check=True)
    os.remove(plain)


def flag(col, yes, no):
    """Polars expression: boolean column -> one of two labels."""
    return pl.when(pl.col(col)).then(pl.lit(yes)).otherwise(pl.lit(no)).alias(col)
