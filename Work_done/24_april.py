#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

BASE = Path("/work/microbiome/shanghai_dogs/resource_generation/MAGs_Onehealth/External_cohorts/Skani_lists")

SKANI_DIR = BASE / "Quality_MAGs" / "Skani_Quality_Results"
META_FILE = BASE / "spire_v1_genome_metadata.tsv"
OUT_FILE = BASE / "SHD_external_MAG_matches.tsv"

print("Loading metadata...")
meta = pd.read_csv(META_FILE, sep="\t", low_memory=False)
meta = meta[["genome_id", "classification"]].rename(columns={"genome_id": "spire_id"})

def clean_id(p):
    name = Path(str(p)).name
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", name)

dfs = []

print("Processing files...")
for f in sorted(SKANI_DIR.glob("*_ani.tsv")):
    name = f.name.replace("_ani.tsv", "")
    parts = name.split("_")

    quality = parts[-1]
    if quality == "ALL":
        continue

    df = pd.read_csv(f, sep="\t")

    df["spire_id"] = df["Ref_file"].apply(clean_id)
    df["shd_id"] = df["Query_file"].apply(clean_id)
    df["cohort"] = "_".join(parts[:-1])
    df["quality"] = quality

    df = df.rename(columns={
        "Align_fraction_ref": "cov_ref",
        "Align_fraction_query": "cov_query"
    })

    df = df[[
        "spire_id",
        "shd_id",
        "cohort",
        "quality",
        "ANI",
        "cov_ref",
        "cov_query"
    ]]

    dfs.append(df)

print("Combining...")
final_df = pd.concat(dfs, ignore_index=True)

print("Merging taxonomy...")
final_df = final_df.merge(meta, on="spire_id", how="left")

print("Saving...")
final_df.to_csv(OUT_FILE, sep="\t", index=False)

print(f"Saved: {OUT_FILE}")
