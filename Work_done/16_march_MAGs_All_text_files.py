import pandas as pd
from pathlib import Path

BASE = Path("/work/microbiome/shanghai_dogs/resource_generation/MAGs_Onehealth/External_cohorts")
METADATA = BASE / "Skani_lists" / "spire_v1_genome_metadata.tsv"
QUALITY_DIR = BASE / "Skani_lists" / "Quality_MAGs"

# Load metadata
meta = pd.read_csv(METADATA, sep="\t", low_memory=False)

# Assign quality
def assign_quality(row):
    comp = float(row["completeness"])
    cont = float(row["contamination"])

    if comp >= 90 and cont < 5:
        return "HQ"
    elif comp >= 50 and cont < 10:
        return "MQ"
    else:
        return "LQ"

meta["quality"] = meta.apply(assign_quality, axis=1)

hq = meta[meta["quality"] == "HQ"]
mq = meta[meta["quality"] == "MQ"]

cohort_dirs = {
    "Coelho_2018_dog": BASE / "Coelho_2018_dog",
    "Wang_2019_dogs": BASE / "Wang_2019_dogs",
    "Yarlagadda_2022_global_dog": BASE / "Yarlagadda_2022_global_dog",
    "Allaway_2020_dogs": BASE / "Allaway_2020_dogs",
    "Liu_2021_Canidae": BASE / "Liu_2021_Canidae",
    "Xu_2019_dogs": BASE / "Xu_2019_dogs",
    "Worsley-Tonks_2020_dog": BASE / "Worsley-Tonks_2020_dog",
}

valid_suffixes = [".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"]

for cohort_name, cohort_root in cohort_dirs.items():
    mags_dir = cohort_root / "mags"

    files = []
    for pattern in ["*.fa", "*.fna", "*.fasta", "*.fa.gz", "*.fna.gz", "*.fasta.gz"]:
        files.extend(mags_dir.rglob(pattern))

    genome_to_path = {}
    for f in files:
        stem = f.name
        for suf in valid_suffixes:
            if stem.endswith(suf):
                stem = stem[: -len(suf)]
                break
        genome_to_path[stem] = str(f.resolve())

    mag_ids = list(genome_to_path.keys())

    cohort_hq = hq[hq["genome_id"].isin(mag_ids)]
    cohort_mq = mq[mq["genome_id"].isin(mag_ids)]

    # ALL = HQ + MQ only
    all_ids = list(dict.fromkeys(
        list(cohort_hq["genome_id"]) + list(cohort_mq["genome_id"])
    ))

    all_file = QUALITY_DIR / f"{cohort_name}_ALL_list.txt"

    with open(all_file, "w") as f:
        for gid in all_ids:
            f.write(genome_to_path[gid] + "\n")

    print(f"{cohort_name}_ALL_list.txt updated: {len(all_ids)} genomes")
