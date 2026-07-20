#!/usr/bin/env python3

import gzip
import glob
import os
import re
from collections import Counter
import pandas as pd

SMORF_FILE = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"

PRODIGAL_DIR = "/work/microbiome/shanghai_dogs/intermediate-outputs/Prodigal"

EGGNOG_DIR = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"


print("Loading smORF catalog...")

smorf_dict = {}

with gzip.open(SMORF_FILE, "rt") as f:

    header = next(f)

    for line in f:

        if line.strip() == "":
            continue

        cols = line.strip().split("\t")

        smorf_id = cols[0]
        sample = cols[1]
        contig = cols[2]

        coords = cols[3]
        strand = cols[4]

        start, end = coords.split("-")

        start = int(start)
        end = int(end)

        key = (sample, contig, start, end, strand)

        smorf_dict[key] = smorf_id

print("smORFs loaded:", len(smorf_dict))

# STEP 2


print("\nReading Prodigal proteins...")

query_to_smorf = {}

faa_files = sorted(glob.glob(PRODIGAL_DIR + "/D*/D*_proteins.faa.gz"))

header_regex = re.compile(
    r'^>(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)'
)

for faa in faa_files:

    sample = os.path.basename(faa).split("_")[0]

    print(sample)

    with gzip.open(faa, "rt") as f:

        for line in f:

            if not line.startswith(">"):
                continue

            m = header_regex.search(line)

            if m is None:
                continue

            query = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))

            strand = "+" if m.group(4) == "1" else "-"

            contig = query.rsplit("_", 1)[0]

            key = (sample, contig, start, end, strand)

            if key in smorf_dict:

                smorf = smorf_dict[key]

                query_to_smorf[(sample, query)] = smorf

print("\nMapped queries:", len(query_to_smorf))

# STEP 3
# Read EggNOG hits

print("\nReading EggNOG hits...")

results = []

hit_files = sorted(glob.glob(EGGNOG_DIR + "/D*/D*.emapper.hits"))

for hitfile in hit_files:

    sample = os.path.basename(hitfile).split(".")[0]

    print(sample)

    with open(hitfile) as f:

        for line in f:

            if line.startswith("#"):
                continue

            cols = line.strip().split()

            query = cols[0]

            key = (sample, query)

            if key not in query_to_smorf:
                continue

            smorf = query_to_smorf[key]

            subject = cols[1]

            identity = float(cols[2])

            aln_length = int(cols[3])

            evalue = float(cols[10].replace("E","e"))

            bitscore = float(cols[11])

            query_cov = float(cols[12])

            subject_cov = float(cols[13])

            results.append({
                "Sample": sample,
                "SmORF": smorf,
                "Query": query,
                "Seed_Ortholog": subject,
                "Identity": identity,
                "Alignment_length": aln_length,
                "Query_coverage": query_cov,
                "Subject_coverage": subject_cov,
                "Evalue": evalue,
                "Bitscore": bitscore
            })

# STEP 4
# Save table

df = pd.DataFrame(results)

print("\nMatched smORFs:", len(df))

df.to_csv(
    "smorf_subject_coverage.tsv",
    sep="\t",
    index=False
)

# STEP 5
# Coverage statistics

cov = df["Subject_coverage"]

print("\nCoverage summary")

print(cov.describe())


# STEP 6
# Coverage bins


bins = Counter()

for x in cov:

    if x < 20:
        bins["<20"] += 1

    elif x < 30:
        bins["20-30"] += 1

    elif x < 50:
        bins["30-50"] += 1

    elif x < 70:
        bins["50-70"] += 1

    elif x < 90:
        bins["70-90"] += 1

    else:
        bins[">=90"] += 1

print("\nCoverage bins")

total = len(cov)

for b in ["<20","20-30","30-50","50-70","70-90",">=90"]:

    n = bins[b]

    pct = (n/total)*100 if total else 0

    print(f"{b:>6} : {n:8d} ({pct:6.2f}%)")

# STEP 7
# Thresholds

print("\nThreshold statistics")

for t in [20,30,40,50,60,70,80,90]:

    n = (cov < t).sum()

    pct = n/len(cov)*100 if len(cov) else 0

    print(f"< {t:2d}% : {n:8d} ({pct:6.2f}%)")

# STEP 8
# Save threshold summary


summary = []

for t in [20,30,40,50,60,70,80,90]:

    n = (cov < t).sum()

    pct = n/len(cov)*100 if len(cov) else 0

    summary.append([t,n,pct])

summary = pd.DataFrame(
    summary,
    columns=["Threshold","Count","Percent"]
)

summary.to_csv(
    "coverage_summary.tsv",
    sep="\t",
    index=False
)

print("\nDone.")
