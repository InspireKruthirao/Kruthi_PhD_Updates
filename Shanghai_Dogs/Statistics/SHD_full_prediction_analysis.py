#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import re


ROOT = Path("/work/microbiome/users/kruthi")
FOLDER_LIST = ROOT / "smorf_gff_folders.txt"

all_rows = []


# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def accuracy(df):
    if len(df) == 0:
        return 0.0

    return 100 * df["correct_bool"].mean()


def get_true_cogs(value):
    """
    Convert true_cogs into a set of COG letters.

    Examples:
    J       -> {'J'}
    JK      -> {'J', 'K'}
    ['J']   -> {'J'}
    None    -> empty set
    """

    if pd.isna(value):
        return set()

    text = str(value)

    if text in ["None", "nan", ""]:
        return set()

    return {
        x for x in text
        if x.isupper() and x.isalpha()
    }


def get_top3(value):
    """
    Example:
    J:60.0%, K:25.0%, E:15.0%
    ->
    {'J', 'K', 'E'}
    """

    if pd.isna(value):
        return set()

    return set(
        re.findall(
            r"([A-Z]):",
            str(value)
        )
    )


def get_prediction_set(value):
    """
    Example:
    JKE
    ->
    {'J', 'K', 'E'}
    """

    if pd.isna(value):
        return set()

    return {
        x for x in str(value)
        if x.isupper() and x.isalpha()
    }


# --------------------------------------------------
# Load all 67 bins
# --------------------------------------------------

with open(FOLDER_LIST) as f:
    folders = [
        x.strip()
        for x in f
        if x.strip()
    ]


for folder in folders:

    infile = (
        ROOT
        / folder
        / "smorf_onehot_predictions_validation.tsv"
    )

    if not infile.exists():
        print(f"MISSING: {folder}")
        continue

    df = pd.read_csv(
        infile,
        sep="\t"
    )

    df["source_bin"] = folder

    df["correct_bool"] = (
        df["correct"]
        .astype(str)
        .str.lower()
        .eq("true")
    )

    all_rows.append(df)


data_all = pd.concat(
    all_rows,
    ignore_index=True
)


# --------------------------------------------------
# Keep only smORFs with known true COG
# --------------------------------------------------

validated = data_all[
    data_all["true_cogs"].notna()
    &
    ~data_all["true_cogs"]
    .astype(str)
    .isin(["None", "nan", ""])
].copy()


print()
print("=" * 70)
print("SHANGHAI DOGS — FULL NEIGHBOURHOOD PREDICTION ANALYSIS")
print("=" * 70)

print()
print("DATASET")
print("-" * 70)

print(f"Bins analysed             : {len(folders):,}")
print(f"All predicted smORFs      : {len(data_all):,}")
print(f"Validated annotated smORFs: {len(validated):,}")
print(
    f"Unknown / not validated   : "
    f"{len(data_all) - len(validated):,}"
)


# ==================================================
# ANALYSIS 1
# OVERALL ACCURACY
# ==================================================

print()
print("=" * 70)
print("1. OVERALL TOP-1 ACCURACY")
print("=" * 70)

correct = int(
    validated["correct_bool"].sum()
)

print(
    f"Validated : {len(validated):,}"
)

print(
    f"Correct   : {correct:,}"
)

print(
    f"Accuracy  : {accuracy(validated):.3f}%"
)


# ==================================================
# ANALYSIS 2
# ACCURACY BY OCCURRENCE NUMBER
# ==================================================

print()
print("=" * 70)
print("2. ACCURACY BY NUMBER OF smORF OCCURRENCES")
print("=" * 70)


occurrence_bins = [
    0,
    5,
    10,
    20,
    50,
    100,
    np.inf
]

occurrence_labels = [
    "1-5",
    "6-10",
    "11-20",
    "21-50",
    "51-100",
    ">100"
]


validated["occurrence_group"] = pd.cut(
    validated["n_occurrences"],
    bins=occurrence_bins,
    labels=occurrence_labels,
    right=True
)


for group in occurrence_labels:

    sub = validated[
        validated["occurrence_group"] == group
    ]

    if len(sub) == 0:
        continue

    print(
        f"{group:8s} | "
        f"N = {len(sub):7,} | "
        f"Correct = {sub['correct_bool'].sum():7,} | "
        f"Accuracy = {accuracy(sub):6.2f}%"
    )


# ==================================================
# ANALYSIS 3
# VOTE SCORE VS ACCURACY
# ==================================================

print()
print("=" * 70)
print("3. VOTE SCORE VS PREDICTION ACCURACY")
print("=" * 70)


vote_bins = [
    0,
    30,
    40,
    50,
    60,
    70,
    80,
    90,
    100.01
]

vote_labels = [
    "0-<30",
    "30-<40",
    "40-<50",
    "50-<60",
    "60-<70",
    "70-<80",
    "80-<90",
    "90-100"
]


validated["vote_group"] = pd.cut(
    validated["vote_score_pct"],
    bins=vote_bins,
    labels=vote_labels,
    right=False,
    include_lowest=True
)


for group in vote_labels:

    sub = validated[
        validated["vote_group"] == group
    ]

    if len(sub) == 0:
        continue

    print(
        f"{group:8s} | "
        f"N = {len(sub):7,} | "
        f"Correct = {sub['correct_bool'].sum():7,} | "
        f"Accuracy = {accuracy(sub):6.2f}%"
    )


# ==================================================
# ANALYSIS 4
# HIGH VOTE SCORE THRESHOLDS
# ==================================================

print()
print("=" * 70)
print("4. ACCURACY ABOVE DIFFERENT VOTE-SCORE THRESHOLDS")
print("=" * 70)


for threshold in [
    40,
    50,
    60,
    70,
    80,
    90
]:

    sub = validated[
        validated["vote_score_pct"]
        >= threshold
    ]

    print(
        f"Vote >= {threshold:2d}% | "
        f"N = {len(sub):7,} | "
        f"Accuracy = {accuracy(sub):6.2f}%"
    )


# ==================================================
# ANALYSIS 5
# CONSERVATION VS ACCURACY
# ==================================================

print()
print("=" * 70)
print("5. NEIGHBOURHOOD CONSERVATION VS ACCURACY")
print("=" * 70)


conservation_bins = [
    0,
    50,
    75,
    90,
    100.01
]

conservation_labels = [
    "0-<50",
    "50-<75",
    "75-<90",
    "90-100"
]


validated["conservation_group"] = pd.cut(
    validated["conservation_pct"],
    bins=conservation_bins,
    labels=conservation_labels,
    right=False,
    include_lowest=True
)


for group in conservation_labels:

    sub = validated[
        validated["conservation_group"]
        == group
    ]

    if len(sub) == 0:
        continue

    print(
        f"{group:8s} | "
        f"N = {len(sub):7,} | "
        f"Correct = {sub['correct_bool'].sum():7,} | "
        f"Accuracy = {accuracy(sub):6.2f}%"
    )


# ==================================================
# ANALYSIS 6
# CONFIDENCE TIER
# ==================================================

print()
print("=" * 70)
print("6. ACCURACY BY CURRENT CONFIDENCE TIER")
print("=" * 70)


for tier in [
    "VERY HIGH",
    "HIGH",
    "MEDIUM",
    "LOW"
]:

    sub = validated[
        validated["confidence_tier"] == tier
    ]

    if len(sub) == 0:
        continue

    print(
        f"{tier:10s} | "
        f"N = {len(sub):7,} | "
        f"Correct = {sub['correct_bool'].sum():7,} | "
        f"Accuracy = {accuracy(sub):6.2f}%"
    )


# ==================================================
# ANALYSIS 7
# TOP-3 RECOVERY
# ==================================================

print()
print("=" * 70)
print("7. TOP-3 PREDICTION RECOVERY")
print("=" * 70)


validated["true_set"] = (
    validated["true_cogs"]
    .apply(get_true_cogs)
)

validated["top3_set"] = (
    validated["top3_candidates"]
    .apply(get_top3)
)


validated["top3_correct"] = validated.apply(
    lambda row:
    len(
        row["true_set"]
        &
        row["top3_set"]
    ) > 0,
    axis=1
)


top3_correct = int(
    validated["top3_correct"].sum()
)

top3_accuracy = (
    100
    * top3_correct
    / len(validated)
)


print(
    f"Top-1 correct : "
    f"{validated['correct_bool'].sum():,} "
    f"({accuracy(validated):.2f}%)"
)

print(
    f"Top-3 correct : "
    f"{top3_correct:,} "
    f"({top3_accuracy:.2f}%)"
)


wrong_top1 = validated[
    ~validated["correct_bool"]
]

rescued_top3 = wrong_top1[
    wrong_top1["top3_correct"]
]


print()
print(
    f"Top-1 incorrect predictions : "
    f"{len(wrong_top1):,}"
)

print(
    f"Recovered within top-3      : "
    f"{len(rescued_top3):,}"
)

if len(wrong_top1) > 0:

    print(
        f"Top-3 rescue rate           : "
        f"{100 * len(rescued_top3) / len(wrong_top1):.2f}%"
    )


# ==================================================
# ANALYSIS 8
# 90% PREDICTION SET
# ==================================================

print()
print("=" * 70)
print("8. TRUE COG PRESENT IN 90% PREDICTION SET")
print("=" * 70)


validated["prediction90_set"] = (
    validated["prediction_set_90pct"]
    .apply(get_prediction_set)
)


validated["prediction90_correct"] = (
    validated.apply(
        lambda row:
        len(
            row["true_set"]
            &
            row["prediction90_set"]
        ) > 0,
        axis=1
    )
)


prediction90_correct = int(
    validated[
        "prediction90_correct"
    ].sum()
)


prediction90_accuracy = (
    100
    * prediction90_correct
    / len(validated)
)


print(
    f"True COG included: "
    f"{prediction90_correct:,} / "
    f"{len(validated):,}"
)

print(
    f"Coverage         : "
    f"{prediction90_accuracy:.2f}%"
)


# ==================================================
# ANALYSIS 9
# VERY HIGH CONFIDENCE DETAILS
# ==================================================

print()
print("=" * 70)
print("9. VERY HIGH CONFIDENCE PREDICTIONS")
print("=" * 70)


vh = validated[
    validated["confidence_tier"]
    == "VERY HIGH"
]


print(
    f"N        : {len(vh):,}"
)

print(
    f"Correct  : "
    f"{vh['correct_bool'].sum():,}"
)

print(
    f"Accuracy : "
    f"{accuracy(vh):.2f}%"
)


if len(vh) > 0:

    print(
        f"Mean vote score    : "
        f"{vh['vote_score_pct'].mean():.2f}%"
    )

    print(
        f"Mean conservation  : "
        f"{vh['conservation_pct'].mean():.2f}%"
    )

    print(
        f"Median occurrences : "
        f"{vh['n_occurrences'].median():.1f}"
    )


# ==================================================
# FINAL SUMMARY
# ==================================================

print()
print("=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(
    f"Overall top-1 accuracy       : "
    f"{accuracy(validated):.2f}%"
)

print(
    f"VERY HIGH accuracy           : "
    f"{accuracy(vh):.2f}%"
)

print(
    f"Top-3 accuracy               : "
    f"{top3_accuracy:.2f}%"
)

print(
    f"90% prediction-set coverage  : "
    f"{prediction90_accuracy:.2f}%"
)

print()
print("=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
