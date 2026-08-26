cat > summarize_all_shd_bins.py <<'PY'
#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re


ROOT = Path("/work/microbiome/users/kruthi")
FOLDER_LIST = ROOT / "smorf_gff_folders.txt"

OUTPUT = ROOT / "SHD_all_bins_prediction_summary.tsv"


def accuracy(df):
    if len(df) == 0:
        return 0.0

    return 100 * df["correct_bool"].sum() / len(df)


rows = []

with open(FOLDER_LIST) as f:
    folders = [line.strip() for line in f if line.strip()]


print(f"Bins to process: {len(folders)}")
print()


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

    # Convert correct column safely to Boolean
    df["correct_bool"] = (
        df["correct"]
        .astype(str)
        .str.lower()
        .eq("true")
    )

    # Validation only includes smORFs
    # with a known true COG annotation
    validated = df[
        df["true_cogs"].notna()
        &
        ~df["true_cogs"]
        .astype(str)
        .isin(["None", "nan", ""])
    ].copy()

    # Extract occurrence range from folder name
    match = re.search(
        r"SmORF_neighbourhoods_(\d+)_(\d+)$",
        folder
    )

    if match:
        occurrence_min = int(match.group(1))
        occurrence_max = int(match.group(2))
    else:
        occurrence_min = None
        occurrence_max = None

    result = {
        "bin": folder,
        "occurrence_min": occurrence_min,
        "occurrence_max": occurrence_max,

        "total_predictions": len(df),

        "annotated_validation_smORFs": len(validated),

        "correct_top1": int(
            validated["correct_bool"].sum()
        ),

        "overall_accuracy_pct": round(
            accuracy(validated),
            3
        ),
    }

    # Accuracy for each confidence tier
    for tier in [
        "VERY HIGH",
        "HIGH",
        "MEDIUM",
        "LOW"
    ]:

        tier_df = validated[
            validated["confidence_tier"] == tier
        ]

        prefix = (
            tier
            .lower()
            .replace(" ", "_")
        )

        result[
            f"{prefix}_n"
        ] = len(tier_df)

        result[
            f"{prefix}_accuracy_pct"
        ] = round(
            accuracy(tier_df),
            3
        )

    rows.append(result)

    print(
        f"{folder}: "
        f"Predicted={len(df):,} | "
        f"Validated={len(validated):,} | "
        f"Accuracy={accuracy(validated):.2f}%"
    )


summary = pd.DataFrame(rows)

summary = summary.sort_values(
    by=[
        "occurrence_min",
        "occurrence_max"
    ]
)

summary.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


print()
print("==========================================")
print("ALL BINS COMPLETE")
print("==========================================")

print(f"Bins processed: {len(summary)}")

total_validated = summary[
    "annotated_validation_smORFs"
].sum()

total_correct = summary[
    "correct_top1"
].sum()

overall_accuracy = (
    100
    * total_correct
    / total_validated
)

print(
    f"Total validated smORFs: "
    f"{total_validated:,}"
)

print(
    f"Total correct: "
    f"{total_correct:,}"
)

print(
    f"Overall accuracy: "
    f"{overall_accuracy:.3f}%"
)

print()
print(f"Saved → {OUTPUT}")
PY
