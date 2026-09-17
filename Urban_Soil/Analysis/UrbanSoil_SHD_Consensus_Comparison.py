#!/usr/bin/env python3
"""Compare Shanghai Dogs immediate-neighbour rules with Urban Soil predictions.

This is a read-only validation experiment. It consumes the existing Urban Soil
prediction TSVs and does not change GFF files or overwrite predictions.
"""

from pathlib import Path
import csv


ROOT = Path(
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil/"
    "Neighbourhood_Analysis_SQL_ge6"
)
BINS = (
    "6_10", "11_15", "16_20", "21_25", "26_30", "31_35",
    "36_40", "41_45", "46_50", "51_55", "56_60",
)
VALID_COGS = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ") - {"R", "S"}


def cog_set(value):
    """Extract informative COG category letters from a GFF Name value."""
    value = (value or "").strip()
    if value in {"", "None", "Unknown"} or "-" not in value:
        return set()
    category_text = value.split("-", 1)[1]
    return {letter for letter in category_text if letter in VALID_COGS}


def true_set(value):
    value = (value or "").strip()
    if value in {"", "None", "nan", "-"}:
        return set()
    return {letter for letter in value.replace(",", "") if letter in VALID_COGS}


def choose_consensus(left, right, weighted):
    """Use agreed neighbour category; use weighted top-1 to break multi-label ties."""
    agreed = left & right
    if not agreed:
        return None
    if weighted in agreed:
        return weighted
    return sorted(agreed)[0]


def pct(numerator, denominator):
    return 100.0 * numerator / denominator if denominator else 0.0


def main():
    counts = {
        name: {"eligible": 0, "correct": 0}
        for name in (
            "weighted_all", "left_majority", "right_majority",
            "neighbour_consensus", "hybrid_consensus_else_weighted",
        )
    }
    annotated = 0
    files = 0

    for label in BINS:
        path = ROOT / f"SmORF_neighbourhoods_{label}" / "smorf_onehot_predictions_validation.tsv"
        if not path.is_file():
            print(f"WARNING: missing {path}")
            continue
        files += 1
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                truth = true_set(row.get("true_cogs"))
                if not truth:
                    continue
                annotated += 1

                weighted = (row.get("predicted_cog_cat") or "").strip()
                left = cog_set(row.get("top_left_neighbor"))
                right = cog_set(row.get("top_right_neighbor"))
                consensus = choose_consensus(left, right, weighted)

                counts["weighted_all"]["eligible"] += 1
                counts["weighted_all"]["correct"] += weighted in truth

                if left:
                    counts["left_majority"]["eligible"] += 1
                    counts["left_majority"]["correct"] += bool(left & truth)

                if right:
                    counts["right_majority"]["eligible"] += 1
                    counts["right_majority"]["correct"] += bool(right & truth)

                if consensus is not None:
                    counts["neighbour_consensus"]["eligible"] += 1
                    counts["neighbour_consensus"]["correct"] += consensus in truth

                hybrid = consensus if consensus is not None else weighted
                counts["hybrid_consensus_else_weighted"]["eligible"] += 1
                counts["hybrid_consensus_else_weighted"]["correct"] += hybrid in truth

    print(f"Prediction files read: {files}/{len(BINS)}")
    print(f"Annotated smORFs: {annotated:,}\n")
    print(f"{'Method':38s} {'Eligible':>10s} {'Coverage':>10s} {'Correct':>10s} {'Accuracy':>10s}")
    print("-" * 84)
    for method, values in counts.items():
        eligible = values["eligible"]
        correct = values["correct"]
        print(
            f"{method:38s} {eligible:10,d} {pct(eligible, annotated):9.3f}% "
            f"{correct:10,d} {pct(correct, eligible):9.3f}%"
        )

    print("\nInterpretation:")
    print("- Compare hybrid accuracy with weighted_all at 100% coverage.")
    print("- Consensus can be more accurate but usually covers fewer smORFs.")
    print("- This is a quick screen using stored majority-neighbour names, not final model tuning.")


if __name__ == "__main__":
    main()

