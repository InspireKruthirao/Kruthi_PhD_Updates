#!/usr/bin/env python3
"""Predict Urban Soil smORF COG categories from conserved gene neighbours.

Expected input layout:

    SmORF_neighbourhoods_26_30/
      US_SM.100AA.*/
        SAMPLE_contig_*/
          contig_*.gff

The GFF attributes are expected to contain fields such as:

    ID=US.ORF...;Name=COG0745-T;target=0;Original_ID=...
    ID=US.ORF...;Name=Unknown;target=1;Original_ID=...;Note=TARGET_smORF

The script writes one prediction table and, when necessary, one failure table
inside the selected SmORF neighbourhood bin.
"""

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_BIN_DIR = Path(
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil/"
    "Neighbourhood_Analysis_SQL_ge6/SmORF_neighbourhoods_26_30"
)
URBAN_SOIL_PREFIX = "US_SM.100AA."
DEFAULT_OUTPUT = "smorf_onehot_predictions_validation.tsv"
DEFAULT_FAILURE_OUTPUT = "smorf_failed_reasons.tsv"
UNINFORMATIVE_COGS = {"R", "S"}


COG_DESC = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division and chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall, membrane and envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion and vesicular transport",
    "V": "Defence mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages and transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}

POSITION_WEIGHTS = {
    -5: 0.2,
    -4: 0.4,
    -3: 0.6,
    -2: 0.8,
    -1: 1.0,
    1: 1.0,
    2: 0.8,
    3: 0.6,
    4: 0.4,
    5: 0.2,
}

COLS = [
    "rank",
    "smorf_id",
    "true_cogs",
    "correct",
    "current_annotation",
    "n_occurrences",
    "conservation_pct",
    "confidence_tier",
    "predicted_cog_cat",
    "vote_score_pct",
    "predicted_function",
    "top3_candidates",
    "prediction_set_90pct",
    "top_left_neighbor",
    "top_left_pct",
    "top_right_neighbor",
    "top_right_pct",
    "all_scores",
]

TIER_ORDER = {"VERY HIGH": 0, "HIGH": 1, "MEDIUM": 2, "LOW": 3}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Predict a COG category for each Urban Soil smORF from the GFF "
            "files in one neighbourhood bin."
        )
    )
    parser.add_argument(
        "bin_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_BIN_DIR,
        help="SmORF_neighbourhoods_* directory to analyse",
    )
    parser.add_argument(
        "--prefix",
        default=URBAN_SOIL_PREFIX,
        help=f"smORF directory prefix (default: {URBAN_SOIL_PREFIX})",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"prediction TSV filename (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--failure-output",
        default=DEFAULT_FAILURE_OUTPUT,
        help=f"failure TSV filename (default: {DEFAULT_FAILURE_OUTPUT})",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="print progress after this many smORFs (default: 100)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="process only the first N smORFs; useful for testing",
    )
    return parser.parse_args()


def parse_attributes(raw_attributes):
    """Convert a GFF3 attribute string to a dictionary."""
    attributes = {}
    for item in raw_attributes.split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, value = item.partition("=")
        if separator:
            attributes[key.strip()] = value.strip()
    return attributes


def extract_cog_categories(name):
    """Extract valid COG category letters from COG0123-KLT-style names."""
    if not name or name == "Unknown" or "-" not in name:
        return None
    category_text = name.split("-", 1)[1]
    categories = [category for category in category_text if category in COG_DESC]
    return "".join(dict.fromkeys(categories)) or None


def parse_gff(gff_file):
    """Read one Urban Soil GFF and return genes sorted by genomic start."""
    genes = []
    with gff_file.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise ValueError(
                    f"{gff_file}: line {line_number} has {len(parts)} columns; expected 9"
                )

            attributes = parse_attributes(parts[8])
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError as error:
                raise ValueError(
                    f"{gff_file}: non-integer coordinates on line {line_number}"
                ) from error

            name = attributes.get("Name", "Unknown")
            is_target = attributes.get("target") == "1"
            if not is_target:
                is_target = "TARGET" in attributes.get("Note", "").upper()

            genes.append(
                {
                    "gene_id": attributes.get("ID", "Unknown"),
                    "original_id": attributes.get("Original_ID", parts[0]),
                    "name": name,
                    "cog_cats": extract_cog_categories(name),
                    "strand": parts[6],
                    "start": start,
                    "end": end,
                    "is_target": is_target,
                }
            )

    genes.sort(key=lambda gene: (gene["start"], gene["end"]))
    return genes


def normalize(genes):
    """Orient a neighbourhood so upstream/downstream follow target direction."""
    target_strand = next(
        (gene["strand"] for gene in genes if gene["is_target"]), "+"
    )
    return list(reversed(genes)) if target_strand == "-" else genes


def informative_categories(category_text):
    if not category_text:
        return []
    return [
        category
        for category in category_text
        if category in COG_DESC and category not in UNINFORMATIVE_COGS
    ]


def extract_true_cogs(name):
    categories = informative_categories(extract_cog_categories(name))
    return set(categories) if categories else None


def vote_one_occurrence(normalized_genes, target_index):
    scores = defaultdict(float)
    for index, gene in enumerate(normalized_genes):
        if gene["is_target"]:
            continue

        relative_position = index - target_index
        position_weight = POSITION_WEIGHTS.get(relative_position, 0.0)
        if position_weight == 0:
            continue

        categories = informative_categories(gene["cog_cats"])
        if not categories:
            continue

        weight_each = position_weight / len(categories)
        for category in categories:
            scores[category] += weight_each
    return scores


def predict_one_smorf(smorf_dir):
    smorf_id = smorf_dir.name
    gff_files = sorted(smorf_dir.glob("*/*.gff"))
    if not gff_files:
        return None, "no_gff_files"

    total_scores = defaultdict(float)
    left_counts = Counter()
    right_counts = Counter()
    target_name = "Unknown"
    true_cogs = None
    parsed_occurrences = 0
    invalid_gffs = 0

    for gff_file in gff_files:
        try:
            genes = parse_gff(gff_file)
        except (OSError, ValueError) as error:
            invalid_gffs += 1
            print(f"    WARNING: {error}", file=sys.stderr, flush=True)
            continue

        target_index = next(
            (index for index, gene in enumerate(genes) if gene["is_target"]), None
        )
        if target_index is None:
            continue

        normalized = normalize(genes)
        normalized_target_index = next(
            index for index, gene in enumerate(normalized) if gene["is_target"]
        )

        if target_name == "Unknown":
            target_name = normalized[normalized_target_index]["name"]
            true_cogs = extract_true_cogs(target_name)

        if normalized_target_index > 0:
            left_counts[normalized[normalized_target_index - 1]["name"]] += 1
        if normalized_target_index < len(normalized) - 1:
            right_counts[normalized[normalized_target_index + 1]["name"]] += 1

        for category, score in vote_one_occurrence(
            normalized, normalized_target_index
        ).items():
            total_scores[category] += score
        parsed_occurrences += 1

    if parsed_occurrences == 0:
        if invalid_gffs:
            return None, "all_gffs_invalid_or_missing_target"
        return None, "target_not_found_in_any_gff"
    if not total_scores:
        return None, "no_informative_neighbours"

    total_vote = sum(total_scores.values())
    if total_vote <= 0:
        return None, "zero_total_vote_weight"

    score_pct = {
        category: round(score / total_vote * 100, 1)
        for category, score in total_scores.items()
    }
    ranked = sorted(score_pct.items(), key=lambda item: (-item[1], item[0]))
    top_category, top_score = ranked[0]

    cumulative_score = 0.0
    prediction_set = []
    for category, percentage in ranked:
        prediction_set.append(category)
        cumulative_score += percentage
        if cumulative_score >= 90:
            break

    left_pct = (
        left_counts.most_common(1)[0][1] / parsed_occurrences * 100
        if left_counts
        else 0.0
    )
    right_pct = (
        right_counts.most_common(1)[0][1] / parsed_occurrences * 100
        if right_counts
        else 0.0
    )
    conservation = (left_pct + right_pct) / 2
    top_left = left_counts.most_common(1)[0][0] if left_counts else "None"
    top_right = right_counts.most_common(1)[0][0] if right_counts else "None"

    if top_score >= 50 and conservation >= 90 and parsed_occurrences >= 15:
        tier = "VERY HIGH"
    elif top_score >= 40 and conservation >= 75 and parsed_occurrences >= 8:
        tier = "HIGH"
    elif top_score >= 30 and conservation >= 50 and parsed_occurrences >= 4:
        tier = "MEDIUM"
    else:
        tier = "LOW"

    correct = None if true_cogs is None else top_category in true_cogs
    all_scores = " | ".join(
        f"{category}:{percentage}%" for category, percentage in ranked
    )

    result = {
        "smorf_id": smorf_id,
        "current_annotation": target_name,
        "true_cogs": ",".join(sorted(true_cogs)) if true_cogs else None,
        "correct": correct,
        "n_occurrences": parsed_occurrences,
        "conservation_pct": round(conservation, 1),
        "confidence_tier": tier,
        "predicted_cog_cat": top_category,
        "vote_score_pct": top_score,
        "predicted_function": COG_DESC[top_category],
        "top3": " | ".join(
            f"{category}:{percentage}%" for category, percentage in ranked[:3]
        ),
        "pred_set_90pct": ",".join(prediction_set),
        "top_left_neighbor": top_left,
        "top_left_pct": round(left_pct, 1),
        "top_right_neighbor": top_right,
        "top_right_pct": round(right_pct, 1),
        "all_scores": all_scores,
    }
    return result, None


def print_validation_summary(results, label=""):
    validated = [result for result in results if result["true_cogs"] is not None]
    header = f"VALIDATION{' - ' + label if label else ''}"
    print(f"\n{header}")
    print("-" * len(header))
    if not validated:
        print(
            "No targets have an existing usable COG annotation. Predictions were "
            "still generated, but accuracy cannot be calculated."
        )
        return

    correct_count = sum(result["correct"] is True for result in validated)
    print(f"Annotated smORFs : {len(validated)}")
    print(f"Correct          : {correct_count}")
    print(f"Accuracy         : {correct_count / len(validated):.3%}")

    print("\nACCURACY BY TIER")
    print("----------------")
    for tier in ("VERY HIGH", "HIGH", "MEDIUM", "LOW"):
        subset = [
            result for result in validated if result["confidence_tier"] == tier
        ]
        if subset:
            tier_correct = sum(result["correct"] is True for result in subset)
            print(
                f"{tier:<10}  N={len(subset):5d}  "
                f"Acc={tier_correct / len(subset):.3%}"
            )


def write_tsv(results, output_file):
    sorted_results = sorted(
        results,
        key=lambda result: (
            TIER_ORDER.get(result["confidence_tier"], 4),
            -result["vote_score_pct"],
            -result["n_occurrences"],
            result["smorf_id"],
        ),
    )

    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLS)
        for rank, result in enumerate(sorted_results, start=1):
            writer.writerow(
                [
                    rank,
                    result["smorf_id"],
                    result["true_cogs"] or "",
                    "" if result["correct"] is None else result["correct"],
                    result["current_annotation"],
                    result["n_occurrences"],
                    result["conservation_pct"],
                    result["confidence_tier"],
                    result["predicted_cog_cat"],
                    result["vote_score_pct"],
                    result["predicted_function"],
                    result["top3"],
                    result["pred_set_90pct"],
                    result["top_left_neighbor"],
                    result["top_left_pct"],
                    result["top_right_neighbor"],
                    result["top_right_pct"],
                    result["all_scores"],
                ]
            )
    return sorted_results


def write_failed_tsv(failures, output_file):
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["smorf_id", "reason"])
        writer.writerows(failures)


def discover_smorf_dirs(bin_dir, prefix):
    """Find Urban Soil smORF directories; GFF absence is reported as a failure."""
    return sorted(
        (
            entry
            for entry in bin_dir.iterdir()
            if entry.is_dir() and entry.name.startswith(prefix)
        ),
        key=lambda entry: entry.name,
    )


def run_one_bin_dir(
    bin_dir,
    prefix,
    output_name,
    failure_output_name,
    progress_every,
    limit=None,
):
    smorf_dirs = discover_smorf_dirs(bin_dir, prefix)
    if limit is not None:
        smorf_dirs = smorf_dirs[:limit]

    if not smorf_dirs:
        print(f"ERROR: no {prefix}* directories found in {bin_dir}", file=sys.stderr)
        return [], [], 2

    print(f"smORFs found: {len(smorf_dirs):,}")
    results = []
    failures = []

    for number, smorf_dir in enumerate(smorf_dirs, start=1):
        result, reason = predict_one_smorf(smorf_dir)
        if result is not None:
            results.append(result)
        else:
            failures.append((smorf_dir.name, reason))

        if number % progress_every == 0 or number == len(smorf_dirs):
            print(
                f"  {number:,}/{len(smorf_dirs):,} processed; "
                f"predicted={len(results):,}; failed={len(failures):,}",
                flush=True,
            )

    prediction_file = bin_dir / output_name
    write_tsv(results, prediction_file)
    print(f"Predictions: {prediction_file}")

    failure_file = bin_dir / failure_output_name
    if failures:
        write_failed_tsv(failures, failure_file)
        print(f"Failures: {failure_file}")
        print("Failure reasons:")
        for reason, count in Counter(reason for _, reason in failures).most_common():
            print(f"  {reason}: {count:,}")
    elif failure_file.exists():
        failure_file.unlink()

    print_validation_summary(results, label=bin_dir.name)
    return results, failures, 0


def main():
    args = parse_args()
    bin_dir = args.bin_dir.resolve()

    if not bin_dir.is_dir():
        print(f"ERROR: bin directory does not exist: {bin_dir}", file=sys.stderr)
        return 1
    if args.progress_every < 1:
        print("ERROR: --progress-every must be at least 1", file=sys.stderr)
        return 1
    if args.limit is not None and args.limit < 1:
        print("ERROR: --limit must be at least 1", file=sys.stderr)
        return 1

    print("=" * 70)
    print(f"Processing: {bin_dir}")
    print(f"Urban Soil prefix: {args.prefix}")
    print("=" * 70)

    results, failures, status = run_one_bin_dir(
        bin_dir=bin_dir,
        prefix=args.prefix,
        output_name=args.output,
        failure_output_name=args.failure_output,
        progress_every=args.progress_every,
        limit=args.limit,
    )
    if status:
        return status

    print("\nSUMMARY")
    print(f"Predicted: {len(results):,}")
    print(f"Failed:    {len(failures):,}")
    print("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

