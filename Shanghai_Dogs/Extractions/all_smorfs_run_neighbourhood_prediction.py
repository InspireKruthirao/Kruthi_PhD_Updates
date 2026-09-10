#!/usr/bin/env python3

import re
import sys
from pathlib import Path
from collections import defaultdict, Counter


# Bin directory passed from PBS
bin_dir = Path(sys.argv[1])


COG_DESC = {
    "C": "Energy production and conversion",
    "D": "Cell cycle control and division",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation and ribosome",
    "K": "Transcription",
    "L": "DNA replication and repair",
    "M": "Cell wall and membrane",
    "N": "Cell motility",
    "O": "Protein turnover and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis",
    "R": "General function prediction",
    "S": "Function unknown",
    "T": "Signal transduction",
    "U": "Intracellular trafficking",
    "V": "Defense mechanisms",
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


def parse_gff(gff_file):
    genes = []

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")

            if len(parts) < 9:
                continue

            attrs = parts[8]

            name_m = re.search(
                r"Name=([^;]+)",
                attrs
            )

            is_target = "target=1" in attrs

            name = (
                name_m.group(1)
                if name_m
                else "Unknown"
            )

            strand = parts[6]
            start = int(parts[3])

            cog_cats = (
                name.split("-")[1]
                if name != "Unknown" and "-" in name
                else None
            )

            genes.append({
                "name": name,
                "cog_cats": cog_cats,
                "strand": strand,
                "start": start,
                "is_target": is_target,
            })

    genes.sort(
        key=lambda gene: gene["start"]
    )

    return genes


def normalize(genes):
    target_strand = next(
        (
            gene["strand"]
            for gene in genes
            if gene["is_target"]
        ),
        "+"
    )

    if target_strand == "-":
        return list(reversed(genes)), True

    return genes, False


def extract_true_cogs(name):
    if not name or "-" not in name:
        return None

    cog_part = name.split("-", 1)[1]

    valid = {
        cog
        for cog in cog_part
        if cog in COG_DESC
        and cog not in ("S", "R")
    }

    return valid if valid else None


def vote_one_occurrence(norm_genes, target_idx):
    scores = defaultdict(float)

    for i, gene in enumerate(norm_genes):
        if gene["is_target"]:
            continue

        rel_pos = i - target_idx

        pos_weight = POSITION_WEIGHTS.get(
            rel_pos,
            0.0
        )

        if pos_weight == 0:
            continue

        cog_cats = gene["cog_cats"]

        if not cog_cats:
            continue

        if set(cog_cats) <= {"S", "R"}:
            continue

        valid = [
            cog
            for cog in cog_cats
            if cog in COG_DESC
            and cog not in ("S", "R")
        ]

        if not valid:
            continue

        weight_each = (
            pos_weight / len(valid)
        )

        for category in valid:
            scores[category] += weight_each

    return scores


def predict_one_smorf(smorf_dir):
    smorf_id = smorf_dir.name

    gff_files = sorted(
        smorf_dir.glob("*/*.gff")
    )

    if not gff_files:
        return None, "no_gff_files"

    total_scores = defaultdict(float)
    left_counts = Counter()
    right_counts = Counter()

    target_name = "Unknown"
    true_cogs = None
    n_parsed = 0

    for gff_file in gff_files:
        genes = parse_gff(gff_file)

        target_idx = next(
            (
                i
                for i, gene in enumerate(genes)
                if gene["is_target"]
            ),
            None
        )

        if target_idx is None:
            continue

        norm_genes, _ = normalize(genes)

        norm_target_idx = next(
            i
            for i, gene in enumerate(norm_genes)
            if gene["is_target"]
        )

        if target_name == "Unknown":
            target_name = (
                norm_genes[norm_target_idx]["name"]
            )

            true_cogs = extract_true_cogs(
                target_name
            )

        if norm_target_idx > 0:
            left_name = (
                norm_genes[
                    norm_target_idx - 1
                ]["name"]
            )

            left_counts[left_name] += 1

        if norm_target_idx < len(norm_genes) - 1:
            right_name = (
                norm_genes[
                    norm_target_idx + 1
                ]["name"]
            )

            right_counts[right_name] += 1

        scores = vote_one_occurrence(
            norm_genes,
            norm_target_idx
        )

        for category, score in scores.items():
            total_scores[category] += score

        n_parsed += 1

    if n_parsed == 0:
        return None, "target_not_found_in_any_gff"

    if not total_scores:
        return None, "no_informative_neighbours"

    total_vote = sum(
        total_scores.values()
    )

    if total_vote <= 0:
        return None, "zero_total_vote_weight"

    score_pct = {
        category: round(
            score / total_vote * 100,
            1
        )
        for category, score
        in total_scores.items()
    }

    ranked = sorted(
        score_pct.items(),
        key=lambda x: -x[1]
    )

    top_cat, top_score = ranked[0]

    top3_str = " | ".join(
        f"{category}:{score}%"
        for category, score
        in ranked[:3]
    )

    cumulative = 0.0
    prediction_set = []

    for category, score in ranked:
        prediction_set.append(category)
        cumulative += score

        if cumulative >= 90:
            break

    left_pct = (
        left_counts.most_common(1)[0][1]
        / n_parsed
        * 100
        if left_counts
        else 0.0
    )

    right_pct = (
        right_counts.most_common(1)[0][1]
        / n_parsed
        * 100
        if right_counts
        else 0.0
    )

    conservation = (
        left_pct + right_pct
    ) / 2

    top_left = (
        left_counts.most_common(1)[0][0]
        if left_counts
        else "None"
    )

    top_right = (
        right_counts.most_common(1)[0][0]
        if right_counts
        else "None"
    )

    if (
        top_score >= 50
        and conservation >= 90
        and n_parsed >= 15
    ):
        tier = "VERY HIGH"

    elif (
        top_score >= 40
        and conservation >= 75
        and n_parsed >= 8
    ):
        tier = "HIGH"

    elif (
        top_score >= 30
        and conservation >= 50
        and n_parsed >= 4
    ):
        tier = "MEDIUM"

    else:
        tier = "LOW"

    correct = (
        true_cogs is not None
        and top_cat in true_cogs
    )

    result = {
        "smorf_id": smorf_id,
        "current_annotation": target_name,
        "true_cogs": (
            ",".join(sorted(true_cogs))
            if true_cogs
            else None
        ),
        "correct": correct,
        "n_occurrences": n_parsed,
        "conservation_pct": round(
            conservation,
            1
        ),
        "confidence_tier": tier,
        "predicted_cog_cat": top_cat,
        "vote_score_pct": top_score,
        "predicted_function": COG_DESC.get(
            top_cat,
            "?"
        ),
        "top3": top3_str,
        "pred_set_90pct": ",".join(
            prediction_set
        ),
        "top_left_neighbor": top_left,
        "top_left_pct": round(
            left_pct,
            1
        ),
        "top_right_neighbor": top_right,
        "top_right_pct": round(
            right_pct,
            1
        ),
        "all_scores": " | ".join(
            f"{category}:{score}%"
            for category, score
            in ranked
        ),
    }

    return result, None


def print_validation_summary(results, label=""):
    validated = [
        result
        for result in results
        if result["true_cogs"] is not None
    ]

    print(f"\nVALIDATION — {label}")
    print("-" * 50)

    if not validated:
        print("No annotated smORFs found.")
        return

    correct = sum(
        1
        for result in validated
        if result["correct"]
    )

    accuracy = correct / len(validated)

    print(
        f"Annotated smORFs : {len(validated)}"
    )

    print(
        f"Correct          : {correct}"
    )

    print(
        f"Accuracy         : {accuracy:.3%}"
    )

    print("\nACCURACY BY TIER")
    print("----------------")

    for tier in [
        "VERY HIGH",
        "HIGH",
        "MEDIUM",
        "LOW"
    ]:
        subset = [
            result
            for result in validated
            if result["confidence_tier"] == tier
        ]

        if not subset:
            continue

        tier_correct = sum(
            1
            for result in subset
            if result["correct"]
        )

        tier_accuracy = (
            tier_correct / len(subset)
        )

        print(
            f"{tier:<10} "
            f"N={len(subset):5d} "
            f"Acc={tier_accuracy:.3%}"
        )


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


TIER_ORDER = {
    "VERY HIGH": 0,
    "HIGH": 1,
    "MEDIUM": 2,
    "LOW": 3,
}


def write_tsv(results, out_file):
    results_sorted = sorted(
        results,
        key=lambda result: (
            TIER_ORDER.get(
                result["confidence_tier"],
                4
            ),
            -result["vote_score_pct"],
            -result["n_occurrences"],
        )
    )

    with out_file.open("w") as f:
        f.write(
            "\t".join(COLS) + "\n"
        )

        for rank, result in enumerate(
            results_sorted,
            1
        ):
            row = [
                rank,
                result["smorf_id"],
                result["true_cogs"],
                result["correct"],
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

            f.write(
                "\t".join(
                    str(value)
                    for value in row
                )
                + "\n"
            )

    return results_sorted


def write_failed_tsv(failed, out_file):
    with out_file.open("w") as f:
        f.write(
            "smorf_id\treason\n"
        )

        for name, reason in failed:
            f.write(
                f"{name}\t{reason}\n"
            )


def run_one_bin_dir(bin_path):
    smorf_dirs = sorted([
        directory
        for directory in bin_path.iterdir()
        if (
            directory.is_dir()
            and directory.name.startswith(
                "SHD1_SM"
            )
            and any(
                directory.glob("*/*.gff")
            )
        )
    ])

    if not smorf_dirs:
        print(
            "[skip] No smORF directories "
            "with GFF files."
        )
        return []

    print(
        f"smORFs found : "
        f"{len(smorf_dirs):,}"
    )

    results = []
    failed = []

    for number, smorf_dir in enumerate(
        smorf_dirs,
        1
    ):
        result, reason = predict_one_smorf(
            smorf_dir
        )

        if result:
            results.append(result)

        else:
            failed.append(
                (smorf_dir.name, reason)
            )

        if (
            number % 1000 == 0
            or number == len(smorf_dirs)
        ):
            print(
                f"{number:>8,} / "
                f"{len(smorf_dirs):,}   "
                f"predicted: {len(results):,}   "
                f"failed: {len(failed):,}",
                flush=True
            )

    out_file = (
        bin_path
        / "smorf_onehot_predictions_validation.tsv"
    )

    write_tsv(
        results,
        out_file
    )

    print(
        f"Saved → {out_file}"
    )

    if failed:
        failed_file = (
            bin_path
            / "smorf_failed_reasons.tsv"
        )

        write_failed_tsv(
            failed,
            failed_file
        )

    print_validation_summary(
        results,
        label=bin_path.name
    )

    return results


def main():
    if not bin_dir.exists():
        print(
            f"ERROR: {bin_dir} does not exist"
        )

        sys.exit(1)

    print("=" * 60)
    print(
        f"Processing: {bin_dir.name}"
    )
    print("=" * 60)

    results = run_one_bin_dir(
        bin_dir
    )

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(
        f"Total smORFs predicted : "
        f"{len(results):,}"
    )


if __name__ == "__main__":
    main()
