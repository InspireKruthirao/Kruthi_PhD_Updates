#!/usr/bin/env python3
"""
Usage:
    python3 cog_neighbourhood_vote.py --base-dir SmORF_neighbourhoods_26_30
"""

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

RANK_WEIGHTS = {1: 1.0, 2: 0.8, 3: 0.6, 4: 0.4, 5: 0.2}  # rank 1 = closest
NAME_RE = re.compile(r"^COG\d+-([A-Za-z]+)$")  # e.g. COG0673-CO
TIE_TOLERANCE = 1e-9
SR = {"S", "R"}


def parse_gff(path):
    """Return the genes in a GFF file, sorted by start coordinate."""
    genes = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) != 9:
                continue
            try:
                start = int(fields[3])
            except ValueError:
                print(f"Warning: bad coordinates in {path}; skipping line")
                continue
            attrs = dict(
                kv.strip().split("=", 1) for kv in fields[8].split(";") if "=" in kv
            )
            genes.append({
                "start": start,
                "strand": fields[6],
                "name": attrs.get("Name"),
                "target": attrs.get("target") == "1",
                "note": attrs.get("Note", ""),
            })
    return sorted(genes, key=lambda gene: gene["start"])


def get_cog_categories(name, exclude=()):
    """COG3332-S -> ['S'];  COG0673-CO -> ['C', 'O']."""
    match = NAME_RE.match(name or "")
    if not match:
        return []
    return [letter for letter in match.group(1) if letter not in exclude]


def find_target(genes):
    """Prefer Note=TARGET_smORF, fall back to target=1."""
    for gene in genes:
        if gene["note"] == "TARGET_smORF":
            return gene
    return next((gene for gene in genes if gene["target"]), None)


def neighbour_genes(genes, target, max_n):
    """Yield (rank, gene) for the closest genes on each side of the target."""
    ordered = genes[::-1] if target["strand"] == "-" else genes
    i = ordered.index(target)
    for side in (ordered[:i][::-1][:max_n], ordered[i + 1:][:max_n]):
        yield from enumerate(side, start=1)


def resolve(votes):
    """Return (winning category, its % of all votes)."""
    if not votes:
        return "NO_EVIDENCE", 0.0
    top = max(votes.values())
    winners = sorted(c for c, score in votes.items() if abs(score - top) < TIE_TOLERANCE)
    category = winners[0] if len(winners) == 1 else "UNRESOLVED"
    return category, round(top / sum(votes.values()) * 100, 1)


def confidence_tier(pct):
    if pct >= 50:
        return "VERY_HIGH"
    if pct >= 40:
        return "HIGH"
    if pct >= 30:
        return "MEDIUM"
    return "LOW"


def evaluate(prediction, known):
    """'True'/'False' if the prediction can be checked, otherwise 'NA'."""
    if not known or prediction in ("UNRESOLVED", "NO_EVIDENCE"):
        return "NA"
    return str(prediction in known)


def process_smorf(smorf_dir, max_n):
    votes = defaultdict(float)         # all categories
    votes_no_sr = defaultdict(float)   # secondary comparison without S and R
    known = set()
    occurrences = [p for p in smorf_dir.iterdir() if p.is_dir()]
    n_with_target = 0

    for occurrence in occurrences:
        gff_files = sorted(occurrence.glob("*.gff"))
        if not gff_files:
            continue
        genes = parse_gff(gff_files[0])
        target = find_target(genes)
        if target is None:
            continue

        n_with_target += 1
        known.update(get_cog_categories(target["name"]))

        for rank, gene in neighbour_genes(genes, target, max_n):
            for tally, exclude in ((votes, ()), (votes_no_sr, SR)):
                categories = get_cog_categories(gene["name"], exclude)
                for category in categories:
                    tally[category] += RANK_WEIGHTS[rank] / len(categories)

    predicted, pct = resolve(votes)
    predicted_no_sr, pct_no_sr = resolve(votes_no_sr)
    known_no_sr = known - SR

    breakdown = ";".join(
        f"{c}:{s:.2f}" for c, s in sorted(votes.items(), key=lambda kv: (-kv[1], kv[0]))
    )

    return {
        "smorf_id": smorf_dir.name,
        "n_occurrences": len(occurrences),
        "n_occurrences_with_target": n_with_target,
        "target_annotation_status": "ANNOTATED" if known else "UNANNOTATED",
        "known_category": ",".join(sorted(known)),
        "predicted_category": predicted,
        "top_score_pct": pct,
        "confidence_tier": confidence_tier(pct),
        "vote_breakdown": breakdown,
        "correct": evaluate(predicted, known),
        "prediction_excl_SR": predicted_no_sr,
        "top_score_pct_excl_SR": pct_no_sr,
        "known_category_excl_SR": ",".join(sorted(known_no_sr)),
        "correct_excl_SR": evaluate(predicted_no_sr, known_no_sr),
    }


def print_table(rows):
    print("\n=== RESULTS (one row per smORF) ===")
    print("\t".join(rows[0]))
    for row in rows:
        print("\t".join(str(value) for value in row.values()))


def report_accuracy(rows, field, label):
    evaluated = [r for r in rows if r[field] in ("True", "False")]
    correct = sum(r[field] == "True" for r in evaluated)
    print(f"\n=== ACCURACY ({label}) ===")
    if evaluated:
        print(f"evaluated: {len(evaluated)}  correct: {correct}  "
              f"accuracy: {correct / len(evaluated):.3f}")
    else:
        print("no annotated smORFs available for evaluation")


def report_counts(rows, field, title):
    print(f"\n=== {title} ===")
    for value, n in sorted(Counter(r[field] for r in rows).items()):
        print(f"  {value:<25} n={n:5d}")
    print(f"  {'TOTAL':<25} n={len(rows):5d}")


def report_confidence_accuracy(rows):
    print("\n=== ACCURACY BY CONFIDENCE TIER ===")
    for tier in ("VERY_HIGH", "HIGH", "MEDIUM", "LOW"):
        subset = [r for r in rows
                  if r["confidence_tier"] == tier and r["correct"] in ("True", "False")]
        if subset:
            correct = sum(r["correct"] == "True" for r in subset)
            print(f"  {tier:<10} n={len(subset):5d}  correct={correct:5d}  "
                  f"accuracy={correct / len(subset):.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", required=True, type=Path,
                        help="Directory containing smORF folders.")
    parser.add_argument("--max-neighbours", type=int, default=5,
                        help="Neighbours per side per occurrence (1-5, default: 5).")
    args = parser.parse_args()

    if not 1 <= args.max_neighbours <= max(RANK_WEIGHTS):
        parser.error(f"--max-neighbours must be between 1 and {max(RANK_WEIGHTS)}")
    if not args.base_dir.is_dir():
        parser.error(f"Base directory not found: {args.base_dir}")

    smorf_dirs = sorted(p for p in args.base_dir.iterdir() if p.is_dir())
    rows = [process_smorf(d, args.max_neighbours) for d in smorf_dirs]
    if not rows:
        parser.error(f"No smORF folders found in {args.base_dir}")

    print_table(rows)
    report_accuracy(rows, "correct", "neighbourhood prediction vs annotated targets")
    report_accuracy(rows, "correct_excl_SR", "excluding S/R")
    report_counts(rows, "predicted_category", "PREDICTED CATEGORY COUNTS")
    report_confidence_accuracy(rows)


if __name__ == "__main__":
    main()
