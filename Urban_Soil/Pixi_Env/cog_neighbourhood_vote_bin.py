#!/usr/bin/env python3
"""
Predict the COG functional category of every smORF from its neighbours --
ONE BIN at a time, so it can be run as a separate PBS job per bin.

Same prediction logic as cog_neighbourhood_vote.py (rank-weighted vote,
1.0/0.8/0.6/0.4/0.2, up to 5 neighbours per side), unchanged. What's
different is purely operational:

  - takes ONE SmORF_neighbourhoods_X_Y folder (not the whole parent dir),
    so several bins can run as separate, simultaneous PBS jobs instead of
    one long serial job over everything
  - prints progress every PROGRESS_EVERY smORFs (flush=True), so a PBS
    log file shows live movement instead of staying silent until the end
    -- watch it with: tail -f logs/<bin>.log
  - writes the full per-smORF table to a TSV file instead of printing it
    (printing millions of rows to a log file is slow and not useful) --
    only the summary blocks (accuracy, counts, confidence tiers) print to
    the log, for a quick glance while it's running or after

Usage:
    python3 cog_neighbourhood_vote_bin.py <bin_dir> [out_dir]

    <bin_dir>  path to one SmORF_neighbourhoods_X_Y folder (required)
    [out_dir]  where to write <bin_name>_results.tsv (default: bin_dir's
               parent / "results")
"""

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

MAX_NEIGHBOURS = 5  # neighbours per side per occurrence (1-5)
RANK_WEIGHTS = {1: 1.0, 2: 0.8, 3: 0.6, 4: 0.4, 5: 0.2}  # rank 1 = closest
NAME_RE = re.compile(r"^COG\d+-([A-Za-z]+)$")  # e.g. COG0673-CO
TIE_TOLERANCE = 1e-9
SR = {"S", "R"}
PROGRESS_EVERY = 500


def parse_gff(path):
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
    ordered = genes[::-1] if target["strand"] == "-" else genes
    i = ordered.index(target)
    for side in (ordered[:i][::-1][:max_n], ordered[i + 1:][:max_n]):
        yield from enumerate(side, start=1)


def resolve(votes):
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
    if not known or prediction in ("UNRESOLVED", "NO_EVIDENCE"):
        return "NA"
    return str(prediction in known)


def process_smorf(smorf_dir):
    votes = defaultdict(float)
    votes_no_sr = defaultdict(float)  # secondary comparison without S and R
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

        for rank, gene in neighbour_genes(genes, target, MAX_NEIGHBOURS):
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


def write_tsv(rows, out_path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("\t".join(rows[0].keys()) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row.values()) + "\n")
    print(f"\nWrote {len(rows):,} rows: {out_path}")


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
    if len(sys.argv) < 2:
        print("Usage: python3 cog_neighbourhood_vote_bin.py <bin_dir> [out_dir]")
        sys.exit(1)

    bin_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else bin_dir.parent / "results"

    if not bin_dir.is_dir():
        print(f"Bin folder not found: {bin_dir}")
        sys.exit(1)

    smorf_dirs = sorted(p for p in bin_dir.iterdir() if p.is_dir())
    total = len(smorf_dirs)
    print(f"[{bin_dir.name}] smORFs to process: {total:,}", flush=True)
    if not smorf_dirs:
        print(f"[{bin_dir.name}] No smORF folders found.")
        return

    rows = []
    for number, d in enumerate(smorf_dirs, 1):
        rows.append(process_smorf(d))
        if number % PROGRESS_EVERY == 0 or number == total:
            print(f"[{bin_dir.name}] {number:,}/{total:,} smORFs processed", flush=True)

    write_tsv(rows, out_dir / f"{bin_dir.name}_results.tsv")

    print(f"\n########## [{bin_dir.name}] SUMMARY ##########")
    report_accuracy(rows, "correct", "neighbourhood prediction vs annotated targets")
    report_accuracy(rows, "correct_excl_SR", "excluding S/R")
    report_counts(rows, "predicted_category", "PREDICTED CATEGORY COUNTS")
    report_confidence_accuracy(rows)
    print(f"[{bin_dir.name}] DONE", flush=True)


if __name__ == "__main__":
    main()
