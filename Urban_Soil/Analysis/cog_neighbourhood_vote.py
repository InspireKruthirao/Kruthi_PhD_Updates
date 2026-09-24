#!/usr/bin/env python3
"""
cog_neighbourhood_vote.py

USAGE
-----
    python3 cog_neighbourhood_vote.py \
        --base-dir SmORF_neighbourhoods_26_30 \
        --out-summary results_26_30_summary.tsv
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

RANK_WEIGHTS = {1: 1.0, 2: 0.8, 3: 0.6, 4: 0.4, 5: 0.2}
NAME_RE = re.compile(r"^COG\d+-([A-Za-z]+)$")
TIE_TOLERANCE = 1e-9
EXCLUDE_FROM_EXCL_SR = {"S", "R"}


def parse_attributes(attr_str):
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        k, v = field.split("=", 1)
        attrs[k] = v
    return attrs


def parse_gff(gff_path):
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            attrs = parse_attributes(f[8])
            genes.append(
                {
                    "seqid": f[0],
                    "start": int(f[3]),
                    "end": int(f[4]),
                    "strand": f[6],
                    "id": attrs.get("ID"),
                    "name": attrs.get("Name"),
                    "target": attrs.get("target") == "1",
                    "note": attrs.get("Note", ""),
                }
            )
    genes.sort(key=lambda g: g["start"])
    return genes


def get_cog_categories(name, exclude=frozenset()):
    """'COG3332-S' -> ['S']; 'COG0673-CO' -> ['C','O'] (minus any excluded letters)."""
    if not name or name == "Unknown":
        return []
    m = NAME_RE.match(name)
    if not m:
        return []
    letters = [c for c in m.group(1) if c not in exclude]
    return letters


def find_target(genes):
    for g in genes:
        if g["note"] == "TARGET_smORF":
            return g
    for g in genes:
        if g["target"]:
            return g
    return None


def strand_relative_order(genes, target):
    """Reverse gene order if target is on '-' strand, so list order always
    runs in the direction of transcription. Distances still use real
    genomic coordinates, unaffected by list order."""
    if target["strand"] == "-":
        return list(reversed(genes))
    return genes


def _distance(target, g):
    """Genomic intergenic gap regardless of list order / strand."""
    if g["start"] >= target["start"]:
        return max(g["start"] - target["end"], 0)
    return max(target["start"] - g["end"], 0)


def ranked_neighbours(genes, target, max_n=5):
    """Per-side adjacency rank (1..max_n each direction), strand-relative."""
    ordered = strand_relative_order(genes, target)
    idx = ordered.index(target)

    out = []
    upstream = ordered[:idx][::-1][:max_n]  # nearest first, reading direction
    for rank, g in enumerate(upstream, start=1):
        dist = _distance(target, g)
        out.append((rank, "upstream", dist, g))

    downstream = ordered[idx + 1:][:max_n]
    for rank, g in enumerate(downstream, start=1):
        dist = _distance(target, g)
        out.append((rank, "downstream", dist, g))
    return out


def find_occurrence_dirs(smorf_dir):
    return sorted(p for p in smorf_dir.iterdir() if p.is_dir())


def resolve(votes):
    if not votes:
        return "NO_EVIDENCE", 0.0
    max_score = max(votes.values())
    top = sorted(c for c, v in votes.items() if abs(v - max_score) < TIE_TOLERANCE)
    if len(top) > 1:
        return "UNRESOLVED", max_score
    return top[0], max_score


def confidence_tier(top_pct):
    if top_pct >= 50:
        return "VERY_HIGH"
    if top_pct >= 40:
        return "HIGH"
    if top_pct >= 30:
        return "MEDIUM"
    return "LOW"


def process_smorf(smorf_dir, max_n=5):
    smorf_id = smorf_dir.name
    occ_dirs = find_occurrence_dirs(smorf_dir)

    votes_all = defaultdict(float)
    votes_excl = defaultdict(float)
    known_all = set()
    n_occurrences_used = 0

    for occ_dir in occ_dirs:
        gff_files = list(occ_dir.glob("*.gff"))
        if not gff_files:
            continue
        genes = parse_gff(gff_files[0])
        target = find_target(genes)
        if target is None:
            continue
        n_occurrences_used += 1

        known_all.update(get_cog_categories(target["name"]))

        for rank, _direction, _dist, g in ranked_neighbours(genes, target, max_n):
            cats_all = get_cog_categories(g["name"])
            cats_excl = get_cog_categories(g["name"], EXCLUDE_FROM_EXCL_SR)
            weight = RANK_WEIGHTS[rank]

            if cats_all:
                share = weight / len(cats_all)
                for c in cats_all:
                    votes_all[c] += share
            if cats_excl:
                share = weight / len(cats_excl)
                for c in cats_excl:
                    votes_excl[c] += share

    prediction_all, score_all = resolve(votes_all)
    prediction_excl, score_excl = resolve(votes_excl)

    total_all = sum(votes_all.values()) or 1.0
    total_excl = sum(votes_excl.values()) or 1.0
    top_pct_all = round(score_all / total_all * 100, 1) if votes_all else 0.0
    top_pct_excl = round(score_excl / total_excl * 100, 1) if votes_excl else 0.0

    known_excl = known_all - EXCLUDE_FROM_EXCL_SR

    def eval_correct(prediction, known_set):
        if not known_set or prediction in ("UNRESOLVED", "NO_EVIDENCE"):
            return "NA"
        return str(prediction in known_set)

    correct_all = eval_correct(prediction_all, known_all)
    correct_excl = eval_correct(prediction_excl, known_excl)

    summary_row = {
        "smorf_id": smorf_id,
        "n_occurrences": len(occ_dirs),
        "n_occurrences_with_target": n_occurrences_used,
        "confidence_tier": confidence_tier(top_pct_all),
        "prediction": prediction_all,
        "top_score_pct": top_pct_all,
        "vote_breakdown": ";".join(f"{c}:{v:.2f}" for c, v in sorted(votes_all.items(), key=lambda x: -x[1])),
        "known_category": ",".join(sorted(known_all)),
        "correct": correct_all,
        "prediction_excl_SR": prediction_excl,
        "top_score_pct_excl_SR": top_pct_excl,
        "known_category_excl_SR": ",".join(sorted(known_excl)),
        "correct_excl_SR": correct_excl,
    }
    return summary_row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--max-neighbours", type=int, default=5)
    args = ap.parse_args()

    smorf_dirs = sorted(p for p in args.base_dir.iterdir() if p.is_dir())
    print(f"Found {len(smorf_dirs)} smORF folders under {args.base_dir}", file=sys.stderr)

    summary_rows = []
    for i, smorf_dir in enumerate(smorf_dirs, start=1):
        if i % 200 == 0:
            print(f"  ... {i}/{len(smorf_dirs)} processed", file=sys.stderr)
        summary_rows.append(process_smorf(smorf_dir, args.max_neighbours))

    if summary_rows:
        with open(args.out_summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)

    def report(field, label):
        evaluated = [r for r in summary_rows if r[field] in ("True", "False")]
        correct = sum(1 for r in evaluated if r[field] == "True")
        print(f"\n=== ACCURACY ({label}) ===", file=sys.stderr)
        if evaluated:
            print(f"evaluated: {len(evaluated)}  correct: {correct}  accuracy: {correct/len(evaluated):.3f}",
                  file=sys.stderr)
        else:
            print("no known categories available to evaluate against", file=sys.stderr)
        return evaluated

    print(f"\nsmORFs total: {len(summary_rows)}", file=sys.stderr)
    report("correct", "ALL categories, matches original 7-step spec")
    report("correct_excl_SR", "excluding S/R - narrower metric, see docstring")

    print("\n=== ACCURACY BY CONFIDENCE TIER (all categories) ===", file=sys.stderr)
    for tier in ["VERY_HIGH", "HIGH", "MEDIUM", "LOW"]:
        subset = [r for r in summary_rows if r["confidence_tier"] == tier and r["correct"] in ("True", "False")]
        if not subset:
            continue
        acc = sum(1 for r in subset if r["correct"] == "True") / len(subset)
        print(f"  {tier:<10} n={len(subset):5d}  accuracy={acc:.3f}", file=sys.stderr)

    print(f"\nWrote: {args.out_summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
