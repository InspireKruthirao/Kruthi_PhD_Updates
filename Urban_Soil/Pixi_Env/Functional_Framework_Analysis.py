#!/usr/bin/env python3
"""
Usage:
    python3 Functional_Framework_Analysis.py                     # submit all bins as PBS jobs
    python3 Functional_Framework_Analysis.py <bin_dir> [out_dir]  # process one bin
"""
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# -config -
PARENT_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL_ge6")
SCRIPT_PATH = Path(__file__).resolve()
OUT_DIR = PARENT_DIR / "results"
LOG_DIR = PARENT_DIR / "logs"
PYTHON_BIN = "/home/n12228516/.conda/envs/kruthi/bin/python3"

MAX_NEIGHBOURS = 5
RANK_WEIGHTS = {1: 1.0, 2: 0.8, 3: 0.6, 4: 0.4, 5: 0.2}  # rank 1 = closest
NAME_RE = re.compile(r"^COG\d+-([A-Za-z]+)$")
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
    match = NAME_RE.match(name or "")
    if not match:
        return []
    return [letter for letter in match.group(1) if letter not in exclude]


def find_target(genes):
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


def accuracy_stats(rows, field):
    evaluated = [r for r in rows if r[field] in ("True", "False")]
    correct = sum(r[field] == "True" for r in evaluated)
    return len(evaluated), correct


# -run one bin -

def run_bin(bin_dir, out_dir):
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

    out_path = out_dir / f"{bin_dir.name}_results.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("\t".join(rows[0].keys()) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row.values()) + "\n")
    print(f"\nWrote {len(rows):,} rows: {out_path}")

    print(f"\n########## [{bin_dir.name}] SUMMARY ##########")
    for field, label in (("correct", "neighbourhood prediction vs annotated targets"),
                          ("correct_excl_SR", "excluding S/R")):
        n, correct = accuracy_stats(rows, field)
        print(f"\n=== ACCURACY ({label}) ===")
        print(f"evaluated: {n}  correct: {correct}  accuracy: {correct / n:.3f}"
              if n else "no annotated smORFs available for evaluation")

    print("\n=== PREDICTED CATEGORY COUNTS ===")
    for value, n in sorted(Counter(r["predicted_category"] for r in rows).items()):
        print(f"  {value:<25} n={n:5d}")
    print(f"  {'TOTAL':<25} n={len(rows):5d}")

    print("\n=== ACCURACY BY CONFIDENCE TIER ===")
    for tier in ("VERY_HIGH", "HIGH", "MEDIUM", "LOW"):
        n, correct = accuracy_stats([r for r in rows if r["confidence_tier"] == tier], "correct")
        if n:
            print(f"  {tier:<10} n={n:5d}  correct={correct:5d}  accuracy={correct / n:.3f}")

    print(f"[{bin_dir.name}] DONE", flush=True)


# -submit one PBS job per bin -

def submit_all():
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    bins = sorted(p for p in PARENT_DIR.glob("SmORF_neighbourhoods_*") if p.is_dir())
    if not bins:
        print(f"No SmORF_neighbourhoods_* folders found under {PARENT_DIR}")
        sys.exit(1)

    print(f"Found {len(bins)} bins. Submitting one job per bin...")
    for bin_path in bins:
        bin_name = bin_path.name
        job_script = f"""#!/bin/bash
#PBS -l select=1:ncpus=2:mem=64gb
#PBS -l walltime=24:00:00
#PBS -q cpu_batch_exec

set -eo pipefail
export PYTHONNOUSERSITE=1

{PYTHON_BIN} -s -u "{SCRIPT_PATH}" "{bin_path}" "{OUT_DIR}"
"""
        result = subprocess.run(
            ["qsub", "-N", f"cog_{bin_name}", "-o", str(LOG_DIR / f"{bin_name}.log"), "-j", "oe"],
            input=job_script, text=True, capture_output=True,
        )
        job_id = result.stdout.strip()
        print(f"  {bin_name} -> {job_id}")

    print("\nAll jobs submitted. Watch progress live with:")
    print(f"  tail -f {LOG_DIR}/*.log")
    print("\nCheck job status with:")
    print("  qstat -u $USER")


def main():
    if len(sys.argv) == 1:
        submit_all()
        return

    bin_dir = Path(sys.argv[1])
    if not bin_dir.is_dir():
        print("Usage:")
        print("  python3 Functional_Framework_Analysis.py                     # submit all bins as PBS jobs")
        print("  python3 Functional_Framework_Analysis.py <bin_dir> [out_dir]  # process one bin")
        sys.exit(1)

    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else OUT_DIR
    run_bin(bin_dir, out_dir)


if __name__ == "__main__":
    main()
