#!/usr/bin/env python3
"""
ONE file: COG category prediction for every smORF, with many method
combinations, from smorfs_all_bins.sqlite. Loads the database ONCE, then
every combination is pure in-memory work -- no small files touched.

    python3 sweep_db.py [db]

    db   default: smorfs_all_bins.sqlite (in the current folder)

Writes:
    sweep_summary.csv      one row per (combination, bin) + an ALL row,
                           best overall accuracy first
    sweep_predictions.tsv  one row per smORF, one column per combination:
                           the predicted category from every method

"linear, 5 nbrs, all-cats (= original)" is exactly what
Functional_Framework_Analysis.py does, so its ALL row should match your
existing numbers. CONTROL rows ignore the neighbours entirely -- a real
method has to beat them to show the neighbourhood carries signal.

To add ideas: edit the loops in build_combinations(). Nothing else changes.
"""

import csv
import json
import random
import sqlite3
import sys
import time
from collections import Counter, defaultdict

TIE_TOLERANCE = 1e-9
TIERS = {"VERY_HIGH": 50, "HIGH": 40, "MEDIUM": 30}  # LOW = below MEDIUM
TIER_NAMES = ("VERY_HIGH", "HIGH", "MEDIUM", "LOW")
SR = {"S", "R"}

WEIGHT_CURVES = {
    "linear":  {1: 1.0, 2: 0.8, 3: 0.6, 4: 0.4, 5: 0.2},
    "steeper": {1: 1.0, 2: 0.5, 3: 0.25, 4: 0.125, 5: 0.0625},
    "flat":    {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0},
}


# ------------------------------------------------------------- load once --

def load_db(db_path):
    """[(bin_name, smorf_id, occurrences)]; each occurrence is
    (neighbours [(rank, cats), ...], target_cats)."""
    conn = sqlite3.connect(db_path)
    data = []
    for bin_name, smorf_id, occ_json in conn.execute(
        "SELECT bin_name, smorf_id, occurrences_json FROM smorfs ORDER BY bin_name, smorf_id"
    ):
        occs = [
            ([(rank, tuple(cats)) for rank, cats in o["neighbours"]], tuple(o["target_categories"]))
            for o in json.loads(occ_json)
        ]
        data.append((bin_name, smorf_id, occs))
    conn.close()
    return data


# --------------------------------------------------------------- methods --
# Every method: (occurrences, **params) -> (predicted, top_score_pct, known)

def resolve(votes):
    if not votes:
        return "NO_EVIDENCE", 0.0
    top = max(votes.values())
    winners = sorted(c for c, s in votes.items() if abs(s - top) < TIE_TOLERANCE)
    return (winners[0] if len(winners) == 1 else "UNRESOLVED"), round(top / sum(votes.values()) * 100, 1)


def known_categories(occs, exclude=()):
    return {c for _, target in occs for c in target if c not in exclude}


def rank_weighted_vote(occs, rank_weights, max_neighbours=5, exclude=()):
    """Each neighbour votes for its category, weighted by closeness,
    split evenly across multi-letter categories."""
    votes = defaultdict(float)
    for neighbours, _ in occs:
        for rank, cats in neighbours:
            if rank > max_neighbours or rank not in rank_weights:
                continue
            cats = [c for c in cats if c not in exclude]
            for c in cats:
                votes[c] += rank_weights[rank] / len(cats)
    return (*resolve(votes), known_categories(occs, exclude))


def majority_of_n(occs, n=5, min_agreement=4, exclude=()):
    """Take the n closest neighbours of every occurrence; predict a
    category only if at least min_agreement/n of them agree on it."""
    counts = defaultdict(int)
    total = 0
    for neighbours, _ in occs:
        for _, cats in sorted(neighbours, key=lambda rc: rc[0])[:n]:
            cats = [c for c in cats if c not in exclude]
            if cats:
                total += 1
                for c in cats:
                    counts[c] += 1
    known = known_categories(occs, exclude)
    if not counts:
        return "NO_EVIDENCE", 0.0, known
    top_cat, top = max(counts.items(), key=lambda kv: kv[1])
    if top / total >= min_agreement / n:
        return top_cat, round(top / total * 100, 1), known
    return "NO_EVIDENCE", 0.0, known


def control_always(occs, category):
    """CONTROL: ignore the neighbours, always predict one category."""
    return category, 0.0, known_categories(occs)


def control_random(occs, categories, weights, rng):
    """CONTROL: ignore the neighbours, pick a random category, drawn in
    proportion to how common each category is among annotated targets."""
    return rng.choices(categories, weights)[0], 0.0, known_categories(occs)


# ---------------------------------------------------------- combinations --

def build_combinations(data):
    """(label, method, params). Every loop multiplies the count."""
    combos = []
    for wname, weights in WEIGHT_CURVES.items():
        for n in (1, 2, 3, 4, 5):
            if n == 1 and wname != "linear":
                continue  # with 1 neighbour every weight curve gives the same answer
            for ex_name, ex in (("all-cats", ()), ("excl-SR", SR)):
                label = f"{wname}, {n} nbrs, {ex_name}"
                if wname == "linear" and n == 5:
                    label += " (= original)"
                combos.append((label, rank_weighted_vote,
                               dict(rank_weights=weights, max_neighbours=n, exclude=ex)))

    for n, k in ((3, 2), (5, 3), (5, 4), (7, 5)):
        for ex_name, ex in (("all-cats", ()), ("excl-SR", SR)):
            combos.append((f"majority {k} of {n}, {ex_name}", majority_of_n,
                           dict(n=n, min_agreement=k, exclude=ex)))

    # controls need the category frequencies of the annotated targets
    freq = Counter(c for _, _, occs in data for c in known_categories(occs))
    if freq:
        cats, weights = zip(*freq.most_common())
        combos.append((f"CONTROL: always '{cats[0]}' (most common)", control_always,
                       dict(category=cats[0])))
        combos.append(("CONTROL: random, by category frequency", control_random,
                       dict(categories=cats, weights=weights, rng=random.Random(0))))
    return combos


# -------------------------------------------------------------- evaluate --

def tier(pct):
    for name in TIER_NAMES[:-1]:
        if pct >= TIERS[name]:
            return name
    return "LOW"


def new_stats():
    return {"n": 0, "evaluated": 0, "correct": 0, **{t: [0, 0] for t in TIER_NAMES}}


def summary_row(label, bin_name, s):
    row = [label, bin_name, s["n"], f"{s['evaluated'] / s['n']:.3f}" if s["n"] else "",
           s["evaluated"], s["correct"],
           f"{s['correct'] / s['evaluated']:.3f}" if s["evaluated"] else ""]
    for t in TIER_NAMES:
        n, c = s[t]
        row += [n, f"{c / n:.3f}" if n else ""]
    return row


def main():
    db = sys.argv[1] if len(sys.argv) > 1 else "smorfs_all_bins.sqlite"

    t0 = time.time()
    data = load_db(db)
    bins = sorted({b for b, _, _ in data})
    print(f"Loaded {len(data):,} smORFs from {len(bins)} bins in {time.time() - t0:.1f}s", flush=True)

    combos = build_combinations(data)
    print(f"Running {len(combos)} combinations\n", flush=True)

    results = []                  # (label, overall accuracy, stats by bin)
    predictions = defaultdict(list)  # label -> predicted category per smORF, in data order
    for label, method, params in combos:
        t1 = time.time()
        stats = defaultdict(new_stats)
        for bin_name, _, occs in data:
            predicted, pct, known = method(occs, **params)
            predictions[label].append(predicted)
            for key in (bin_name, "ALL"):
                s = stats[key]
                s["n"] += 1
                if known and predicted not in ("UNRESOLVED", "NO_EVIDENCE"):
                    ok = predicted in known
                    s["evaluated"] += 1
                    s["correct"] += ok
                    s[tier(pct)][0] += 1
                    s[tier(pct)][1] += ok
        a = stats["ALL"]
        acc = a["correct"] / a["evaluated"] if a["evaluated"] else 0.0
        print(f"  {label:<44} coverage={a['evaluated'] / a['n']:.3f}  accuracy={acc:.3f}  "
              f"({time.time() - t1:.1f}s)", flush=True)
        results.append((label, acc, stats))

    results.sort(key=lambda r: r[1], reverse=True)
    with open("sweep_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["combination", "bin", "n_smorfs", "coverage", "evaluated", "correct", "accuracy"]
                   + [f"{t.lower()}_{x}" for t in TIER_NAMES for x in ("n", "acc")])
        for label, _, stats in results:
            w.writerow(summary_row(label, "ALL", stats["ALL"]))
            for b in bins:
                w.writerow(summary_row(label, b, stats[b]))

    labels = [label for label, _, _ in combos]
    with open("sweep_predictions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bin", "smorf_id", "known_category"] + labels)
        for i, (bin_name, smorf_id, occs) in enumerate(data):
            w.writerow([bin_name, smorf_id, ",".join(sorted(known_categories(occs)))]
                       + [predictions[label][i] for label in labels])

    print(f"\nTotal {time.time() - t0:.1f}s. Saved:")
    print("  sweep_summary.csv      (accuracy per combination and bin, best first)")
    print("  sweep_predictions.tsv  (every smORF's prediction from every combination)")


if __name__ == "__main__":
    main()
