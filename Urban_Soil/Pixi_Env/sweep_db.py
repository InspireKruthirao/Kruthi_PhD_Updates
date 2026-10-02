#!/usr/bin/env python3
"""
COG function prediction for smORFs: 100 rule combinations, one file.

PART 1 - PREDICT   every combination predicts a COG category for ALL smORFs
                   (annotated or not)            -> predictions_all_smorfs.tsv
PART 2 - EVALUATE  each combination is scored ONLY on smORFs that already
                   have a known COG category     -> combinations_results.csv
                                                    results_by_bin.csv

    python3 sweep_db.py [db]        (default db: smorfs_all_bins.sqlite)

All results are percentages. combinations_results.csv is the sheet: one
row per rule (C001...C100) with what it does and how it scored.
C027 (all categories) and C028 (excluding S and R) are the original method.

To add a rule: add it in build_combinations(). Nothing else changes.
"""

import csv
import json
import random
import sqlite3
import sys
import time
from array import array
from collections import Counter, defaultdict

SR = ("S", "R")
NO_CALL = {"NO_EVIDENCE", "UNRESOLVED", "LOW_CONFIDENCE", "TOO_FEW_OCCURRENCES"}
TIERS = (("VERY_HIGH", 50), ("HIGH", 40), ("MEDIUM", 30), ("LOW", 0))
TIE_TOLERANCE = 1e-9

MAX_RANK = 10
WEIGHTS = {  # weight of the neighbour at each distance (rank 1 = closest)
    "flat":       {r: 1.0 for r in range(1, MAX_RANK + 1)},
    "linear":     {r: max(0.0, 1.0 - 0.2 * (r - 1)) for r in range(1, MAX_RANK + 1)},  # 1, .8, .6, .4, .2
    "steeper":    {r: 0.5 ** (r - 1) for r in range(1, MAX_RANK + 1)},                # 1, .5, .25, ...
    "very_steep": {r: 0.25 ** (r - 1) for r in range(1, MAX_RANK + 1)},               # 1, .25, .06, ...
}


# =========================================================== load once ==

def load_db(db_path):
    """[(bin, smorf_id, occs)]; occs = [(neighbours [(rank, cats)], target_cats)]"""
    conn = sqlite3.connect(db_path)
    data = []
    for bin_name, smorf_id, occ_json in conn.execute(
        "SELECT bin_name, smorf_id, occurrences_json FROM smorfs ORDER BY bin_name, smorf_id"
    ):
        occs = [([(rank, tuple(cats)) for rank, cats in o["neighbours"]], tuple(o["target_categories"]))
                for o in json.loads(occ_json)]
        data.append((bin_name, smorf_id, occs))
    conn.close()
    return data


def known_categories(occs, exclude=()):
    return {c for _, target in occs for c in target if c not in exclude}


# ========================================================= the methods ==
# every method: (occs, **settings) -> (predicted_category, top_score_pct)

def resolve(votes):
    if not votes:
        return "NO_EVIDENCE", 0.0
    top = max(votes.values())
    winners = sorted(c for c, s in votes.items() if abs(s - top) < TIE_TOLERANCE)
    return (winners[0] if len(winners) == 1 else "UNRESOLVED"), round(top / sum(votes.values()) * 100, 1)


def tally(neighbours, weights, n, exclude, multi, background, votes):
    """Add one occurrence's neighbour votes into votes; True if it voted."""
    voted = False
    for rank, cats in neighbours:
        if rank > n:
            continue
        cats = [c for c in cats if c not in exclude]
        if not cats or (multi == "ignore" and len(cats) > 1):
            continue
        share = weights[rank] if multi == "full" else weights[rank] / len(cats)
        for c in cats:
            votes[c] += share / background[c] if background else share
        voted = True
    return voted


def weighted_vote(occs, curve="linear", n=5, exclude=(), multi="split",
                  background=None, cutoff=None, min_occ=None):
    """Neighbours vote for their category, closer neighbours count more.
    Options: drop categories, handle multi-letter categories, down-weight
    common categories, require a minimum winning share, or require a
    minimum number of occurrences that contribute votes."""
    votes = defaultdict(float)
    used = sum(tally(nb, WEIGHTS[curve], n, exclude, multi, background, votes) for nb, _ in occs)
    if min_occ and used < min_occ:
        return "TOO_FEW_OCCURRENCES", 0.0
    pred, pct = resolve(votes)
    if cutoff and pred not in NO_CALL and pct < cutoff:
        return "LOW_CONFIDENCE", pct
    return pred, pct


def majority(occs, n=5, k=4, exclude=()):
    """Take the n closest neighbours of every occurrence; predict only if
    at least k out of n of them agree on one category."""
    counts, total = defaultdict(int), 0
    for neighbours, _ in occs:
        for _, cats in sorted(neighbours, key=lambda rc: rc[0])[:n]:
            cats = [c for c in cats if c not in exclude]
            if cats:
                total += 1
                for c in cats:
                    counts[c] += 1
    if not counts:
        return "NO_EVIDENCE", 0.0
    top_cat, top = max(counts.items(), key=lambda kv: kv[1])
    if top / total >= k / n:
        return top_cat, round(top / total * 100, 1)
    return "NO_EVIDENCE", 0.0


def consensus(occs, curve="linear", n=5, exclude=()):
    """Each occurrence (genome) makes its own prediction; then the
    occurrences vote, one vote each."""
    winners = Counter()
    for neighbours, _ in occs:
        votes = defaultdict(float)
        tally(neighbours, WEIGHTS[curve], n, exclude, "split", None, votes)
        pred, _ = resolve(votes)
        if pred not in NO_CALL:
            winners[pred] += 1
    return resolve(winners)


def control_always(occs, category, exclude=()):
    """CONTROL: ignore the neighbours, always predict one category."""
    return category, 0.0


def control_random(occs, categories, freqs, rng, exclude=()):
    """CONTROL: ignore the neighbours, random category weighted by how
    common each category is among annotated smORFs."""
    return rng.choices(categories, freqs)[0], 0.0


def control_shuffled(occs, pool, rng, exclude=()):
    """CONTROL: run the original method on the neighbours of a DIFFERENT,
    randomly chosen smORF."""
    return weighted_vote(pool[rng.randrange(len(pool))])


# ================================================= the 100 combinations ==

def build_combinations(data):
    combos = []

    def add(group, rule, fn, **settings):
        combos.append(dict(id=f"C{len(combos) + 1:03d}", group=group, rule=rule, fn=fn,
                           settings=settings, exclude=settings.get("exclude", ())))

    cats_opts = (("all categories", ()), ("excluding S and R", SR))
    curve_name = {"flat": "flat", "linear": "linear", "steeper": "steeper", "very_steep": "very steep"}

    # A. rank-weighted vote: weighting x neighbours x categories (34)
    for n in (1, 2, 3, 4, 5):
        for curve in ("linear", "steeper", "flat", "very_steep"):
            if n == 1 and curve != "linear":
                continue  # with 1 neighbour every weighting gives the same answer
            for cname, ex in cats_opts:
                extra = " (ORIGINAL METHOD)" if curve == "linear" and n == 5 else ""
                add("A. Weighted vote", f"{curve_name[curve]} weights, {n} neighbour(s), {cname}{extra}",
                    weighted_vote, curve=curve, n=n, exclude=ex)

    # B. majority rule (10)
    for n, k in ((3, 2), (5, 3), (5, 4), (7, 5), (9, 6)):
        for cname, ex in cats_opts:
            add("B. Majority rule", f"at least {k} of the {n} nearest agree, {cname}",
                majority, n=n, k=k, exclude=ex)

    # C. confidence cutoff (20)
    for curve in ("linear", "steeper"):
        for cut in (30, 40, 50, 60, 70):
            for cname, ex in cats_opts:
                add("C. Confidence cutoff", f"{curve} weights, 5 neighbours, predict only if winner has "
                    f">= {cut}% of the vote, {cname}", weighted_vote, curve=curve, n=5, exclude=ex, cutoff=cut)

    # D. down-weight categories that are common everywhere (6)
    background = Counter()
    for _, _, occs in data:
        for neighbours, _ in occs:
            for rank, cats in neighbours:
                if rank <= 5:
                    for c in cats:
                        background[c] += 1 / len(cats)
    mean = sum(background.values()) / len(background) if background else 1
    background = {c: v / mean for c, v in background.items()}
    for curve, n in (("linear", 5), ("steeper", 5), ("steeper", 3)):
        for cname, ex in cats_opts:
            add("D. Down-weight common categories", f"{curve} weights, {n} neighbours, votes divided by "
                f"how common the category is overall, {cname}", weighted_vote,
                curve=curve, n=n, exclude=ex, background=background)

    # E. consensus across genomes (8)
    for curve in ("linear", "steeper"):
        for n in (3, 5):
            for cname, ex in cats_opts:
                add("E. Consensus across genomes", f"each genome predicts with {curve} weights and {n} "
                    f"neighbours, then genomes vote, {cname}", consensus, curve=curve, n=n, exclude=ex)

    # F. minimum evidence (6)
    for m in (2, 3, 5):
        for cname, ex in cats_opts:
            add("F. Minimum evidence", f"steeper weights, 5 neighbours, predict only if >= {m} genomes "
                f"contribute votes, {cname}", weighted_vote, curve="steeper", n=5, exclude=ex, min_occ=m)

    # G. more neighbours (4)
    for n in (7, 10):
        for cname, ex in cats_opts:
            add("G. More neighbours", f"steeper weights, {n} neighbours, {cname}",
                weighted_vote, curve="steeper", n=n, exclude=ex)

    # H. multi-letter categories such as 'EG' (4)
    for multi, text in (("full", "count every letter fully"), ("ignore", "ignore multi-letter genes")):
        for cname, ex in cats_opts:
            add("H. Multi-letter categories", f"steeper weights, 5 neighbours, {text}, {cname}",
                weighted_vote, curve="steeper", n=5, exclude=ex, multi=multi)

    # I. other exclusions (4)
    for curve in ("linear", "steeper"):
        for letter in ("S", "R"):
            add("I. Other exclusions", f"{curve} weights, 5 neighbours, excluding {letter} only",
                weighted_vote, curve=curve, n=5, exclude=(letter,))

    # J. controls (4)
    freq = Counter(c for _, _, occs in data for c in known_categories(occs))
    freq_no_sr = Counter({c: v for c, v in freq.items() if c not in SR})
    if freq:
        top = freq.most_common(1)[0][0]
        add("J. Control", f"always predict '{top}' (most common category), all categories",
            control_always, category=top)
    if freq_no_sr:
        top = freq_no_sr.most_common(1)[0][0]
        add("J. Control", f"always predict '{top}' (most common category), excluding S and R",
            control_always, category=top, exclude=SR)
    if freq:
        cats, counts = zip(*freq.most_common())
        add("J. Control", "random category, weighted by how common each is, all categories",
            control_random, categories=cats, freqs=counts, rng=random.Random(0))
    add("J. Control", "original method run on the neighbours of a different random smORF",
        control_shuffled, pool=[occs for _, _, occs in data], rng=random.Random(1))
    return combos


# ===================================================== PART 1: PREDICT ==

def predict(data, combos):
    """Every combination predicts for every smORF, annotated or not."""
    preds, pcts = {}, {}
    for c in combos:
        t = time.time()
        p, s = [], array("d")
        for _, _, occs in data:
            pred, pct = c["fn"](occs, **c["settings"])
            p.append(pred)
            s.append(pct)
        preds[c["id"]], pcts[c["id"]] = p, s
        made = sum(x not in NO_CALL for x in p)
        print(f"  {c['id']}  predicted {made / len(data) * 100:5.1f}% of smORFs  "
              f"({time.time() - t:.1f}s)  {c['rule']}", flush=True)

    with open("predictions_all_smorfs.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bin", "smorf_id", "annotated", "known_category"] + [c["id"] for c in combos])
        for i, (b, sid, occs) in enumerate(data):
            known = ",".join(sorted(known_categories(occs)))
            w.writerow([b, sid, "yes" if known else "no", known] + [preds[c["id"]][i] for c in combos])
    return preds, pcts


# ==================================================== PART 2: EVALUATE ==

def pct(a, b):
    return f"{a / b * 100:.1f}" if b else ""


def tier_of(score):
    return next(name for name, low in TIERS if score >= low)


def evaluate(data, combos, preds, pcts):
    """Score each combination ONLY on smORFs with a known COG category."""
    annotated = [bool(known_categories(occs)) for _, _, occs in data]
    bins = sorted({b for b, _, _ in data})
    rows, by_bin = [], []

    for c in combos:
        p, s = preds[c["id"]], pcts[c["id"]]
        st = defaultdict(lambda: {"n": 0, "pred": 0, "annot": 0, "eval": 0, "correct": 0,
                                  **{t: [0, 0] for t, _ in TIERS}})
        for i, (b, _, occs) in enumerate(data):
            made = p[i] not in NO_CALL
            known = known_categories(occs, c["exclude"]) if annotated[i] else set()
            for key in (b, "ALL"):
                x = st[key]
                x["n"] += 1
                x["pred"] += made
                x["annot"] += annotated[i]
                if made and known:
                    ok = p[i] in known
                    x["eval"] += 1
                    x["correct"] += ok
                    x[tier_of(s[i])][0] += 1
                    x[tier_of(s[i])][1] += ok
        a = st["ALL"]
        row = {"ID": c["id"], "Group": c["group"], "Rule": c["rule"],
               "smORFs predicted (% of all)": pct(a["pred"], a["n"]),
               "Coverage (% of annotated smORFs evaluated)": pct(a["eval"], a["annot"]),
               "Evaluated (n)": a["eval"], "Correct (n)": a["correct"],
               "Accuracy (%)": pct(a["correct"], a["eval"])}
        for t, _ in TIERS:
            row[f"{t} n"] = a[t][0]
            row[f"{t} accuracy (%)"] = pct(a[t][1], a[t][0])
        rows.append(row)
        for b in bins:
            x = st[b]
            by_bin.append({"ID": c["id"], "Rule": c["rule"], "Bin": b,
                           "Coverage (%)": pct(x["eval"], x["annot"]), "Evaluated (n)": x["eval"],
                           "Correct (n)": x["correct"], "Accuracy (%)": pct(x["correct"], x["eval"])})

    ranked = sorted(rows, key=lambda r: float(r["Accuracy (%)"] or 0), reverse=True)
    for rank, r in enumerate(ranked, 1):
        r["Accuracy rank"] = rank
    cols = ["ID", "Accuracy rank", "Group", "Rule", "Accuracy (%)",
            "Coverage (% of annotated smORFs evaluated)", "smORFs predicted (% of all)",
            "Evaluated (n)", "Correct (n)"] + [f"{t} {x}" for t, _ in TIERS for x in ("n", "accuracy (%)")]
    for name, out, fields in (("combinations_results.csv", rows, cols),
                              ("results_by_bin.csv", by_bin, list(by_bin[0]))):
        with open(name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(out)

    print(f"\nAnnotated smORFs used for evaluation: {sum(annotated):,} of {len(data):,}")
    print("\nTop 10 by accuracy:")
    for r in ranked[:10]:
        print(f"  {r['ID']}  accuracy {r['Accuracy (%)']:>5}%  coverage "
              f"{r['Coverage (% of annotated smORFs evaluated)']:>5}%  {r['Rule']}")
    orig = next(r for r in rows if "ORIGINAL" in r["Rule"] and "all categories" in r["Rule"])
    print(f"\nOriginal method ({orig['ID']}): accuracy {orig['Accuracy (%)']}%  "
          f"coverage {orig['Coverage (% of annotated smORFs evaluated)']}%")


def main():
    db = sys.argv[1] if len(sys.argv) > 1 else "smorfs_all_bins.sqlite"
    t0 = time.time()
    data = load_db(db)
    print(f"Loaded {len(data):,} smORFs in {time.time() - t0:.1f}s", flush=True)
    combos = build_combinations(data)

    print(f"\nPART 1 - PREDICT: {len(combos)} combinations on all smORFs", flush=True)
    preds, pcts = predict(data, combos)

    print("\nPART 2 - EVALUATE: annotated smORFs only", flush=True)
    evaluate(data, combos, preds, pcts)

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Saved:")
    print("  combinations_results.csv   the sheet: one row per rule, accuracy in %")
    print("  results_by_bin.csv         the same, per bin")
    print("  predictions_all_smorfs.tsv every smORF x every rule (C001...)")


if __name__ == "__main__":
    main()
