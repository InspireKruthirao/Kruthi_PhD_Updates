#!/usr/bin/env python3
"""
Tests suggested at the supervisor meeting of 28 October 2026 (Luis, Ben),
on the same TRAIN/TEST split of the annotated smORFs as cascade.py.

    python3 feedback_tests.py smorf_neighbourhoods.sqlite RESULTS_FOLDER [--clusters FILE.clstr]

RESULTS_FOLDER is a sweep_db.py --save output folder (combinations_results.csv
is read from it); everything is written there:

    feedback_report.txt         all the tables below
    one_factor_comparisons.tsv  every pair of rules that differ in one setting
    ml_predictions.tsv          random forest prediction + confidence per smORF

1. Baselines on TEST: always the most common category (Luis: the baseline
   to beat, rather than a random category).
2. Exact against approximate nearest neighbours (Ben/Luis: 150k smORFs is
   small enough to search exactly): how many of the exact neighbours Annoy
   and Voyager find, and the TEST accuracy of each.
3. Random forest (Ben; the AE rules of sweep_db.py) on one-hot COG letters
   of each neighbour position (5 per side), with strand, distance and the
   number of occurrences, trained on TRAIN and scored on TEST. Accuracy at each confidence level, so every prediction can be
   reported with an estimated accuracy (Luis). With --clusters (CD-HIT
   .clstr of the smORFs at 90%), also with TRAIN/TEST split by 90% cluster,
   so near-identical smORFs are never on both sides.
4. One factor at a time (Ben): every pair of rules that differ in exactly
   one setting, with the change in accuracy and coverage, summarised per
   setting. Uses the accuracies on all annotated smORFs from the sweep.

A prediction is correct if it is any of the smORF's known categories (the
excluded ones not counted, for rules "excluding S and R"). Coverage is the
share of the TEST smORFs that get a prediction and are evaluated.
"""

import csv
import inspect
import random
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl

import sweep_db as sw

SEED = 0  # as cascade.py
CONFIDENCE_STEPS = (0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
CALIBRATION_BINS = (0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01)

REPORT = None  # feedback_report.txt, written as it goes


def say(text=""):
    print(text, flush=True)
    if REPORT is not None:
        with open(REPORT, "a") as fh:
            fh.write(text + "\n")


def pct(a, b):
    return f"{100 * a / b:.1f}" if b else ""


# --------------------------------------------------------------- splits ----

def read_clusters(path, wanted):
    """smorf_id -> 90% cluster number, for the smorf_ids in wanted"""
    found, cluster = {}, -1
    pattern = re.compile(r">(\S+?)\.\.\.")
    with open(path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                cluster += 1
                continue
            m = pattern.search(line)
            if m and m.group(1) in wanted:
                found[m.group(1)] = cluster
    return found


def cluster_split(d, path):
    """TRAIN and TEST, whole 90% clusters at a time (about half each)"""
    annotated = np.flatnonzero(d.smorfs["annotated"].to_numpy())
    ids = d.smorfs["smorf_id"].to_list()
    found = read_clusters(path, {ids[i] for i in annotated})
    cl = np.array([found.get(ids[i], -1 - i) for i in annotated])  # not found: own cluster
    say(f"  90% clusters: {len(found):,} of {len(annotated):,} annotated smORFs found, "
        f"{len(np.unique(cl)):,} clusters")
    groups = sorted(set(cl.tolist()))
    random.Random(SEED).shuffle(groups)
    size = dict(zip(*np.unique(cl, return_counts=True)))
    train_groups, n = set(), 0
    for g in groups:
        if n >= len(annotated) // 2:
            break
        train_groups.add(g)
        n += size[g]
    in_train = np.array([g in train_groups for g in cl])
    return annotated[in_train], annotated[~in_train]


# -------------------------------------------------------------- scoring ----

def score(d, pred, rows, k):
    """(evaluated, correct, coverage %, accuracy %) of pred on rows"""
    p = pred[rows]
    ev = (p < d.n_cats) & k[rows].any(axis=1)
    ok = ev & k[rows, np.minimum(p, d.n_cats - 1)]
    return int(ev.sum()), int(ok.sum()), pct(ev.sum(), len(rows)), pct(ok.sum(), ev.sum())


def run_rule(d, combo):
    return combo["fn"](d, **combo["settings"])[0]


# ---------------------------------------------------------- 1 baselines ----

def baselines(d, combos, test, ks):
    say("\n1. BASELINES on TEST (always the most common category)")
    say("rule\tcoverage %\taccuracy %\tdescription")
    for c in combos:
        if c["fn"] is sw.control_always:
            n, ok, cov, acc = score(d, run_rule(d, c), test, ks[tuple(c["exclude"])])
            say(f"{c['id']}\t{cov}\t{acc}\t{c['rule']}")


# ------------------------------------------------- 2 exact vs approximate ----

def exact_vs_approximate(d, combos, test, ks):
    say("\n2. EXACT AGAINST APPROXIMATE NEAREST NEIGHBOURS (COG IDs within 3 genes, 512 dims)")
    for k in (10, 25):
        t = time.time()
        exact = sw.nearest(d, "exact", k, 3, 512)
        if exact is None:
            say("  no gene information in this database: skipped")
            return
        who, nn_x, sim_x = exact
        say(f"  exact search, k = {k}: {time.time() - t:.0f}s for {len(who):,} smORFs")
        # the k-th most similar other smORF, exactly (ties count as found)
        other = (nn_x >= 0) & (nn_x != who[:, None])
        kth = np.array([s[o][:k][-1] if o.any() else 1.0 for s, o in zip(sim_x, other)])
        for engine, opts in (("annoy", {"n_trees": 20}), ("voyager", {"M": 12, "ef_construction": 200,
                                                                     "query_ef": 100})):
            t = time.time()
            found = sw.nearest(d, engine, k, 3, 512, **opts)
            if found is None:
                continue
            w2, nn_a, sim_a = found
            assert np.array_equal(w2, who)
            o = (nn_a >= 0) & (nn_a != who[:, None])
            o &= np.cumsum(o, axis=1) <= k
            same = np.array([len(set(a[oa]) & set(x[ox][:k])) for a, oa, x, ox in zip(nn_a, o, nn_x, other)])
            good = ((sim_a >= kth[:, None] - 1e-5) & o).sum(axis=1)
            want = np.minimum(other.sum(axis=1), k)
            say(f"  {engine:8s} k = {k}: finds {pct(same.sum(), want.sum())}% of the exact neighbours; "
                f"{pct(good.sum(), want.sum())}% of its neighbours are as similar as the exact k-th "
                f"(ties count) ({time.time() - t:.0f}s)")
    say("\n  TEST accuracy of the similar-neighbourhood rules")
    say("rule\tcoverage %\taccuracy %\tdescription")
    for c in combos:
        if c["fn"] is sw.similar_vote:
            n, ok, cov, acc = score(d, run_rule(d, c), test, ks[tuple(c["exclude"])])
            say(f"{c['id']}\t{cov}\t{acc}\t{c['rule']}")


# ------------------------------------------------------ 3 random forest ----

def forest_table(d, prob, test, k):
    """TEST accuracy at each confidence level, and per confidence bin"""
    conf = prob.max(axis=1)
    pred = np.where(conf > 0, prob.argmax(axis=1), d.code("NO_EVIDENCE"))
    say("  predict only if confidence >=\tcoverage %\taccuracy %")
    for step in CONFIDENCE_STEPS:
        n, ok, cov, acc = score(d, np.where(conf >= step, pred, d.code("LOW_CONFIDENCE")), test, k)
        say(f"  {step:.1f}\t{cov}\t{acc}")
    say("  confidence bin\tTEST smORFs\taccuracy %   (estimated accuracy per prediction)")
    pt, kt = pred[test], k[test]
    ev = (pt < d.n_cats) & kt.any(axis=1)
    ok = ev & kt[np.arange(len(test)), np.minimum(pt, d.n_cats - 1)]
    for lo, hi in zip(CALIBRATION_BINS[:-1], CALIBRATION_BINS[1:]):
        b = ev & (conf[test] >= lo) & (conf[test] < hi)
        say(f"  {lo:.1f}-{min(hi, 1):.1f}\t{b.sum():,}\t{pct((ok & b).sum(), b.sum())}")
    return pred, conf


def machine_learning(d, train, test, ks, clusters, folder):
    """The random forest rules of sweep_db.py (random_forest): accuracy by
    confidence on TEST, the most important features, and the same with
    TRAIN/TEST split by 90% cluster"""
    say("\n3. RANDOM FOREST (sweep_db.py rules AE): TRAIN -> TEST, by confidence")
    x = sw.ml_features(d)
    names = sw.feature_names(d)
    say(f"  {x.shape[1]} features; {x[:, :-1].any(axis=1).sum():,} of {d.smorfs.height:,} smORFs have at least one")
    cascade = folder / "cascade_results.tsv"
    if cascade.exists():  # the cascade's last step, for comparison
        steps = [l.split("\t") for l in cascade.read_text().splitlines() if l[:1].isdigit()]
        if steps:
            say(f"  (cascade.py on the same TEST half: {steps[-1][5]}% accuracy at {steps[-1][4]}% coverage)")

    splits = [("random split (as cascade.py)", train, test)]
    if clusters:
        splits.append(("split by 90% cluster", *cluster_split(d, clusters)))
    for label, tr, te in splits:
        say(f"\n  {label}: TRAIN {len(tr):,}, TEST {len(te):,} annotated smORFs")
        for cname, ex in (("all categories", ()), ("excluding S and R", sw.SR)):
            t = time.time()
            model = sw.fit_forest(x, ks[ex], tr)
            prob = sw.forest_probabilities(model, x, d.n_cats)
            say(f"\n  random forest, {cname} ({time.time() - t:.0f}s)")
            pred, conf = forest_table(d, prob, te, ks[ex])
            if cname == "all categories":
                imp = np.argsort(-model.feature_importances_)[:15]
                say("  most important features: " + ", ".join(
                    f"{names[i]} ({model.feature_importances_[i]:.3f})" for i in imp))
                if te is test:
                    labels = pl.Series(d.labels)
                    part = np.full(d.smorfs.height, "", dtype=object)
                    part[train], part[test] = "TRAIN", "TEST"
                    (d.smorfs.select("bin", "smorf_id", "known_category")
                     .with_columns(split=pl.Series(part.tolist()), prediction=labels.gather(pred),
                                   confidence=pl.Series(np.round(conf, 3)))
                     .write_csv(folder / "ml_predictions.tsv", separator="\t"))


# -------------------------------------------------- 4 one factor at a time ----

def full_settings(c):
    """The rule's settings with the function's defaults filled in"""
    s = {k: p.default for k, p in inspect.signature(c["fn"]).parameters.items()
         if p.default is not inspect.Parameter.empty}
    s.update(c["settings"])
    return s


def same_value(a, b):
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray) or isinstance(a, random.Random):
        return a is b
    try:
        return bool(a == b)
    except Exception:
        return a is b


def shown(v):
    if isinstance(v, np.ndarray):
        return "divided by how common"
    if isinstance(v, tuple) and not v:
        return "none"
    if isinstance(v, tuple) and all(isinstance(x, str) for x in v):
        return ",".join(v)
    return str(v)


def one_factor(combos, folder):
    say("\n4. ONE FACTOR AT A TIME: rule pairs that differ in exactly one setting")
    say("   (accuracy and coverage on all annotated smORFs, from combinations_results.csv)")
    with open(folder / "combinations_results.csv", newline="") as fh:
        res = {r["ID"]: r for r in csv.DictReader(fh)}
    cov_col = "Coverage (% of annotated smORFs evaluated)"
    pairs = []
    by_fn = defaultdict(list)
    for c in combos:
        if c["id"] in res and res[c["id"]]["Accuracy (%)"] and not c["group"].startswith("J."):
            by_fn[c["fn"].__name__].append(c)
    for fn, cs in by_fn.items():
        full = [full_settings(c) for c in cs]
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                keys = set(full[i]) | set(full[j])
                diff = [k for k in keys if not same_value(full[i].get(k), full[j].get(k))]
                if len(diff) != 1:
                    continue
                f = diff[0]
                a, b = (i, j)
                va, vb = full[a].get(f), full[b].get(f)
                if isinstance(va, (int, float)) and isinstance(vb, (int, float)) and va > vb:
                    a, b, va, vb = b, a, vb, va
                elif va is None and vb is not None:
                    pass
                elif vb is None and va is not None:
                    a, b, va, vb = b, a, vb, va
                ra, rb = res[cs[a]["id"]], res[cs[b]["id"]]
                pairs.append({"function": fn, "setting": f, "from": shown(va), "to": shown(vb),
                              "rule A": cs[a]["id"], "rule B": cs[b]["id"],
                              "accuracy A": float(ra["Accuracy (%)"]), "accuracy B": float(rb["Accuracy (%)"]),
                              "coverage A": float(ra[cov_col]), "coverage B": float(rb[cov_col])})
    for p in pairs:
        p["accuracy change"] = round(p["accuracy B"] - p["accuracy A"], 1)
        p["coverage change"] = round(p["coverage B"] - p["coverage A"], 1)
    if not pairs:
        say("  no pairs found")
        return
    with open(folder / "one_factor_comparisons.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    say(f"  {len(pairs)} pairs (all in one_factor_comparisons.tsv). Summary, median change from A to B:")
    say("method\tsetting\tfrom -> to\tpairs\taccuracy change\tcoverage change")
    groups = defaultdict(list)
    for p in pairs:
        groups[(p["function"], p["setting"], p["from"], p["to"])].append(p)
    for (fn, f, va, vb), ps in sorted(groups.items()):
        say(f"{fn}\t{f}\t{va} -> {vb}\t{len(ps)}\t"
            f"{np.median([p['accuracy change'] for p in ps]):+.1f}\t{np.median([p['coverage change'] for p in ps]):+.1f}")
    say("  (setting 'exclude': the two rules are scored on different smORFs, so compare with care)")


# ------------------------------------------------------------------ main ----

def main():
    args = sys.argv[1:]
    clusters = None
    if "--clusters" in args:
        i = args.index("--clusters")
        clusters = args[i + 1]
        del args[i:i + 2]
    if len(args) != 2:
        sys.exit(__doc__)
    db, folder = args[0], Path(args[1])
    global REPORT
    REPORT = folder / "feedback_report.txt"
    REPORT.write_text("")
    t0 = time.time()
    d = sw.load_db(db)
    say(f"Loaded {d.smorfs.height:,} smORFs in {time.time() - t0:.0f}s")
    combos = sw.build_combinations(d)
    train, test = sw.halves(d)
    say(f"annotated smORFs: TRAIN {len(train):,}, TEST {len(test):,} (random split, seed {SEED}, as cascade.py)")
    ks = {ex: sw.known_matrix(d, ex) for ex in ((), sw.SR, ("S",), ("R",))}

    baselines(d, combos, test, ks)
    one_factor(combos, folder)
    exact_vs_approximate(d, combos, test, ks)
    machine_learning(d, train, test, ks, clusters, folder)

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Wrote feedback_report.txt, one_factor_comparisons.tsv, "
          f"ml_predictions.tsv in {folder}")


if __name__ == "__main__":
    main()
