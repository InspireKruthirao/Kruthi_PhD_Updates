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
3. Machine learning (Ben): random forest and logistic regression on one-hot
   COG letters of each neighbour position (5 per side), with strand,
   distance and the number of occurrences, trained on TRAIN and scored on
   TEST. Accuracy at each confidence level, so every prediction can be
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
N_POS = 5  # neighbour positions per side used as ML features

REPORT = None  # feedback_report.txt, written as it goes


def say(text=""):
    print(text, flush=True)
    if REPORT is not None:
        with open(REPORT, "a") as fh:
            fh.write(text + "\n")


def pct(a, b):
    return f"{100 * a / b:.1f}" if b else ""


# --------------------------------------------------------------- splits ----

def random_split(d):
    """TRAIN and TEST smORF indices, exactly as cascade.py"""
    annotated = [int(i) for i in np.flatnonzero(d.smorfs["annotated"].to_numpy())]
    random.Random(SEED).shuffle(annotated)
    half = len(annotated) // 2
    return np.array(sorted(annotated[:half])), np.array(sorted(annotated[half:]))


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

def known_matrix(d, exclude=()):
    """smORF x category: is it a known category (not excluded)?"""
    k = np.zeros((d.smorfs.height, d.n_cats), dtype=bool)
    t = d.targets.filter(~pl.col("cat").is_in(list(exclude)))
    code = t["cat"].cast(pl.Enum(d.labels[:d.n_cats])).to_physical().to_numpy()
    k[t["smorf"].to_numpy(), code] = True
    return k


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


# ------------------------------------------------------ 3 machine learning ----

def ml_features(d):
    """smORF x features: for each side and position 1-5, the share of
    occurrences whose gene there has each COG letter (multi-letter genes
    split), is on the smORF's strand, and its log distance; plus the log
    number of occurrences. Also the feature names."""
    nb = d.nb
    n_smorfs = d.smorfs.height
    n_occ = np.maximum(np.bincount(d.occ_smorf, minlength=n_smorfs), 1)
    use = nb["rank"] <= N_POS
    side = (d.info["side"][nb["ent"]].astype(np.int64) if d.info is not None
            else np.zeros(len(nb["ent"]), dtype=np.int64))
    slot = side * N_POS + nb["rank"].astype(np.int64) - 1  # 0..2*N_POS-1
    n_slot = 2 * N_POS

    letters = np.bincount((nb["smorf"].astype(np.int64) * n_slot + slot)[use] * d.n_cats + nb["code"][use],
                          weights=1 / nb["ncat"][use], minlength=n_smorfs * n_slot * d.n_cats)
    feats = [letters.reshape(n_smorfs, n_slot * d.n_cats)]
    names = [f"{'up' if s < N_POS else 'down'}{s % N_POS + 1}_{c}"
             for s in range(n_slot) for c in d.labels[:d.n_cats]]
    if d.info is not None:
        first = sw.run_starts(nb["ent"])
        first = first[use[first]]
        key = nb["smorf"][first].astype(np.int64) * n_slot + slot[first]
        size = n_smorfs * n_slot
        gene = nb["ent"][first]
        for name, value in (("same_strand", d.info["same"][gene].astype(float)),
                            ("log_distance", np.log1p(np.maximum(d.info["dist"][gene], 0))),
                            ("overlaps", (d.info["gap"][gene] < 0).astype(float))):
            feats.append(np.bincount(key, weights=value, minlength=size).reshape(n_smorfs, n_slot))
            names += [f"{'up' if s < N_POS else 'down'}{s % N_POS + 1}_{name}" for s in range(n_slot)]
    x = np.hstack(feats) / n_occ[:, None]
    x = np.hstack([x, np.log(n_occ)[:, None]])
    names.append("log_occurrences")
    return x.astype(np.float32), names


def train_rows(d, rows, k):
    """One training row per (smORF, known category), weight 1/categories"""
    r, c = np.nonzero(k[rows])
    n = k[rows].sum(axis=1)
    return rows[r], c, 1.0 / n[r]


def ml_test(d, x, has, train, test, ks, label, models, write=None):
    """Fit each model on TRAIN, report TEST accuracy by confidence"""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    say(f"\n  {label}: TRAIN {len(train):,}, TEST {len(test):,} annotated smORFs")
    for cname, ex in (("all categories", ()), ("excluding S and R", sw.SR)):
        k = ks[ex]
        tr = train[has[train] & k[train].any(axis=1)]
        rows, y, w = train_rows(d, tr, k)
        for mname in models:
            t = time.time()
            model = (RandomForestClassifier(n_estimators=300, min_samples_leaf=3, n_jobs=sw.cpus(),
                                            random_state=SEED) if mname == "random forest" else
                     make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000)))
            fit_kw = {"sample_weight": w} if mname == "random forest" else {"logisticregression__sample_weight": w}
            model.fit(x[rows], y, **fit_kw)
            classes = model.classes_
            prob = np.zeros((d.smorfs.height, d.n_cats))
            go = np.flatnonzero(has)
            prob[np.ix_(go, classes)] = model.predict_proba(x[go])
            pred = np.where(has, prob.argmax(axis=1), d.code("NO_EVIDENCE"))
            conf = prob.max(axis=1)
            say(f"\n  {mname}, {cname} ({time.time() - t:.0f}s)")
            say("  predict only if confidence >=\tcoverage %\taccuracy %")
            for step in CONFIDENCE_STEPS:
                p = np.where(conf >= step, pred, d.code("LOW_CONFIDENCE"))
                n, ok, cov, acc = score(d, p, test, k)
                say(f"  {step:.1f}\t{cov}\t{acc}")
            say("  confidence bin\tTEST smORFs\taccuracy %   (estimated accuracy per prediction)")
            pt, kt = pred[test], k[test]
            ev = (pt < d.n_cats) & kt.any(axis=1)
            ok = ev & kt[np.arange(len(test)), np.minimum(pt, d.n_cats - 1)]
            for lo, hi in zip(CALIBRATION_BINS[:-1], CALIBRATION_BINS[1:]):
                b = ev & (conf[test] >= lo) & (conf[test] < hi)
                say(f"  {lo:.1f}-{min(hi, 1):.1f}\t{b.sum():,}\t{pct((ok & b).sum(), b.sum())}")
            if mname == "random forest" and cname == "all categories":
                imp = np.argsort(-model.feature_importances_)[:15]
                say("  most important features: " + ", ".join(f"{NAMES[i]} ({model.feature_importances_[i]:.3f})"
                                                              for i in imp))
                if write is not None:
                    write(pred, conf)


def machine_learning(d, train, test, ks, clusters, folder):
    global NAMES
    say("\n3. MACHINE LEARNING on one-hot COG letters per neighbour position (TRAIN -> TEST)")
    x, NAMES = ml_features(d)
    has = x[:, :-1].any(axis=1)
    say(f"  {x.shape[1]} features; {has.sum():,} of {d.smorfs.height:,} smORFs have at least one")
    say("  (cascade.py on the same TEST half: 56.1% accuracy at 97.8% coverage, all categories)")

    def write(pred, conf):
        labels = pl.Series(d.labels)
        part = np.full(d.smorfs.height, "", dtype=object)
        part[train], part[test] = "TRAIN", "TEST"
        (d.smorfs.select("bin", "smorf_id", "known_category")
         .with_columns(split=pl.Series(part.tolist()), prediction=labels.gather(pred),
                       confidence=pl.Series(np.round(conf, 3)))
         .write_csv(folder / "ml_predictions.tsv", separator="\t"))

    ml_test(d, x, has, train, test, ks, "random split (as cascade.py)", ("random forest", "logistic regression"),
            write)
    if clusters:
        tr, te = cluster_split(d, clusters)
        ml_test(d, x, has, tr, te, ks, "split by 90% cluster", ("random forest",))


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
    train, test = random_split(d)
    say(f"annotated smORFs: TRAIN {len(train):,}, TEST {len(test):,} (random split, seed {SEED}, as cascade.py)")
    ks = {(): known_matrix(d), sw.SR: known_matrix(d, sw.SR), ("S",): known_matrix(d, ("S",)),
          ("R",): known_matrix(d, ("R",))}

    baselines(d, combos, test, ks)
    one_factor(combos, folder)
    exact_vs_approximate(d, combos, test, ks)
    machine_learning(d, train, test, ks, clusters, folder)

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Wrote feedback_report.txt, one_factor_comparisons.tsv, "
          f"ml_predictions.tsv in {folder}")


if __name__ == "__main__":
    main()
