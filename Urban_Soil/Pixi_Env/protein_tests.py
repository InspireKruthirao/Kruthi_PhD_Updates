#!/usr/bin/env python3
"""
Can a model trained on NORMAL PROTEINS predict smORF functions? (Then no
annotated smORFs are needed.) And how well do the rules do on proteins?

    python3 protein_tests.py protein_neighbourhoods.sqlite smorf_neighbourhoods.sqlite OUT_FOLDER \
        [--protein-rules combinations_results.csv] [--smorf-rules combinations_results.csv]

Random forest on the same features as feedback_tests.py (one-hot COG letters
of 5 genes per side, strand, distance, overlap; without the number of
occurrences, as every protein has one), each scored by how often the
predicted category is one of the known ones:

  A. proteins -> proteins        trained on half the proteins, tested on the other half
  B. proteins -> smORFs          trained on all proteins, tested on the TEST half of the
                                 annotated smORFs (as cascade.py) and on all of them.
                                 Two ways: (1) features averaged over a smORF's
                                 occurrences; (2) each occurrence predicted on its
                                 own (like a protein) and the probabilities averaged
  C. smORFs -> smORFs            trained on the TRAIN half (the reference)
  D. proteins + smORFs -> smORFs trained on all proteins and the TRAIN half

With --protein-rules and --smorf-rules (sweep_db.py results on each), also the
rules side by side: rules_proteins_vs_smorfs.tsv.

Writes protein_tests_report.txt (as it goes) in OUT_FOLDER.
"""

import csv
import random
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

import sweep_db as sw

SEED = 0
N_POS = 5
STEPS = (0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
REPORT = None


def say(text=""):
    print(text, flush=True)
    if REPORT is not None:
        with open(REPORT, "a") as fh:
            fh.write(text + "\n")


def pct(a, b):
    return f"{100 * a / b:.1f}" if b else ""


def features(d, cats, by="smorf"):
    """Rows = smORFs (by="smorf", averaged over occurrences) or occurrences
    (by="occ"). Columns: per side and position 1-5, the COG letters (in the
    order of cats; multi-letter genes split), same strand, log distance and
    overlap."""
    nb = d.nb
    to_cat = np.array([cats.index(c) for c in d.labels[:d.n_cats]])
    nc, n_slot = len(cats), 2 * N_POS
    group = nb[by].astype(np.int64)
    n_rows = d.smorfs.height if by == "smorf" else len(d.occ_smorf)
    use = nb["rank"] <= N_POS
    side = d.info["side"][nb["ent"]].astype(np.int64)
    slot = side * N_POS + nb["rank"].astype(np.int64) - 1
    key = group * n_slot + slot
    letters = np.bincount((key[use]) * nc + to_cat[nb["code"][use]], weights=1 / nb["ncat"][use],
                          minlength=n_rows * n_slot * nc).reshape(n_rows, n_slot * nc)
    feats = [letters]
    first = sw.run_starts(nb["ent"])
    first = first[use[first]]
    gene = nb["ent"][first]
    for value in (d.info["same"][gene].astype(float), np.log1p(np.maximum(d.info["dist"][gene], 0)),
                  (d.info["gap"][gene] < 0).astype(float)):
        feats.append(np.bincount(key[first], weights=value, minlength=n_rows * n_slot).reshape(n_rows, n_slot))
    x = np.hstack(feats)
    if by == "smorf":
        x /= np.maximum(np.bincount(d.occ_smorf, minlength=n_rows), 1)[:, None]
    return x.astype(np.float32)


def known(d, cats, exclude=()):
    """smORF x cats: is it a known category (not excluded)?"""
    k = np.zeros((d.smorfs.height, len(cats)), dtype=bool)
    t = d.targets.filter(~pl.col("cat").is_in(list(exclude)))
    k[t["smorf"].to_numpy(), [cats.index(c) for c in t["cat"].to_list()]] = True
    return k


def train(x, k, rows):
    """Random forest on rows (one training row per known category)"""
    from sklearn.ensemble import RandomForestClassifier
    rows = rows[k[rows].any(axis=1) & x[rows].any(axis=1)]
    r, c = np.nonzero(k[rows])
    w = 1.0 / k[rows].sum(axis=1)[r]
    model = RandomForestClassifier(n_estimators=300, min_samples_leaf=3, n_jobs=sw.cpus(), random_state=SEED)
    model.fit(x[rows[r]], c, sample_weight=w)
    return model


def probabilities(model, x, n_cats):
    """rows x categories (0 where the row has no features)"""
    p = np.zeros((len(x), n_cats))
    go = np.flatnonzero(x.any(axis=1))
    for b in range(0, len(go), 200_000):
        g = go[b:b + 200_000]
        p[np.ix_(g, model.classes_)] = model.predict_proba(x[g])
    return p


def report(label, prob, k, rows):
    """Accuracy at each confidence level, on rows"""
    conf = prob[rows].max(axis=1)
    pred = prob[rows].argmax(axis=1)
    has = conf > 0
    kr = k[rows]
    ev_all = has & kr.any(axis=1)
    ok_all = ev_all & kr[np.arange(len(rows)), pred]
    say(f"\n  {label}  ({len(rows):,} smORFs/proteins)")
    say("  confidence >=\tcoverage %\taccuracy %")
    for s in STEPS:
        ev = ev_all & (conf >= s)
        say(f"  {s:.1f}\t{pct(ev.sum(), len(rows))}\t{pct((ok_all & ev).sum(), ev.sum())}")
    return pct(ok_all.sum(), ev_all.sum())


def occurrence_average(prob_occ, d):
    """Average the occurrence probabilities of each smORF"""
    n = d.smorfs.height
    has = prob_occ.any(axis=1)
    out = np.zeros((n, prob_occ.shape[1]))
    np.add.at(out, d.occ_smorf[has], prob_occ[has])
    return out / np.maximum(np.bincount(d.occ_smorf[has], minlength=n), 1)[:, None]


def smorf_split(s):
    """TRAIN and TEST halves of the annotated smORFs, exactly as cascade.py"""
    annotated = [int(i) for i in np.flatnonzero(s.smorfs["annotated"].to_numpy())]
    random.Random(SEED).shuffle(annotated)
    half = len(annotated) // 2
    return np.array(sorted(annotated[:half])), np.array(sorted(annotated[half:]))


def rules_side_by_side(protein_csv, smorf_csv, folder):
    say("\nRULES ON PROTEINS AND ON smORFs (accuracy % / coverage %; sweep_db.py on each database)")
    cov = "Coverage (% of annotated smORFs evaluated)"
    read = lambda p: {r["ID"]: r for r in csv.DictReader(open(p, newline=""))}
    pr, sm = read(protein_csv), read(smorf_csv)
    rows = []
    for i, r in sm.items():
        if i in pr:
            rows.append({"ID": i, "Group": r["Group"], "Rule": r["Rule"],
                         "smORF accuracy %": r["Accuracy (%)"], "smORF coverage %": r[cov],
                         "protein accuracy %": pr[i]["Accuracy (%)"], "protein coverage %": pr[i][cov]})
    with open(folder / "rules_proteins_vs_smorfs.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    a = lambda r, x: float(r[x] or 0)
    with_both = [r for r in rows if a(r, "smORF coverage %") >= 5 and a(r, "protein coverage %") >= 5
                 and not r["Group"].startswith("J.")]
    if len(with_both) > 2:
        corr = np.corrcoef([a(r, "smORF accuracy %") for r in with_both],
                           [a(r, "protein accuracy %") for r in with_both])[0, 1]
        say(f"  {len(with_both)} rules with >= 5% coverage on both; correlation of their accuracies: {corr:.2f}")
    say("  rule\tsmORFs\tproteins\tdescription")
    for i in ("C172", "C166", "C128", "C027", "C037", "C229", "C233", "C097", "C098"):
        if i in sm and i in pr:
            say(f"  {i}\t{sm[i]['Accuracy (%)']} / {sm[i][cov]}\t{pr[i]['Accuracy (%)']} / {pr[i][cov]}\t"
                f"{sm[i]['Rule'][:100]}")
    say("  (rules that need >= 2 occurrences predict nothing for proteins: each protein is one occurrence)")


def main():
    global REPORT
    args = sys.argv[1:]
    opts = {}
    for name in ("--protein-rules", "--smorf-rules"):
        if name in args:
            i = args.index(name)
            opts[name] = args[i + 1]
            del args[i:i + 2]
    if len(args) != 3:
        sys.exit(__doc__)
    folder = Path(args[2])
    folder.mkdir(parents=True, exist_ok=True)
    REPORT = folder / "protein_tests_report.txt"
    REPORT.write_text("")
    t0 = time.time()

    p = sw.load_db(args[0])
    s = sw.load_db(args[1])
    say(f"Loaded {p.smorfs.height:,} proteins and {s.smorfs.height:,} smORFs ({time.time() - t0:.0f}s)")
    if p.info is None or s.info is None:
        sys.exit("both databases need gene strands and distances")
    cats = sorted(set(p.labels[:p.n_cats]) | set(s.labels[:s.n_cats]))
    xp = features(p, cats)
    xs = features(s, cats)
    xs_occ = features(s, cats, by="occ")
    s_train, s_test = smorf_split(s)
    annotated = np.flatnonzero(s.smorfs["annotated"].to_numpy())
    say(f"{xp.shape[1]} features; annotated smORFs: {len(annotated):,} (TRAIN {len(s_train):,}, "
        f"TEST {len(s_test):,})")
    counts = p.targets["cat"].value_counts().sort("count", descending=True)
    say("most common protein categories: " + ", ".join(f"{c} {n:,}" for c, n in counts.head(6).iter_rows()))

    summary = []
    for cname, ex in (("all categories", ()), ("excluding S and R", sw.SR)):
        say(f"\n==================== {cname} ====================")
        kp, ks = known(p, cats, ex), known(s, cats, ex)

        say("\nA. proteins -> proteins (half / half)")
        rows = np.flatnonzero(kp.any(axis=1))
        rng = np.random.default_rng(SEED)
        half = rng.permutation(rows)
        tr, te = np.sort(half[:len(half) // 2]), np.sort(half[len(half) // 2:])
        t = time.time()
        m = train(xp, kp, tr)
        acc = report(f"trained on {len(tr):,} proteins, tested on the other {len(te):,} ({time.time() - t:.0f}s)",
                     probabilities(m, xp, len(cats)), kp, te)
        summary.append((cname, "A. proteins -> proteins", acc))

        say("\nB. proteins -> smORFs (no annotated smORFs used for training)")
        t = time.time()
        m = train(xp, kp, rows)
        say(f"  trained on {len(rows):,} proteins ({time.time() - t:.0f}s)")
        prob_avg = probabilities(m, xs, len(cats))
        prob_occ = occurrence_average(probabilities(m, xs_occ, len(cats)), s)
        acc1 = report("(1) features averaged over occurrences, TEST half", prob_avg, ks, s_test)
        acc2 = report("(2) each occurrence predicted, probabilities averaged, TEST half", prob_occ, ks, s_test)
        report("(2) each occurrence predicted, probabilities averaged, ALL annotated smORFs", prob_occ, ks, annotated)
        summary.append((cname, "B1. proteins -> smORFs (averaged features)", acc1))
        summary.append((cname, "B2. proteins -> smORFs (per occurrence)", acc2))
        if cname == "all categories":
            imp = np.argsort(-m.feature_importances_)[:10]
            names = [f"{'up' if q < N_POS else 'down'}{q % N_POS + 1}_{c}" for q in range(2 * N_POS) for c in cats]
            names += [f"{'up' if q < N_POS else 'down'}{q % N_POS + 1}_{f}" for f in ("same_strand", "log_distance",
                      "overlaps") for q in range(2 * N_POS)]
            say("  most important features (protein model): "
                + ", ".join(f"{names[i]} ({m.feature_importances_[i]:.3f})" for i in imp))

        say("\nC. smORFs -> smORFs (TRAIN -> TEST, the reference)")
        t = time.time()
        m = train(xs, ks, s_train)
        acc = report(f"trained on the smORF TRAIN half ({time.time() - t:.0f}s), TEST half",
                     probabilities(m, xs, len(cats)), ks, s_test)
        summary.append((cname, "C. smORFs -> smORFs", acc))

        say("\nD. proteins + smORFs -> smORFs")
        t = time.time()
        xb = np.vstack([xp, xs])
        kb = np.vstack([kp, ks])
        m = train(xb, kb, np.r_[rows, len(xp) + s_train])
        prob = probabilities(m, xs, len(cats))
        acc = report(f"trained on all proteins and the smORF TRAIN half ({time.time() - t:.0f}s), TEST half",
                     prob, ks, s_test)
        summary.append((cname, "D. proteins + smORFs -> smORFs", acc))

    say("\nSUMMARY: accuracy % with every prediction kept (confidence >= 0)")
    for cname, what, acc in summary:
        say(f"  {cname:18s}  {what:45s}  {acc}")

    if "--protein-rules" in opts and "--smorf-rules" in opts:
        rules_side_by_side(opts["--protein-rules"], opts["--smorf-rules"], folder)
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min: {REPORT}")


if __name__ == "__main__":
    main()
