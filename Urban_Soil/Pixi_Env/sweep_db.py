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
import random
import sqlite3
import sys
import time
from collections import Counter
from dataclasses import dataclass

import numpy as np
import polars as pl

SR = ("S", "R")
NO_CALL = ["NO_EVIDENCE", "UNRESOLVED", "LOW_CONFIDENCE", "TOO_FEW_OCCURRENCES"]
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

# a neighbour is stored as [rank, [categories]]; rewrite it as an object so
# that polars can parse the JSON with a fixed type
NEIGHBOUR = r"\[(\d+), \[([^\]]*)\]\]"
NEIGHBOUR_AS_OBJECT = r'{"rank":${1},"cats":[${2}]}'
OCCS_TYPE = pl.List(pl.Struct({
    "neighbours": pl.List(pl.Struct({"rank": pl.Int64, "cats": pl.List(pl.String)})),
    "target_categories": pl.List(pl.String),
}))


@dataclass
class Data:
    smorfs: pl.DataFrame  # smorf (row number), bin, smorf_id, known_category, annotated
    targets: pl.DataFrame  # smorf, cat: known categories, without duplicates
    labels: list          # categories, then NO_CALL; methods return indices into it
    n_cats: int
    occ_smorf: np.ndarray  # smORF of each occurrence
    nb: dict              # numpy arrays, one entry per neighbour category, in file
                          # order: smorf, occ, ent (neighbour gene), rank, code (of
                          # the category), ncat (number of categories of the
                          # gene), ord (position of the gene in its occurrence
                          # when sorted by rank, counting genes without category)
    nb_by_rank: dict      # nb, with the genes of each occurrence sorted by rank

    def code(self, label):
        return self.labels.index(label)

    def excluded(self, codes, exclude):
        return np.isin(codes, [i for i, c in enumerate(self.labels[:self.n_cats]) if c in exclude])


def load_db(db_path):
    conn = sqlite3.connect(db_path)
    raw = pl.DataFrame(conn.execute(
        "SELECT bin_name, smorf_id, occurrences_json FROM smorfs ORDER BY bin_name, smorf_id"
    ).fetchall(), schema=["bin", "smorf_id", "json"], orient="row")
    conn.close()
    raw = raw.with_row_index("smorf")

    step = -(-raw.height // 64)  # parse the JSON in chunks, in parallel
    occs = (pl.concat([raw.slice(i, step).lazy().select(
        "smorf", pl.col("json").str.replace_all(NEIGHBOUR, NEIGHBOUR_AS_OBJECT)
        .str.json_decode(OCCS_TYPE).alias("o")) for i in range(0, raw.height, step)], parallel=True)
            .collect()
            .explode("o", empty_as_null=True).drop_nulls("o").with_row_index("occ").unnest("o"))

    nb = (occs.select("smorf", "occ", "neighbours")
          .explode("neighbours", empty_as_null=True).drop_nulls("neighbours")
          .with_row_index("ent").unnest("neighbours")
          .with_columns(ncat=pl.col("cats").list.len(),
                        ord=pl.col("rank").rank("ordinal").over("occ") - 1)  # ties: file order
          .explode("cats", empty_as_null=True).rename({"cats": "cat"}))
    targets = (occs.select("smorf", cat="target_categories").explode("cat", empty_as_null=True)
               .drop_nulls("cat").unique(maintain_order=True))
    cats = sorted(set(nb["cat"].drop_nulls().unique()) | set(targets["cat"].unique()))
    nb = (nb.drop_nulls("cat")  # (ord counts the genes without category too)
          .select("smorf", "occ", "ent", "ncat", "ord", rank=pl.col("rank").cast(pl.Int32),
                  code=pl.col("cat").cast(pl.Enum(cats)).to_physical().cast(pl.Int32)))

    known = targets.group_by("smorf").agg(known_category=pl.col("cat").sort().str.join(","))
    smorfs = (raw.select("smorf", "bin", "smorf_id")
              .join(known, on="smorf", how="left", maintain_order="left")
              .with_columns(pl.col("known_category").fill_null(""))
              .with_columns(annotated=pl.col("known_category") != ""))
    as_numpy = lambda df: {c: df[c].to_numpy() for c in df.columns}
    return Data(smorfs, targets, cats + NO_CALL, len(cats), occs["smorf"].to_numpy(),
                as_numpy(nb), as_numpy(nb.sort("occ", "ord", maintain_order=True)))


# ======================================================= vote helpers ==
# The methods vote for all smORFs at once, with numpy arrays, and give the
# same results as voting one smORF at a time in Python: a category's votes
# are added in the same order (np.add.at), totals are summed like Python's
# sum() (py_sum) and percentages rounded with Python's round() (round1).

def py_sum(group, values, n_groups):
    """Python's sum() of the values of each group (groups contiguous and in
    order). Since Python 3.12, sum() uses Neumaier's compensated summation."""
    starts = run_starts(group)
    size = np.diff(np.r_[starts, len(group)])
    total = np.zeros(n_groups)
    comp = np.zeros(n_groups)
    for j in range(size.max() if len(size) else 0):
        i = starts[size > j] + j  # the j-th value of each group
        g, x = group[i], values[i]
        if sys.version_info < (3, 12):
            total[g] += x
            continue
        f = total[g]
        t = f + x
        comp[g] += np.where(np.abs(f) >= np.abs(x), (f - t) + x, (x - t) + f)
        total[g] = t
    return np.where((comp != 0) & np.isfinite(comp), total + comp, total)


def run_starts(a):
    """Where each run of equal values in a starts"""
    return np.flatnonzero(np.r_[True, a[1:] != a[:-1]]) if len(a) else np.zeros(0, dtype=int)


def round1(pct):
    """Python's round(x, 1) (numpy rounds differently on ties)"""
    u, inv = np.unique(pct, return_inverse=True)
    return np.array([round(x, 1) for x in u.tolist()])[inv]


def resolve(d, group, code, v, n_groups, with_pct=True):
    """The category with most votes in each group: (code, pct) per group,
    pct not yet rounded. The votes v (for category code) are in the order they
    are cast; the total is over categories in the order of their first vote."""
    key = group.astype(np.int64) * d.n_cats + code
    votes = np.zeros(n_groups * d.n_cats)
    np.add.at(votes, key, v)
    if with_pct:
        first = np.full(n_groups * d.n_cats, len(key))
        np.minimum.at(first, key, np.arange(len(key)))
        keys = key[first[key] == np.arange(len(key))]  # by group, then first vote
    else:
        present = np.zeros(n_groups * d.n_cats, dtype=bool)
        present[key] = True
        keys = np.flatnonzero(present)  # by group, then category
    g, x = keys // d.n_cats, votes[keys]
    starts = run_starts(g)
    top = np.maximum.reduceat(x, starts) if len(g) else np.zeros(0)
    win = np.abs(x - np.repeat(top, np.diff(np.r_[starts, len(g)]))) < TIE_TOLERANCE
    nwin = np.add.reduceat(win.astype(int), starts) if len(g) else np.zeros(0, int)
    winner = np.maximum.reduceat(np.where(win, keys % d.n_cats, -1), starts) if len(g) else nwin

    voted = g[starts]
    pred = np.full(n_groups, d.code("NO_EVIDENCE"))
    pred[voted] = np.where(nwin == 1, winner, d.code("UNRESOLVED"))
    if not with_pct:
        return pred, None
    pct = np.zeros(n_groups)
    pct[voted] = top / py_sum(g, x, n_groups)[voted] * 100
    return pred, pct


def select(arrays, keep, columns):
    if keep.all():
        return {c: arrays[c] for c in columns}
    return {c: arrays[c][keep] for c in columns}


def tally(d, weights, n, exclude, multi, background):
    """The neighbour votes: arrays smorf, occ, ent, code, v with one entry per
    (neighbour, category) that votes, in the order they are cast."""
    nb = d.nb
    keep = (nb["rank"] <= n) & ~d.excluded(nb["code"], exclude)
    x = select(nb, keep, ("smorf", "occ", "ent", "rank", "code"))
    k = np.bincount(x["ent"])[x["ent"]]  # categories left for this neighbour
    if multi == "ignore":
        x, k = select(x, k == 1, x.keys()), 1
    w = np.array([np.nan] + [weights[r] for r in range(1, MAX_RANK + 1)])
    v = w[x["rank"]]
    if multi != "full":
        v = v / k
    if background:
        v = v / np.array([background.get(c, np.nan) for c in d.labels[:d.n_cats]])[x["code"]]
    x["v"] = v
    return x


# ========================================================= the methods ==
# every method: (data, **settings) -> (predicted_category, top_score_pct),
# two arrays with one value per smORF; categories are indices into d.labels

def weighted_vote(d, curve="linear", n=5, exclude=(), multi="split",
                  background=None, cutoff=None, min_occ=None):
    """Neighbours vote for their category, closer neighbours count more.
    Options: drop categories, handle multi-letter categories, down-weight
    common categories, require a minimum winning share, or require a
    minimum number of occurrences that contribute votes."""
    x = tally(d, WEIGHTS[curve], n, exclude, multi, background)
    n_smorfs = d.smorfs.height
    pred, pct = resolve(d, x["smorf"], x["code"], x["v"], n_smorfs)
    pct = round1(pct)
    if cutoff:
        pred = np.where((pred < d.n_cats) & (pct < cutoff), d.code("LOW_CONFIDENCE"), pred)
    if min_occ:
        few = np.bincount(x["smorf"][run_starts(x["occ"])], minlength=n_smorfs) < min_occ
        pred = np.where(few, d.code("TOO_FEW_OCCURRENCES"), pred)
        pct = np.where(few, 0.0, pct)
    return pred, pct


def majority(d, n=5, k=4, exclude=()):
    """Take the n closest neighbours of every occurrence; predict only if
    at least k out of n of them agree on one category."""
    nb = d.nb_by_rank
    keep = (nb["ord"] < n) & ~d.excluded(nb["code"], exclude)
    smorf, ent, code = select(nb, keep, ("smorf", "ent", "code")).values()
    n_smorfs = d.smorfs.height
    total = np.bincount(smorf[run_starts(ent)], minlength=n_smorfs)
    key = smorf.astype(np.int64) * d.n_cats + code
    counts = np.bincount(key, minlength=n_smorfs * d.n_cats).reshape(n_smorfs, d.n_cats)
    first = np.full(n_smorfs * d.n_cats, len(key))
    np.minimum.at(first, key, np.arange(len(key)))
    top = counts.max(axis=1)
    # the most common category; on ties the one counted first
    top_cat = np.where(counts == top[:, None], first.reshape(n_smorfs, d.n_cats), len(key)).argmin(axis=1)
    share = np.divide(top, total, out=np.zeros(n_smorfs), where=total > 0)
    agree = (total > 0) & (share >= k / n)
    return (np.where(agree, top_cat, d.code("NO_EVIDENCE")),
            round1(np.where(agree, share * 100, 0.0)))


def consensus(d, curve="linear", n=5, exclude=()):
    """Each occurrence (genome) makes its own prediction; then the
    occurrences vote, one vote each."""
    x = tally(d, WEIGHTS[curve], n, exclude, "split", None)
    occ_pred, _ = resolve(d, x["occ"], x["code"], x["v"], len(d.occ_smorf), with_pct=False)
    occs = np.flatnonzero(occ_pred < d.n_cats)
    pred, pct = resolve(d, d.occ_smorf[occs], occ_pred[occs], np.ones(len(occs)), d.smorfs.height)
    return pred, round1(pct)


def control_always(d, category, exclude=()):
    """CONTROL: ignore the neighbours, always predict one category."""
    n = d.smorfs.height
    return np.full(n, d.code(category)), np.zeros(n)


def control_random(d, categories, freqs, rng, exclude=()):
    """CONTROL: ignore the neighbours, random category weighted by how
    common each category is among annotated smORFs."""
    n = d.smorfs.height
    return np.array([d.code(c) for c in rng.choices(categories, freqs, k=n)]), np.zeros(n)


def control_shuffled(d, rng, exclude=()):
    """CONTROL: run the original method on the neighbours of a DIFFERENT,
    randomly chosen smORF."""
    n = d.smorfs.height
    other = np.array([rng.randrange(n) for _ in range(n)])
    pred, pct = weighted_vote(d)
    return pred[other], pct[other]


# ================================================= the 100 combinations ==

def build_combinations(d):
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
    near = d.nb["rank"] <= 5
    code = d.nb["code"][near]
    background = np.zeros(d.n_cats)
    np.add.at(background, code, 1 / d.nb["ncat"][near])
    first = np.full(d.n_cats, len(code))
    np.minimum.at(first, code, np.arange(len(code)))
    background = {d.labels[c]: float(background[c]) for c in sorted(np.flatnonzero(first < len(code)),
                                                             key=lambda c: first[c])}
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
    freq = Counter(d.targets["cat"].to_list())
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
        control_shuffled, rng=random.Random(1))
    return combos


# ===================================================== PART 1: PREDICT ==

def predict(d, combos):
    """Every combination predicts for every smORF, annotated or not."""
    preds, pcts = {}, {}
    for c in combos:
        t = time.time()
        p, s = c["fn"](d, **c["settings"])
        preds[c["id"]], pcts[c["id"]] = pl.Series(d.labels).gather(p), pl.Series(s)
        made = (p < d.n_cats).sum()
        print(f"  {c['id']}  predicted {made / d.smorfs.height * 100:5.1f}% of smORFs  "
              f"({time.time() - t:.1f}s)  {c['rule']}", flush=True)

    (d.smorfs.select("bin", "smorf_id",
                     annotated=pl.when("annotated").then(pl.lit("yes")).otherwise(pl.lit("no")),
                     known_category="known_category")
     .with_columns(preds[c["id"]].alias(c["id"]) for c in combos)
     .write_csv("predictions_all_smorfs.tsv", separator="\t", line_terminator="\r\n",
                quote_style="never"))
    return preds, pcts


# ==================================================== PART 2: EVALUATE ==

def pct(a, b):
    return f"{a / b * 100:.1f}" if b else ""


def evaluate(d, combos, preds, pcts):
    """Score each combination ONLY on smORFs with a known COG category."""
    bins = sorted(set(d.smorfs["bin"].to_list()))
    rows, by_bin = [], []
    known = {}  # exclude -> (smorf, known category) pairs

    tier = pl.lit(None, dtype=pl.String)
    for t, low in reversed(TIERS):
        tier = pl.when(pl.col("pct") >= low).then(pl.lit(t)).otherwise(tier)
    made = ~pl.col("pred").is_in(NO_CALL)
    ev = made & pl.col("has_known")
    ok = ev & pl.col("ok")
    counts = [pl.len().alias("n"), made.sum().alias("pred"), pl.col("annotated").sum().alias("annot"),
              ev.sum().alias("eval"), ok.sum().alias("correct")]
    for t, _ in TIERS:
        counts += [(ev & (tier == t)).sum().alias(f"{t} n"), (ok & (tier == t)).sum().alias(f"{t} ok")]

    for c in combos:
        ex = tuple(c["exclude"])
        if ex not in known:
            known[ex] = d.targets.filter(~pl.col("cat").is_in(list(ex)))
        k = known[ex]
        df = (d.smorfs.select("smorf", "bin", "annotated").with_columns(pred=preds[c["id"]], pct=pcts[c["id"]])
              .join(k.select("smorf", pred="cat", ok=pl.lit(True)), on=["smorf", "pred"], how="left",
                    maintain_order="left")
              .join(k.select("smorf").unique().with_columns(has_known=pl.lit(True)), on="smorf", how="left",
                    maintain_order="left")
              .with_columns(pl.col("ok", "has_known").fill_null(False)))
        st = {r["bin"]: r for r in df.group_by("bin").agg(counts).iter_rows(named=True)}
        a = df.select(counts).row(0, named=True)
        row = {"ID": c["id"], "Group": c["group"], "Rule": c["rule"],
               "smORFs predicted (% of all)": pct(a["pred"], a["n"]),
               "Coverage (% of annotated smORFs evaluated)": pct(a["eval"], a["annot"]),
               "Evaluated (n)": a["eval"], "Correct (n)": a["correct"],
               "Accuracy (%)": pct(a["correct"], a["eval"])}
        for t, _ in TIERS:
            row[f"{t} n"] = a[f"{t} n"]
            row[f"{t} accuracy (%)"] = pct(a[f"{t} ok"], a[f"{t} n"])
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

    print(f"\nAnnotated smORFs used for evaluation: {d.smorfs['annotated'].sum():,} of {d.smorfs.height:,}")
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
    d = load_db(db)
    print(f"Loaded {d.smorfs.height:,} smORFs in {time.time() - t0:.1f}s", flush=True)
    combos = build_combinations(d)

    print(f"\nPART 1 - PREDICT: {len(combos)} combinations on all smORFs", flush=True)
    preds, pcts = predict(d, combos)

    print("\nPART 2 - EVALUATE: annotated smORFs only", flush=True)
    evaluate(d, combos, preds, pcts)

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Saved:")
    print("  combinations_results.csv   the sheet: one row per rule, accuracy in %")
    print("  results_by_bin.csv         the same, per bin")
    print("  predictions_all_smorfs.tsv every smORF x every rule (C001...)")


if __name__ == "__main__":
    main()
