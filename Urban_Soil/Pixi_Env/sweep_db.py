#!/usr/bin/env python3
"""
COG function prediction for smORFs: 100 rule combinations, one file.

PART 1 - PREDICT   every combination predicts a COG category for ALL smORFs
                   (annotated or not)            -> predictions_all_smorfs.tsv
PART 2 - EVALUATE  each combination is scored ONLY on smORFs that already
                   have a known COG category     -> combinations_results.csv
                                                    results_by_bin.csv

    python3 sweep_db.py [db]        (default db: smorfs_all_bins.sqlite)

The input is the database written by build_all_bins_db.py: one row per
(bin, smORF), with a JSON list of the smORF's occurrences (one per genome
it was found in). Each occurrence lists its neighbour genes as
[rank, [COG letters]], where rank is the distance in genes from the smORF
(1 = adjacent; up to two genes per rank, one on each side), and the
smORF's own COG letters (target_categories), which are what is predicted.
A multi-letter category such as 'EG' is stored as its letters ['E', 'G'].

All results are percentages. combinations_results.csv is the sheet: one
row per rule (C001...C100) with what it does and how it scored.
C027 (all categories) and C028 (excluding S and R) are the original method.

To add a rule: add it at the end of build_combinations(). Nothing else
changes.
"""

import csv
import random
import sqlite3
import sys
import time
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
import polars as pl

SR = ("S", "R")
NO_CALL = ["NO_EVIDENCE", "UNRESOLVED", "LOW_CONFIDENCE", "TOO_FEW_OCCURRENCES"]
TIERS = (("VERY_HIGH", 50), ("HIGH", 40), ("MEDIUM", 30), ("LOW", 0))
TIE_TOLERANCE = 1e-9
PCT_TOLERANCE = 1e-9  # scores are not rounded: exactly 50% may come out as 49.99999999999999

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
    """Everything loaded from the database (see load_db).

    Categories are handled as integer codes, indices into labels. The
    neighbours are stored as flat numpy arrays (nb), one entry per
    (neighbour gene, category), so a gene with categories 'EG' has two
    entries. The methods (left, voters, counts) select and group these
    entries, and cache the results as many rules share them."""
    smorfs: pl.DataFrame  # smorf (row number), bin, smorf_id, known_category, annotated
    targets: pl.DataFrame  # smorf, cat: known categories, without duplicates
    labels: list          # categories, then NO_CALL; methods return indices into it
    n_cats: int
    occ_smorf: np.ndarray  # smORF of each occurrence
    nb: dict              # numpy arrays, one entry per neighbour category, in file
                          # order: smorf, occ, ent (neighbour gene), rank, code (of
                          # the category), key (smorf * n_cats + code), ncat (number
                          # of categories of the gene), ord (position of the gene in
                          # its occurrence when sorted by rank, counting genes
                          # without category)
    nb_by_rank: dict      # nb, with the genes of each occurrence sorted by rank
    _left: dict = field(default_factory=dict)    # cache for left()
    _voters: dict = field(default_factory=dict)  # cache for voters()
    _counts: dict = field(default_factory=dict)  # cache for counts()

    def code(self, label):
        """The code of a category or NO_CALL label"""
        return self.labels.index(label)

    def excluded(self, codes, exclude):
        """For each category code in codes: is its category in exclude?"""
        return np.array([c in exclude for c in self.labels[:self.n_cats]])[codes]

    def left(self, exclude):
        """For each entry of nb, when the categories in exclude are dropped:

        keep  is the entry kept (its category is not excluded)?
        k     how many categories its neighbour gene has left
        vote  where its vote is in the table of votes by (rank, k, code)
              (see vote_table)"""
        if exclude not in self._left:
            nb = self.nb
            keep = ~self.excluded(nb["code"], exclude)
            starts = run_starts(nb["ent"])  # the categories of a gene are next to each other
            k = np.repeat(np.add.reduceat(keep.astype(np.int32), starts), np.diff(np.r_[starts, len(keep)]))
            vote = (nb["rank"].astype(np.int32) * (nb["ncat"].max() + 1) + k) * self.n_cats + nb["code"]
            self._left[exclude] = keep, k, vote
        return self._left[exclude]

    def voters(self, n, exclude, multi, columns):
        """The entries of nb that vote: those at rank <= n whose category is
        not excluded and, if multi is "ignore", whose gene has a single
        category left. Returns the given columns of nb for these entries,
        plus vote (see left), as a new dict that the caller may modify."""
        sel = (n, exclude, multi == "ignore", columns)
        if sel not in self._voters:
            keep, k, vote = self.left(exclude)
            keep = keep & (self.nb["rank"] <= n)
            if multi == "ignore":
                keep &= k == 1
            if len(self._voters) == 2:  # keep the last two: consecutive rules often share them
                del self._voters[next(iter(self._voters))]
            self._voters[sel] = select({**self.nb, "vote": vote}, keep, columns + ("vote",))
        return dict(self._voters[sel])

    def counts(self, exclude):
        """The entries of nb kept by left(exclude), counted by (key, rank, k):
        arrays key (sorted), rank, k, vote (as in left) and n (the count).

        All the entries in a group cast the same vote, so a rule can add up
        its votes from these counts, without going through every entry."""
        if exclude not in self._counts:
            keep, _, vote = self.left(exclude)
            ks = self.nb["ncat"].max() + 1
            rk = vote[keep] // self.n_cats  # rank * ks + k
            both, n = np.unique(self.nb["key"][keep] * ((MAX_RANK + 1) * ks) + rk, return_counts=True)
            key, rk = np.divmod(both, (MAX_RANK + 1) * ks)
            rank, k = np.divmod(rk, ks)
            self._counts[exclude] = {"key": key, "rank": rank, "k": k,
                                     "vote": rk * self.n_cats + key % self.n_cats, "n": n}
        return self._counts[exclude]


def load_db(db_path):
    """Read every smORF from the database at db_path into a Data.

    smORFs are numbered (smorf) in (bin, smorf_id) order, occurrences (occ)
    and neighbour genes (ent) in file order. Neighbours without a category
    are dropped from nb, but still count for ord. A smORF is annotated if
    any of its occurrences has target categories (known_category lists
    them all, sorted and comma-separated)."""
    conn = sqlite3.connect(db_path)
    raw = pl.DataFrame(conn.execute(
        "SELECT bin_name, smorf_id, occurrences_json FROM smorfs ORDER BY bin_name, smorf_id"
    ).fetchall(), schema=["bin", "smorf_id", "json"], orient="row")
    conn.close()
    raw = raw.with_row_index("smorf")

    step = -(-raw.height // 64)  # parse the JSON in chunks, in parallel
    occs = (pl.concat([raw.slice(i, step).lazy().select(
        "smorf", pl.col("json").str.replace_all(NEIGHBOUR, NEIGHBOUR_AS_OBJECT)
        .str.json_decode(OCCS_TYPE).alias("o")).explode("o", empty_as_null=True).drop_nulls("o")
        for i in range(0, raw.height, step)], parallel=True)
            .collect()
            .with_row_index("occ").unnest("o"))

    nb = (occs.lazy().select("smorf", "occ", "neighbours")
          .explode("neighbours", empty_as_null=True).drop_nulls("neighbours")
          .with_row_index("ent").unnest("neighbours")
          .collect(engine="streaming")  # (explodes in parallel)
          .with_columns(ncat=pl.col("cats").list.len(),
                        ord=pl.col("rank").rank("ordinal").over("occ") - 1)  # ties: file order
          .explode("cats", empty_as_null=True).rename({"cats": "cat"}))
    targets = (occs.select("smorf", cat="target_categories").explode("cat", empty_as_null=True)
               .drop_nulls("cat").unique(maintain_order=True))
    cats = sorted(set(nb["cat"].drop_nulls().unique()) | set(targets["cat"].unique()))
    nb = (nb.drop_nulls("cat")  # (ord counts the genes without category too)
          .select("smorf", "occ", "ent", "ncat", "ord", rank=pl.col("rank").cast(pl.UInt8),
                  code=pl.col("cat").cast(pl.Enum(cats)).to_physical().cast(pl.UInt8))
          .with_columns(key=pl.col("smorf").cast(pl.Int64) * len(cats) + pl.col("code")))

    known = targets.group_by("smorf").agg(known_category=pl.col("cat").sort().str.join(","))
    smorfs = (raw.select("smorf", "bin", "smorf_id")
              .join(known, on="smorf", how="left", maintain_order="left")
              .with_columns(pl.col("known_category").fill_null(""))
              .with_columns(annotated=pl.col("known_category") != ""))
    as_numpy = lambda df: {c: df[c].to_numpy() for c in df.columns}
    return Data(smorfs, targets, cats + NO_CALL, len(cats), occs["smorf"].to_numpy(),
                as_numpy(nb), as_numpy(nb.sort("occ", "ord", maintain_order=True)))


# ======================================================= vote helpers ==

def run_starts(a):
    """Where each run of equal values in a starts"""
    return np.flatnonzero(np.r_[True, a[1:] != a[:-1]]) if len(a) else np.zeros(0, dtype=int)


def resolve(d, key, v, n_groups):
    """The category with most votes in each group, for groups 0...n_groups-1.

    key is group * n_cats + code for every vote, and v its weight. Returns
    two arrays with one value per group: the winning code (UNRESOLVED on a
    tie, NO_EVIDENCE if the group got no votes) and the winner's share of
    the group's votes, in % (0 if no votes)."""
    votes = np.bincount(key, weights=v, minlength=n_groups * d.n_cats)
    present = np.zeros(n_groups * d.n_cats, dtype=bool)
    present[key] = True
    keys = np.flatnonzero(present)  # by group, then category
    return pick(d, keys, votes[keys], n_groups)


def pick(d, keys, x, n_groups):
    """As resolve, but from votes already added up: x is the total vote of
    each key in keys (sorted, without duplicates)"""
    g = keys // d.n_cats
    starts = run_starts(g)
    if not len(g):
        return np.full(n_groups, d.code("NO_EVIDENCE")), np.zeros(n_groups)
    top = np.maximum.reduceat(x, starts)
    win = np.abs(x - np.repeat(top, np.diff(np.r_[starts, len(g)]))) < TIE_TOLERANCE
    nwin = np.add.reduceat(win.astype(int), starts)
    winner = np.maximum.reduceat(np.where(win, keys % d.n_cats, -1), starts)

    voted = g[starts]
    pred = np.full(n_groups, d.code("NO_EVIDENCE"))
    pred[voted] = np.where(nwin == 1, winner, d.code("UNRESOLVED"))
    pct = np.zeros(n_groups)
    pct[voted] = top / np.add.reduceat(x, starts) * 100
    return pred, pct


def select(arrays, keep, columns):
    """The given columns of the dict of arrays, only where the mask keep is True"""
    if keep.all():
        return {c: arrays[c] for c in columns}
    keep = np.flatnonzero(keep)  # (faster than indexing every column with a mask)
    return {c: arrays[c].take(keep) for c in columns}


def vote_table(d, weights, multi, background):
    """The vote of a neighbour category, by (rank, k, code), flattened so
    that it can be indexed by the vote arrays of Data.left; k is the number
    of categories left for the neighbour gene.

    The vote is the weight of the rank, divided by k if multi is "split",
    and divided by background[code] if a background is given."""
    v = np.array([np.nan] + [weights[r] for r in range(1, MAX_RANK + 1)])[:, None, None]
    ks = d.nb["ncat"].max() + 1
    if multi == "split":  # (with "ignore", k is 1)
        v = v / np.maximum(np.arange(ks), 1)[None, :, None]  # (k is never 0)
    if background is not None:
        v = v / background[None, None, :]
    return np.broadcast_to(v, (MAX_RANK + 1, ks, d.n_cats)).ravel()


def tally(d, weights, n, exclude, multi, background, columns=("key",)):
    """The neighbour votes: the given columns of nb and v (the vote, see
    vote_table), with one entry per (neighbour, category) that votes (see
    Data.voters), in file order."""
    x = d.voters(n, tuple(exclude), multi, columns)
    x["v"] = vote_table(d, weights, multi, background).take(x.pop("vote"))
    return x


# ========================================================= the methods ==
# every method: (data, **settings) -> (predicted_category, top_score_pct),
# two arrays with one value per smORF; categories are indices into d.labels
# the controls take exclude only so that they are evaluated with it

def weighted_vote(d, curve="linear", n=5, exclude=(), multi="split",
                  background=None, cutoff=None, min_occ=None):
    """Neighbours vote for their category, closer neighbours count more;
    the votes of all occurrences of a smORF are added up.

    curve       the weight of each rank, a key of WEIGHTS
    n           only neighbours with rank <= n vote
    exclude     categories that get no votes (the gene's other categories
                still do)
    multi       a gene with several categories: "split" its vote between
                them, give each the "full" vote, or "ignore" the gene
    background  if given, the vote for each category is divided by it
    cutoff      predict LOW_CONFIDENCE if the winner has < cutoff % of the
                vote
    min_occ     predict TOO_FEW_OCCURRENCES if fewer than min_occ
                occurrences contribute votes"""
    # the votes are added up from the counts of identical votes (see Data.counts)
    c = d.counts(tuple(exclude))
    use = c["rank"] <= n
    if multi == "ignore":
        use &= c["k"] == 1
    c = select(c, use, ("key", "vote", "n"))
    v = vote_table(d, WEIGHTS[curve], multi, background).take(c["vote"]) * c["n"]
    starts = run_starts(c["key"])
    n_smorfs = d.smorfs.height
    pred, pct = pick(d, c["key"][starts], np.add.reduceat(v, starts) if len(v) else v, n_smorfs)
    if cutoff:
        # predict only if pct >= cutoff; the tolerance makes a winner with exactly
        # the cutoff share pass however rounding errors in the sums fall
        low = pct < cutoff - PCT_TOLERANCE
        pred = np.where((pred < d.n_cats) & low, d.code("LOW_CONFIDENCE"), pred)
    if min_occ:
        x = d.voters(n, tuple(exclude), multi, ("smorf", "occ"))
        few = np.bincount(x["smorf"][run_starts(x["occ"])], minlength=n_smorfs) < min_occ
        pred = np.where(few, d.code("TOO_FEW_OCCURRENCES"), pred)
        pct = np.where(few, 0.0, pct)
    return pred, pct


def majority(d, n=5, k=4, exclude=()):
    """Take the n closest neighbours of every occurrence; predict the most
    common category only if at least a fraction k/n of these neighbours
    have it.

    The neighbours of all occurrences are pooled. Neighbours without a
    category still take one of the n places but are not counted, nor are
    those with only excluded categories. A gene with several categories
    counts for each. On ties, the category counted first wins."""
    nb = d.nb_by_rank
    keep = (nb["ord"] < n) & ~d.excluded(nb["code"], exclude)
    key, ent = select(nb, keep, ("key", "ent")).values()
    n_smorfs = d.smorfs.height
    total = np.bincount(key[run_starts(ent)] // d.n_cats, minlength=n_smorfs)
    counts = np.bincount(key, minlength=n_smorfs * d.n_cats).reshape(n_smorfs, d.n_cats)
    first = np.full(n_smorfs * d.n_cats, len(key))
    np.minimum.at(first, key, np.arange(len(key)))
    top = counts.max(axis=1)
    # the most common category; on ties the one counted first
    top_cat = np.where(counts == top[:, None], first.reshape(n_smorfs, d.n_cats), len(key)).argmin(axis=1)
    share = np.divide(top, total, out=np.zeros(n_smorfs), where=total > 0)
    agree = (total > 0) & (share >= k / n)
    return (np.where(agree, top_cat, d.code("NO_EVIDENCE")),
            np.where(agree, share * 100, 0.0))


def consensus(d, curve="linear", n=5, exclude=()):
    """Each occurrence (genome) makes its own prediction, by weighted vote
    of its neighbours; then the occurrences vote, one vote each (those
    with a tie do not vote). pct is the winner's share of these votes."""
    x = tally(d, WEIGHTS[curve], n, exclude, "split", None, ("occ", "code"))
    new = np.r_[True, x["occ"][1:] != x["occ"][:-1]][:len(x["occ"])]  # number the occurrences that vote 0, 1, ...
    occ_pred, _ = resolve(d, (np.cumsum(new) - 1) * d.n_cats + x["code"], x["v"], new.sum())
    made = occ_pred < d.n_cats
    smorf = d.occ_smorf[x["occ"][new][made]]
    return resolve(d, smorf.astype(np.int64) * d.n_cats + occ_pred[made], np.ones(made.sum()),
                   d.smorfs.height)


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
    """CONTROL: run the original method (C027) on the neighbours of a
    DIFFERENT, randomly chosen smORF (drawn with replacement, so it is
    occasionally the same one)."""
    n = d.smorfs.height
    other = np.array([rng.randrange(n) for _ in range(n)])
    pred, pct = weighted_vote(d)
    return pred[other], pct[other]


# ================================================= the 100 combinations ==

def build_combinations(d):
    """The list of rules, each a dict: id (C001...), group, rule (the text
    in the results), fn (the method), settings (its arguments) and exclude
    (the categories ignored, also when evaluating).

    IDs are given in order, so a rule added anywhere but at the end changes
    the IDs of the rules after it."""
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
    background = np.bincount(d.nb["code"][near], weights=1 / d.nb["ncat"][near], minlength=d.n_cats)
    seen = np.bincount(d.nb["code"][near], minlength=d.n_cats) > 0
    background = background / (background[seen].mean() if seen.any() else 1)
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
    """Every combination predicts for every smORF, annotated or not.

    Writes predictions_all_smorfs.tsv (one row per smORF, one column per
    rule) and returns the predictions and their pct, as dicts by rule ID."""
    preds, pcts = {}, {}
    for c in combos:
        t = time.time()
        preds[c["id"]], pcts[c["id"]] = p, s = c["fn"](d, **c["settings"])
        made = (p < d.n_cats).sum()
        print(f"  {c['id']}  predicted {made / d.smorfs.height * 100:5.1f}% of smORFs  "
              f"({time.time() - t:.1f}s)  {c['rule']}", flush=True)

    labels = pl.Series(d.labels)
    (d.smorfs.select("bin", "smorf_id",
                     annotated=pl.when("annotated").then(pl.lit("yes")).otherwise(pl.lit("no")),
                     known_category="known_category")
     .with_columns(labels.gather(preds[c["id"]]).alias(c["id"]) for c in combos)
     .write_csv("predictions_all_smorfs.tsv", separator="\t", line_terminator="\r\n",
                quote_style="never"))
    return preds, pcts


# ==================================================== PART 2: EVALUATE ==

def pct(a, b):
    """a / b as a percentage string, or "" if b is 0"""
    return f"{a / b * 100:.1f}" if b else ""


def evaluate(d, combos, preds, pcts):
    """Score each combination ONLY on smORFs with a known COG category.

    A smORF is evaluated if the rule predicts a category for it and it has
    a known category not excluded by the rule; the prediction is correct if
    it is any of its known categories. Coverage is the share of annotated
    smORFs that are evaluated. Accuracy is also given per confidence tier
    (TIERS, by pct). Writes combinations_results.csv and results_by_bin.csv."""
    bins = sorted(set(d.smorfs["bin"].to_list()))
    bin_of = d.smorfs["bin"].cast(pl.Enum(bins)).to_physical().to_numpy()
    annotated = d.smorfs["annotated"].to_numpy()
    annot_by_bin = np.bincount(bin_of[annotated], minlength=len(bins))
    smorf = d.targets["smorf"].to_numpy()
    cat = d.targets["cat"].cast(pl.Enum(d.labels[:d.n_cats])).to_physical().to_numpy()
    rows, by_bin = [], []
    known = {}  # exclude -> (smorf x category: is it a known category, smorf: has one)

    for c in combos:
        ex = tuple(c["exclude"])
        if ex not in known:
            k = np.zeros((d.smorfs.height, d.n_cats), dtype=bool)
            keep = ~d.excluded(cat, ex)
            k[smorf[keep], cat[keep]] = True
            known[ex] = k, k.any(axis=1)
        k, has_known = known[ex]
        p, s = preds[c["id"]], pcts[c["id"]]
        made = p < d.n_cats
        ev = made & has_known
        ok = ev & k[np.arange(len(p)), np.minimum(p, d.n_cats - 1)]
        tier = np.full(len(p), -1)
        for i, (t, low) in reversed(list(enumerate(TIERS))):
            tier = np.where(~(s < low - PCT_TOLERANCE), i, tier)  # (NaN would count as >=, as in polars)

        a = {"n": len(p), "pred": made.sum(), "annot": annotated.sum(), "eval": ev.sum(), "correct": ok.sum()}
        for i, (t, _) in enumerate(TIERS):
            a[f"{t} n"], a[f"{t} ok"] = (ev & (tier == i)).sum(), (ok & (tier == i)).sum()
        a = {x: int(y) for x, y in a.items()}
        row = {"ID": c["id"], "Group": c["group"], "Rule": c["rule"],
               "smORFs predicted (% of all)": pct(a["pred"], a["n"]),
               "Coverage (% of annotated smORFs evaluated)": pct(a["eval"], a["annot"]),
               "Evaluated (n)": a["eval"], "Correct (n)": a["correct"],
               "Accuracy (%)": pct(a["correct"], a["eval"])}
        for t, _ in TIERS:
            row[f"{t} n"] = a[f"{t} n"]
            row[f"{t} accuracy (%)"] = pct(a[f"{t} ok"], a[f"{t} n"])
        rows.append(row)
        ev_by_bin = np.bincount(bin_of[ev], minlength=len(bins))
        ok_by_bin = np.bincount(bin_of[ok], minlength=len(bins))
        for i, b in enumerate(bins):
            n_ev, n_ok, n_annot = int(ev_by_bin[i]), int(ok_by_bin[i]), int(annot_by_bin[i])
            by_bin.append({"ID": c["id"], "Rule": c["rule"], "Bin": b,
                           "Coverage (%)": pct(n_ev, n_annot), "Evaluated (n)": n_ev,
                           "Correct (n)": n_ok, "Accuracy (%)": pct(n_ok, n_ev)})

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
