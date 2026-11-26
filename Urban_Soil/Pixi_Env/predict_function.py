#!/usr/bin/env python3
"""
Predict the COG category of every smORF from its neighbour genes with 490
rules, and score each rule on the smORFs whose category is known.

    python3 predict_function.py smorf_neighbourhoods.sqlite [--save] [--report]

--save    predictions_all_smorfs.tsv (smORF x rule), combinations_results.csv
          (accuracy and coverage per rule), results_by_bin.csv
--report  report.txt (random forest by confidence; rules compared one setting
          at a time), random_forest_predictions.tsv, one_factor_comparisons.tsv

Input (from smorf_neighbourhoods_sqlite.py): one row per smORF with a JSON
list of occurrences: neighbour genes as [rank, [COG letters]], neighbour_info
(side, same strand, distance, gap, COG ID) and target_categories.
C027/C028 are the original method. Add new rules at the end of build_combinations().
"""

import contextlib
import inspect
import os
import random
import sqlite3
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
import polars as pl

SR = ("S", "R")
NO_CALL = ["NO_EVIDENCE", "UNRESOLVED", "LOW_CONFIDENCE", "TOO_FEW_OCCURRENCES"]
# groups by the winner's share of the vote: name, lowest share in %
TIERS = (("vote >=50%", 50), ("vote 40-50%", 40), ("vote 30-40%", 30), ("vote <30%", 0))
TIE_TOLERANCE = 1e-9
PCT_TOLERANCE = 1e-9  # 50% may come out as 49.99999999999999
SIM_TOLERANCE = 1e-6  # similarity 0 may come out as 3e-08 (Annoy rounds)

MAX_RANK = 10
WEIGHTS = {  # weight of the neighbour at each rank (1 = closest)
    "flat":       {r: 1.0 for r in range(1, MAX_RANK + 1)},
    "linear":     {r: max(0.0, 1.0 - 0.2 * (r - 1)) for r in range(1, MAX_RANK + 1)},  # 1, .8, .6, .4, .2
    "steeper":    {r: 0.5 ** (r - 1) for r in range(1, MAX_RANK + 1)},                # 1, .5, .25, ...
    "very_steep": {r: 0.25 ** (r - 1) for r in range(1, MAX_RANK + 1)},               # 1, .25, .06, ...
}

# JSON rewritten as objects so that polars can parse it
NEIGHBOUR = (r"\[(\d+), \[([^\]]*)\]\]", r'{"rank":${1},"cats":[${2}]}')
OCCS_TYPE = pl.List(pl.Struct({
    "neighbours": pl.List(pl.Struct({"rank": pl.Int64, "cats": pl.List(pl.String)})),
    "target_categories": pl.List(pl.String),
}))
INFO = (r'\["(up|down)", (\d+), (-?\d+), (-?\d+), ("[^"]*")\]',
        r'{"side":"${1}","same":${2},"dist":${3},"gap":${4},"cog":${5}}')
INFO_TYPE = pl.List(pl.Struct({"neighbour_info": pl.List(pl.Struct(
    {"side": pl.String, "same": pl.Int8, "dist": pl.Int64, "gap": pl.Int64, "cog": pl.String}))}))


# ---------------------------------------------------------------- data ----

@dataclass
class Data:
    """The database in memory. Categories are integer codes (indices into
    labels). nb: one entry per (neighbour gene, category) with smorf, occ, ent
    (gene), rank, code, key (smorf * n_cats + code), ncat (categories of the
    gene) and ord (position of the gene by rank)."""
    smorfs: pl.DataFrame  # smorf, bin, smorf_id, known_category, annotated
    targets: pl.DataFrame  # smorf, cat: known categories
    labels: list          # categories, then NO_CALL
    n_cats: int
    occ_smorf: np.ndarray  # smORF of each occurrence
    nb: dict
    nb_by_rank: dict      # nb sorted by rank within each occurrence
    db_path: str = ""
    info: dict = None     # per gene: occ, side, same, dist, gap, cog (None if not in the database)
    cache: dict = field(default_factory=dict)
    _left: dict = field(default_factory=dict)
    _voters: dict = field(default_factory=dict)
    _counts: dict = field(default_factory=dict)

    @property
    def n(self):
        return self.smorfs.height

    def code(self, label):
        return self.labels.index(label)

    def excluded(self, codes, exclude):
        return np.array([c in exclude for c in self.labels[:self.n_cats]])[codes]

    def left(self, exclude):
        """Per nb entry with exclude dropped: kept?, categories left on the
        gene (k), index in vote_table (vote)"""
        if exclude not in self._left:
            nb = self.nb
            keep = ~self.excluded(nb["code"], exclude)
            starts = run_starts(nb["ent"])
            k = np.repeat(np.add.reduceat(keep.astype(np.int32), starts), np.diff(np.r_[starts, len(keep)]))
            vote = (nb["rank"].astype(np.int32) * (nb["ncat"].max() + 1) + k) * self.n_cats + nb["code"]
            self._left[exclude] = keep, k, vote
        return self._left[exclude]

    def voters(self, n, exclude, multi, columns):
        """columns (plus vote) of the nb entries that vote; the last two are cached"""
        sel = (n, exclude, multi == "ignore", columns)
        if sel not in self._voters:
            keep, k, vote = self.left(exclude)
            keep = keep & (self.nb["rank"] <= n)
            if multi == "ignore":
                keep &= k == 1
            if len(self._voters) == 2:
                del self._voters[next(iter(self._voters))]
            self._voters[sel] = select({**self.nb, "vote": vote}, keep, columns + ("vote",))
        return dict(self._voters[sel])

    def counts(self, exclude):
        """Kept nb entries counted by (key, rank, k)"""
        if exclude not in self._counts:
            keep, _, vote = self.left(exclude)
            ks = self.nb["ncat"].max() + 1
            rk = vote[keep] // self.n_cats
            both, n = np.unique(self.nb["key"][keep] * ((MAX_RANK + 1) * ks) + rk, return_counts=True)
            key, rk = np.divmod(both, (MAX_RANK + 1) * ks)
            rank, k = np.divmod(rk, ks)
            self._counts[exclude] = {"key": key, "rank": rank, "k": k,
                                     "vote": rk * self.n_cats + key % self.n_cats, "n": n}
        return self._counts[exclude]


def parse(raw, pattern, dtype):
    """Parse the JSON of raw in 64 chunks, in parallel"""
    step = -(-raw.height // 64)
    return pl.concat([raw.slice(i, step).lazy().select(
        "smorf", pl.col("json").str.replace_all(*pattern).str.json_decode(dtype).alias("o"))
        .explode("o", empty_as_null=True).drop_nulls("o")
        for i in range(0, raw.height, step)], parallel=True).collect()


def load_db(db_path):
    conn = sqlite3.connect(db_path)
    raw = pl.DataFrame(conn.execute(
        "SELECT bin_name, smorf_id, occurrences_json FROM smorfs ORDER BY bin_name, smorf_id"
    ).fetchall(), schema=["bin", "smorf_id", "json"], orient="row").with_row_index("smorf")
    conn.close()

    occs = parse(raw, NEIGHBOUR, OCCS_TYPE).with_row_index("occ").unnest("o")
    nb = (occs.lazy().select("smorf", "occ", "neighbours")
          .explode("neighbours", empty_as_null=True).drop_nulls("neighbours")
          .with_row_index("ent").unnest("neighbours")
          .collect(engine="streaming")
          .with_columns(ncat=pl.col("cats").list.len(),
                        ord=pl.col("rank").rank("ordinal").over("occ") - 1)  # ties: file order
          .explode("cats", empty_as_null=True).rename({"cats": "cat"}))
    targets = (occs.select("smorf", cat="target_categories").explode("cat", empty_as_null=True)
               .drop_nulls("cat").unique(maintain_order=True))
    cats = sorted(set(nb["cat"].drop_nulls().unique()) | set(targets["cat"].unique()))
    nb = (nb.drop_nulls("cat")  # (ord still counts genes without a category)
          .select("smorf", "occ", "ent", "ncat", "ord", rank=pl.col("rank").cast(pl.UInt8),
                  code=pl.col("cat").cast(pl.Enum(cats)).to_physical().cast(pl.UInt8))
          .with_columns(key=pl.col("smorf").cast(pl.Int64) * len(cats) + pl.col("code")))

    known = targets.group_by("smorf").agg(known_category=pl.col("cat").sort().str.join(","))
    smorfs = (raw.select("smorf", "bin", "smorf_id")
              .join(known, on="smorf", how="left", maintain_order="left")
              .with_columns(pl.col("known_category").fill_null(""))
              .with_columns(annotated=pl.col("known_category") != ""))
    as_numpy = lambda df: {c: df[c].to_numpy() for c in df.columns}
    d = Data(smorfs, targets, cats + NO_CALL, len(cats), occs["smorf"].to_numpy(),
             as_numpy(nb), as_numpy(nb.sort("occ", "ord", maintain_order=True)), db_path)

    # gene strands, distances and COG IDs
    if raw.height and len(d.nb["ent"]):
        info = (parse(raw, INFO, INFO_TYPE).with_row_index("occ").unnest("o")
                .explode("neighbour_info", empty_as_null=True).drop_nulls("neighbour_info")
                .unnest("neighbour_info"))
        if info.height > d.nb["ent"].max():
            x = {c: info[c].to_numpy() for c in ("occ", "same", "dist", "gap")}
            if np.array_equal(x["occ"][d.nb["ent"]], d.nb["occ"]):  # same genes, same order
                x["side"] = (info["side"] == "down").to_numpy().astype(np.int8)
                x["cog"] = info["cog"]
                d.info = x
    if d.info is None:
        print("  (no gene distances in this database: strand/distance rules predict nothing)", flush=True)
    return d


# -------------------------------------------------------- vote helpers ----

def run_starts(a):
    """Where each run of equal values starts"""
    return np.flatnonzero(np.r_[True, a[1:] != a[:-1]]) if len(a) else np.zeros(0, dtype=int)


def per_smorf(d, key, weights=None):
    """smORF x category sums at key (smorf * n_cats + code)"""
    return np.bincount(key, weights=weights, minlength=d.n * d.n_cats).reshape(d.n, d.n_cats)


def resolve(d, key, v, n_groups):
    """Winner per group from votes v at key (group * n_cats + code): (code,
    winner's %); UNRESOLVED on a tie, NO_EVIDENCE without votes"""
    votes = np.bincount(key, weights=v, minlength=n_groups * d.n_cats)
    present = np.zeros(n_groups * d.n_cats, dtype=bool)
    present[key] = True
    keys = np.flatnonzero(present)
    return pick(d, keys, votes[keys], n_groups)


def pick(d, keys, x, n_groups):
    """As resolve, from totals x at sorted unique keys"""
    g = keys // d.n_cats
    if not len(g):
        return np.full(n_groups, d.code("NO_EVIDENCE")), np.zeros(n_groups)
    starts = run_starts(g)
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


def margin(d, keys, x, n_groups):
    """Winner's minus runner-up's % of the votes (totals x at sorted unique keys)"""
    out = np.zeros(n_groups)
    if not len(keys):
        return out
    g = keys // d.n_cats
    order = np.lexsort((-x, g))
    g, x = g[order], x[order]
    starts = run_starts(g)
    has_second = np.r_[starts[1:], len(g)] > starts + 1
    second = np.where(has_second, x[np.minimum(starts + 1, len(x) - 1)], 0.0)
    out[g[starts]] = (x[starts] - second) / np.add.reduceat(x, starts) * 100
    return out


def break_ties(d, keys, x, score, pred):
    """Where pred is UNRESOLVED: the tied category with the highest score, if only one"""
    if not len(keys):
        return pred
    g = keys // d.n_cats
    starts = run_starts(g)
    size = np.diff(np.r_[starts, len(g)])
    tied = np.abs(x - np.repeat(np.maximum.reduceat(x, starts), size)) < TIE_TOLERANCE
    s = np.where(tied, score, -np.inf)
    win = tied & (np.abs(s - np.repeat(np.maximum.reduceat(s, starts), size)) < TIE_TOLERANCE)
    nwin = np.add.reduceat(win.astype(int), starts)
    winner = np.maximum.reduceat(np.where(win, keys % d.n_cats, -1), starts)
    voted = g[starts]
    fix = (pred[voted] == d.code("UNRESOLVED")) & (nwin == 1)
    pred = pred.copy()
    pred[voted[fix]] = winner[fix]
    return pred


def select(arrays, keep, columns):
    if keep.all():
        return {c: arrays[c] for c in columns}
    keep = np.flatnonzero(keep)
    return {c: arrays[c].take(keep) for c in columns}


def vote_table(d, weights, multi, background):
    """Vote by (rank, k, code): rank weight, / k if "split", / background"""
    v = np.array([np.nan] + [weights[r] for r in range(1, MAX_RANK + 1)])[:, None, None]
    ks = d.nb["ncat"].max() + 1
    if multi == "split":
        v = v / np.maximum(np.arange(ks), 1)[None, :, None]
    if background is not None:
        v = v / background[None, None, :]
    return np.broadcast_to(v, (MAX_RANK + 1, ks, d.n_cats)).ravel()


def tally(d, weights, n, exclude, multi, background, columns=("key",)):
    x = d.voters(n, tuple(exclude), multi, columns)
    x["v"] = vote_table(d, weights, multi, background).take(x.pop("vote"))
    return x


def always(d, label="NO_EVIDENCE"):
    return np.full(d.n, d.code(label)), np.zeros(d.n)


def low_confidence(d, pred, low):
    return np.where((pred < d.n_cats) & low, d.code("LOW_CONFIDENCE"), pred)


# --------------------------------------------------------- the methods ----
# each returns (predicted code, winner's %) per smORF

def weighted_vote(d, curve="linear", n=5, exclude=(), multi="split",
                  background=None, cutoff=None, min_occ=None, margin_cutoff=None, tiebreak=None):
    """Neighbours (rank <= n) of all occurrences vote, weighted by curve.
    multi: "split", "full" or "ignore" multi-letter genes. cutoff: winner
    needs >= cutoff %. min_occ: votes from >= min_occ occurrences.
    margin_cutoff: winner needs >= this % more than the runner-up. tiebreak:
    "closest" (closest gene's category) or "common" (most common among neighbours)."""
    c = d.counts(tuple(exclude))
    use = c["rank"] <= n
    if multi == "ignore":
        use &= c["k"] == 1
    c = select(c, use, ("key", "vote", "n"))
    v = vote_table(d, WEIGHTS[curve], multi, background).take(c["vote"]) * c["n"]
    starts = run_starts(c["key"])
    keys, totals = c["key"][starts], np.add.reduceat(v, starts) if len(v) else v
    pred, pct = pick(d, keys, totals, d.n)
    if tiebreak == "closest":
        x = d.voters(n, tuple(exclude), multi, ("key", "rank"))
        closest = np.full(d.n * d.n_cats, MAX_RANK + 1)
        np.minimum.at(closest, x["key"], x["rank"])
        pred = break_ties(d, keys, totals, -closest[keys].astype(float), pred)
    elif tiebreak == "common":
        common = np.bincount(d.nb["code"][d.nb["rank"] <= n], minlength=d.n_cats).astype(float)
        pred = break_ties(d, keys, totals, common[keys % d.n_cats], pred)
    if cutoff:
        pred = low_confidence(d, pred, pct < cutoff - PCT_TOLERANCE)
    if margin_cutoff:
        pred = low_confidence(d, pred, margin(d, keys, totals, d.n) < margin_cutoff - PCT_TOLERANCE)
    if min_occ:
        x = d.voters(n, tuple(exclude), multi, ("smorf", "occ"))
        few = np.bincount(x["smorf"][run_starts(x["occ"])], minlength=d.n) < min_occ
        pred = np.where(few, d.code("TOO_FEW_OCCURRENCES"), pred)
        pct = np.where(few, 0.0, pct)
    return pred, pct


def majority(d, n=5, k=4, exclude=()):
    """Most common category among the n closest neighbours of all
    occurrences, if >= k/n of them have it (ties: first counted)"""
    nb = d.nb_by_rank
    key, ent = select(nb, (nb["ord"] < n) & ~d.excluded(nb["code"], exclude), ("key", "ent")).values()
    total = np.bincount(key[run_starts(ent)] // d.n_cats, minlength=d.n)
    counts = per_smorf(d, key)
    first = np.full(d.n * d.n_cats, len(key))
    np.minimum.at(first, key, np.arange(len(key)))
    top = counts.max(axis=1)
    top_cat = np.where(counts == top[:, None], first.reshape(d.n, d.n_cats), len(key)).argmin(axis=1)
    share = np.divide(top, total, out=np.zeros(d.n), where=total > 0)
    agree = (total > 0) & (share >= k / n)
    return np.where(agree, top_cat, d.code("NO_EVIDENCE")), np.where(agree, share * 100, 0.0)


def consensus(d, curve="linear", n=5, exclude=()):
    """Each occurrence predicts by weighted vote, then occurrences vote"""
    x = tally(d, WEIGHTS[curve], n, exclude, "split", None, ("occ", "code"))
    new = np.r_[True, x["occ"][1:] != x["occ"][:-1]][:len(x["occ"])]
    occ_pred, _ = resolve(d, (np.cumsum(new) - 1) * d.n_cats + x["code"], x["v"], new.sum())
    made = occ_pred < d.n_cats
    smorf = d.occ_smorf[x["occ"][new][made]]
    return resolve(d, smorf.astype(np.int64) * d.n_cats + occ_pred[made], np.ones(made.sum()), d.n)


def control_always(d, category, exclude=()):
    return always(d, category)


def control_random(d, categories, freqs, rng, exclude=()):
    """A random category, weighted by how common it is"""
    return np.array([d.code(c) for c in rng.choices(categories, freqs, k=d.n)]), np.zeros(d.n)


def control_shuffled(d, rng, exclude=()):
    """The original method on another random smORF's neighbours"""
    other = np.array([rng.randrange(d.n) for _ in range(d.n)])
    pred, pct = weighted_vote(d)
    return pred[other], pct[other]


def presence(d, n, exclude):
    """smORF x category: occurrences with the category within rank n; and
    occurrences per smORF (cached)"""
    key = ("presence", n, tuple(exclude))
    if key not in d.cache:
        nb = d.nb
        keep = (nb["rank"] <= n) & ~d.excluded(nb["code"], exclude)
        pair = np.unique(nb["occ"][keep].astype(np.int64) * d.n_cats + nb["code"][keep])
        sk = d.occ_smorf[pair // d.n_cats].astype(np.int64) * d.n_cats + pair % d.n_cats
        d.cache[key] = per_smorf(d, sk), np.bincount(d.occ_smorf, minlength=d.n)
    return d.cache[key]


def top_of(d, score, eligible):
    """Highest-scoring eligible category per row: (code, column)"""
    s = np.where(eligible, score, -np.inf)
    top = np.where(eligible.any(axis=1), s.max(axis=1), 0.0)
    win = eligible & (np.abs(s - top[:, None]) < TIE_TOLERANCE)
    pred = np.where(win.sum(axis=1) == 1, s.argmax(axis=1), d.code("UNRESOLVED"))
    return np.where(eligible.any(axis=1), pred, d.code("NO_EVIDENCE")), s.argmax(axis=1)


def conserved(d, n=3, exclude=(), min_frac=None):
    """Category within rank n in the most occurrences (if in >= min_frac)"""
    count, n_occ = presence(d, n, exclude)
    pred, col = top_of(d, count, count > 0)
    frac = count[np.arange(d.n), col] / np.maximum(n_occ, 1)
    pct = np.where(pred < d.n_cats, frac * 100, 0.0)
    if min_frac:
        pred = low_confidence(d, pred, frac < min_frac - PCT_TOLERANCE)
    return pred, pct


def same_cog_id_occurrences(d):
    """Per occurrence: are the two rank-1 genes (one per side) in the same COG ID?"""
    nb = d.nb
    first = run_starts(nb["ent"])  # one entry per gene
    r1 = first[nb["rank"][first] == 1]
    occ = nb["occ"][r1].astype(np.int64)
    cog = d.info["cog"].gather(nb["ent"][r1]).to_numpy()
    order = np.argsort(occ, kind="stable")
    occ, cog = occ[order], cog[order]
    starts = run_starts(occ)
    two = starts[np.diff(np.r_[starts, len(occ)]) == 2]  # a rank-1 gene on both sides
    ok = np.zeros(len(d.occ_smorf), dtype=bool)
    ok[occ[two][(cog[two] == cog[two + 1]) & (cog[two] != "")]] = True
    return ok


def sandwiched(d, exclude=(), min_genomes=1, min_frac=None, same_strand=False, max_dist=None,
               same_cog_id=False, tiebreak=None):
    """Closest gene on both sides shares a category, in >= min_genomes
    occurrences (and >= min_frac of them). Optionally both on the target's
    strand, within max_dist bp, in the same COG ID. tiebreak="closest": on a
    tie, the category whose sandwiching genes are closer on average."""
    nb = d.nb
    keep = (nb["rank"] == 1) & ~d.excluded(nb["code"], exclude)
    if same_strand or max_dist is not None or same_cog_id or tiebreak:
        if d.info is None:
            return always(d)
        ok = np.ones(len(d.info["occ"]), dtype=bool)
        if same_strand:
            ok &= d.info["same"] == 1
        if max_dist is not None:
            ok &= d.info["dist"] <= max_dist
        keep &= ok[nb["ent"]]
        if same_cog_id:
            keep &= same_cog_id_occurrences(d)[nb["occ"]]
    pair, inverse, k = np.unique(nb["occ"][keep].astype(np.int64) * d.n_cats + nb["code"][keep],
                                 return_inverse=True, return_counts=True)
    both = pair[k >= 2]  # on both rank-1 genes
    key = d.occ_smorf[both // d.n_cats].astype(np.int64) * d.n_cats + both % d.n_cats
    count = per_smorf(d, key)
    n_occ = np.bincount(d.occ_smorf, minlength=d.n)
    eligible = count >= max(min_genomes, 1)
    pred, col = top_of(d, count, eligible)
    if tiebreak == "closest":
        dist = np.bincount(inverse, weights=d.info["dist"][nb["ent"][keep]].astype(float), minlength=len(pair))
        mean = np.divide(per_smorf(d, key, dist[k >= 2]), 2 * count, out=np.full(count.shape, np.inf),
                         where=count > 0)
        tied = pred == d.code("UNRESOLVED")
        top = np.where(eligible, count, -1).max(axis=1)
        m = np.where(eligible & (count == top[:, None]) & tied[:, None], mean, np.inf)
        best = m.min(axis=1)
        close = np.isfinite(m) & (np.abs(m - np.where(np.isfinite(best), best, 0.0)[:, None]) < TIE_TOLERANCE)
        fix = tied & (close.sum(axis=1) == 1) & np.isfinite(best)
        col = np.where(fix, m.argmin(axis=1), col)
        pred = np.where(fix, col, pred)
    frac = count[np.arange(d.n), col] / np.maximum(n_occ, 1)
    if min_frac:
        pred = low_confidence(d, pred, frac < min_frac - PCT_TOLERANCE)
    return pred, np.where(pred < d.n_cats, frac * 100, 0.0)


def enriched(d, n=3, exclude=(), min_genomes=2):
    """Category most enriched within rank n compared with all occurrences"""
    count, n_occ = presence(d, n, exclude)
    frac = count / np.maximum(n_occ, 1)[:, None]
    background = count.sum(axis=0) / max(n_occ.sum(), 1)
    enrichment = np.divide(frac, background, out=np.zeros_like(frac), where=background > 0)
    pred, col = top_of(d, enrichment, count >= min_genomes)
    return pred, np.where(pred < d.n_cats, frac[np.arange(d.n), col] * 100, 0.0)


def conserved_vote(d, curve="very_steep", n=3, exclude=()):
    """Weighted vote times the share of occurrences that have the category"""
    x = tally(d, WEIGHTS[curve], n, exclude, "split", None)
    count, n_occ = presence(d, n, exclude)
    scaled = (per_smorf(d, x["key"], x["v"]) * (count / np.maximum(n_occ, 1)[:, None])).ravel()
    keys = np.flatnonzero(count.ravel() > 0)
    return pick(d, keys, scaled[keys], d.n)


def agree(d, first, second, exclude=()):
    """Predict only where two (method, settings) agree"""
    p1, s1 = first[0](d, **first[1])
    p2, _ = second[0](d, **second[1])
    ok = (p1 < d.n_cats) & (p1 == p2)
    return np.where(ok, p1, d.code("NO_EVIDENCE")), np.where(ok, s1, 0.0)


def in_operon(d, max_gap):
    """Per gene: it and every gene up to the target are on the target's
    strand with <= max_gap bp between them"""
    info = d.info
    fail = ~((info["same"] == 1) & (info["gap"] <= max_gap))
    run = np.r_[True, (info["occ"][1:] != info["occ"][:-1]) | (info["side"][1:] != info["side"][:-1])]
    failed = np.cumsum(fail)
    starts = np.flatnonzero(run)
    return failed - np.repeat(failed[starts] - fail[starts], np.diff(np.r_[starts, len(fail)])) == 0


def rank_weight(d, curve):
    return np.array([0.0] + [WEIGHTS[curve][r] for r in range(1, MAX_RANK + 1)])[d.nb["rank"]]


def vote_where(d, gene_ok, weight=None, curve="linear", n=5, exclude=(), cutoff=None):
    """Weighted vote (split) of the genes with gene_ok and rank <= n"""
    if d.info is None:
        return always(d)
    weight = rank_weight(d, curve) if weight is None else weight
    keep, k, _ = d.left(tuple(exclude))
    use = keep & (d.nb["rank"] <= n) & gene_ok(d.info)[d.nb["ent"]]
    pred, pct = resolve(d, d.nb["key"][use], weight[use] / k[use], d.n)
    if cutoff:
        pred = low_confidence(d, pred, pct < cutoff - PCT_TOLERANCE)
    return pred, pct


def operon_vote(d, max_gap=100, curve="linear", n=5, exclude=(), cutoff=None):
    """Only genes in the target gene's operon vote"""
    return vote_where(d, lambda i: in_operon(d, max_gap), curve=curve, n=n, exclude=exclude, cutoff=cutoff)


def side_vote(d, where, curve="linear", n=5, exclude=()):
    """Only upstream ("up"), downstream ("down") or same-strand genes vote"""
    ok = {"up": lambda i: i["side"] == 0, "down": lambda i: i["side"] == 1,
          "same_strand": lambda i: i["same"] == 1}[where]
    return vote_where(d, ok, curve=curve, n=n, exclude=exclude)


def close_vote(d, side="down", max_dist=100, curve="linear", n=5, exclude=()):
    """Only genes on one side, on the target's strand and within max_dist bp vote"""
    return vote_where(d, lambda i: (i["side"] == (side == "down")) & (i["same"] == 1) & (i["dist"] <= max_dist),
                      curve=curve, n=n, exclude=exclude)


def distance_vote(d, half=500, n=5, exclude=(), max_gap=None):
    """Vote halves every `half` bp from the target gene; max_gap: operon genes only"""
    if d.info is None:
        return always(d)
    weight = 0.5 ** (np.maximum(d.info["dist"][d.nb["ent"]], 0) / half)
    ok = (lambda i: np.ones(len(i["occ"]), dtype=bool)) if max_gap is None else (lambda i: in_operon(d, max_gap))
    return vote_where(d, ok, weight, n=n, exclude=exclude)


# ----------------------------------------------- similar neighbourhoods ----

def context_vectors(d, max_rank=3, dims=512):
    """smORF x COG ID: share of occurrences with the COG ID within max_rank,
    for the dims most common COG IDs"""
    key = ("ctx", max_rank, dims)
    if key not in d.cache:
        if d.info is None:
            return None
        nb = d.nb
        first = run_starts(nb["ent"])
        use = first[nb["rank"][first] <= max_rank]
        cog = d.info["cog"].gather(nb["ent"][use]).to_numpy()
        ok = cog != ""
        ids, code = np.unique(cog[ok], return_inverse=True)
        pair = np.unique(nb["occ"][use][ok].astype(np.int64) * len(ids) + code)
        sm, c = d.occ_smorf[pair // len(ids)].astype(np.int64), pair % len(ids)
        top = np.argsort(-np.bincount(np.unique(sm * len(ids) + c) % len(ids), minlength=len(ids)))[:dims]
        col = np.full(len(ids), -1)
        col[top] = np.arange(len(top))
        keep = col[c] >= 0
        print(f"  context vectors: {len(top)} of {len(ids)} COG IDs, {100 * keep.mean():.1f}% of "
              f"(occurrence, COG ID) pairs", flush=True)
        x = np.bincount(sm[keep] * len(top) + col[c[keep]],
                        minlength=d.n * len(top)).reshape(d.n, len(top)).astype(np.float32)
        x /= np.maximum(np.bincount(d.occ_smorf, minlength=d.n), 1)[:, None]
        d.cache[key] = x
    return d.cache[key]


def letter_vectors(d, curve="linear", n=5):
    """smORF x (side, COG letter): letters of the n closest genes per side,
    weighted by curve, multi-letter genes split, averaged over occurrences"""
    key = ("letters", curve, n)
    if key not in d.cache:
        nb = d.nb
        use = nb["rank"] <= n
        side = (d.info["side"][nb["ent"][use]].astype(np.int64) if d.info is not None
                else np.zeros(use.sum(), dtype=np.int64))
        width = 2 * d.n_cats
        x = np.bincount(nb["smorf"][use].astype(np.int64) * width + side * d.n_cats + nb["code"][use],
                        weights=rank_weight(d, curve)[use] / nb["ncat"][use],
                        minlength=d.n * width).reshape(d.n, width)
        x /= np.maximum(np.bincount(d.occ_smorf, minlength=d.n), 1)[:, None]
        d.cache[key] = x.astype(np.float32)
    return d.cache[key]


def exact_nearest(x, who, m, block=1024):
    """For each row of x[who]: its m most cosine-similar rows (indices into x)
    and their similarity, by brute force. Ties go to the lower index (a tiny
    offset per column; np.argpartition is very slow with many ties)."""
    xs = x[who] / np.maximum(np.linalg.norm(x[who], axis=1, keepdims=True), 1e-12)
    m = min(m, len(who))
    nn = np.empty((len(who), m), dtype=np.int64)
    sim = np.empty((len(who), m))
    offset = np.arange(len(who)) * 1e-13  # far below float32 precision

    def one(b):
        s = (xs[b:b + block] @ xs.T).astype(np.float64) - offset
        top = np.argpartition(-s, m - 1, axis=1)[:, :m]
        ts = np.take_along_axis(s, top, axis=1)
        order = np.argsort(-ts, axis=1, kind="stable")
        top = np.take_along_axis(top, order, axis=1)
        nn[b:b + block] = who[top]
        sim[b:b + block] = np.take_along_axis(ts, order, axis=1) + offset[top]

    try:
        from threadpoolctl import threadpool_limits
        one_blas_thread = threadpool_limits(1)  # the blocks run in parallel instead
    except ImportError:
        one_blas_thread = contextlib.nullcontext()
    with one_blas_thread, ThreadPoolExecutor(min(cpus(), 16)) as pool:
        list(pool.map(one, range(0, len(who), block)))
    return nn, sim


def target_genes(d):
    """(smorf, gene, n_genes) of target genes (sample|Original_ID), sorted;
    None without a neighbours table"""
    if "targets" not in d.cache:
        d.cache["targets"] = None
        conn = sqlite3.connect(d.db_path)
        try:
            rows = conn.execute("SELECT occ_id, original_id FROM neighbours WHERE is_target = 1").fetchall()
        except sqlite3.Error:
            rows = []
            print("  (no neighbours table: smORFs on the same gene may vote for each other)", flush=True)
        conn.close()
        if rows:  # occ_id is smorf_id|sample|contig|start
            t = (pl.DataFrame(rows, schema=["occ_id", "original_id"], orient="row")
                 .with_columns(parts=pl.col("occ_id").str.split("|"))
                 .select(smorf_id=pl.col("parts").list.get(0),
                         gene=pl.col("parts").list.get(1) + "|" + pl.col("original_id"))
                 .join(d.smorfs.select("smorf_id", smorf=pl.col("smorf").cast(pl.Int64)), on="smorf_id"))
            gene = t["gene"].cast(pl.Categorical).to_physical().to_numpy().astype(np.int64)
            pair = np.unique(t["smorf"].to_numpy() * (gene.max() + 1) + gene)
            d.cache["targets"] = (pair // (gene.max() + 1), pair % (gene.max() + 1), gene.max() + 1)
    return d.cache["targets"]


def share_gene(i, j, targets, block=1_000_000):
    """Do smORFs i[p] and j[p] share a target gene?"""
    smorf, gene, n_genes = targets
    key = smorf * n_genes + gene
    count = np.bincount(smorf, minlength=max(int(smorf.max()) + 1, int(i.max()) + 1, int(j.max()) + 1))
    start = np.r_[0, np.cumsum(count)][:-1]
    out = np.zeros(len(i), dtype=bool)
    for b in range(0, len(i), block):
        ii, jj = i[b:b + block], j[b:b + block]
        n = count[ii]
        pair = np.repeat(np.arange(len(ii)), n)
        g = gene[np.repeat(start[ii], n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))]
        want = jj[pair] * n_genes + g
        hit = key[np.minimum(np.searchsorted(key, want), len(key) - 1)] == want
        out[b:b + block] = np.bincount(pair[hit], minlength=len(ii)) > 0
    return out


def cpus():
    """CPUs for this job (PBS NCPUS if set)"""
    if os.environ.get("NCPUS", "").isdigit():
        return max(int(os.environ["NCPUS"]), 1)
    try:
        return max(len(os.sched_getaffinity(0)), 1)
    except AttributeError:
        return os.cpu_count() or 1


def nearest(d, engine, k, max_rank, dims, vectors="cog_ids", **opts):
    """For each smORF with a vector (who): its k + 1 most similar smORFs (nn,
    -1 if none) and cosine similarity (sim). engine: "annoy", "voyager"
    (approximate), "exact" or "faiss". vectors: "cog_ids" or "letters". Cached."""
    x = context_vectors(d, max_rank, dims) if vectors == "cog_ids" else letter_vectors(d)
    if x is None:
        return None
    key = (engine, k, max_rank, dims, vectors) + tuple(sorted(opts.items()))
    if key in d.cache:
        return d.cache[key]
    who = np.flatnonzero(x.any(axis=1))
    m = k + 1
    nn = np.full((len(who), m), -1, dtype=np.int64)
    sim = np.zeros((len(who), m))
    if engine == "annoy":
        try:
            from annoy import AnnoyIndex
        except ImportError:
            print("  (python-annoy not installed: the Annoy rules predict nothing)", flush=True)
            return None
        ikey = ("annoy index", max_rank, dims, vectors, opts["n_trees"])
        if ikey not in d.cache:
            index = AnnoyIndex(x.shape[1], "angular")
            index.set_seed(0)
            for i in who:
                index.add_item(int(i), x[i])
            index.build(opts["n_trees"], n_jobs=1)
            d.cache[ikey] = index
        index = d.cache[ikey]

        def search(rows):
            for row in rows:
                got, dj = index.get_nns_by_item(int(who[row]), m, include_distances=True)
                nn[row, :len(got)] = got
                sim[row, :len(got)] = 1.0 - np.square(dj) / 2  # angular distance -> cosine

        with ThreadPoolExecutor(cpus()) as pool:  # same results as one thread
            list(pool.map(search, np.array_split(np.arange(len(who)), 16 * cpus())))
    elif engine == "voyager":
        try:
            from voyager import Index, Space
        except ImportError:
            print("  (voyager not installed: the Voyager rules predict nothing)", flush=True)
            return None
        ikey = ("voyager index", max_rank, dims, vectors, opts["M"], opts["ef_construction"])
        if ikey not in d.cache and len(who):
            index = Index(Space.Cosine, num_dimensions=x.shape[1], M=opts["M"],
                          ef_construction=opts["ef_construction"], random_seed=1, max_elements=len(who))
            index.add_items(x[who], ids=who.tolist(), num_threads=1)  # one thread: same index every run
            d.cache[ikey] = index
        if len(who):
            m = min(m, len(who))
            got, dj = d.cache[ikey].query(x[who], k=m, query_ef=max(opts["query_ef"], m))
            nn[:, :m] = got.astype(np.int64)
            sim[:, :m] = 1.0 - dj.astype(np.float64)  # cosine distance -> cosine
    elif engine in ("exact", "faiss"):
        # searched once for the largest k used (25), then cut
        ikey = (f"{engine} search", max_rank, dims, vectors)
        if ikey not in d.cache or d.cache[ikey][0].shape[1] < m:
            if engine == "exact":
                d.cache[ikey] = exact_nearest(x, who, max(m, 26)) if len(who) else (nn, sim)
            else:
                try:
                    import faiss
                except ImportError:
                    print("  (faiss not installed: the FAISS rules predict nothing)", flush=True)
                    return None
                if len(who):  # inner product of normalised vectors = cosine
                    faiss.omp_set_num_threads(cpus())
                    xs = np.ascontiguousarray(x[who], dtype=np.float32)
                    faiss.normalize_L2(xs)
                    index = faiss.IndexFlatIP(xs.shape[1])
                    index.add(xs)
                    s, got = index.search(xs, min(max(m, 26), len(who)))
                    d.cache[ikey] = (np.where(got >= 0, who[np.maximum(got, 0)], -1), s.astype(np.float64))
                else:
                    d.cache[ikey] = (nn, sim)
        got, s = d.cache[ikey]
        w = min(m, got.shape[1])
        nn[:, :w], sim[:, :w] = got[:, :w], s[:, :w]
    else:
        raise ValueError(f"unknown search engine: {engine}")
    d.cache[key] = (who, nn, np.where(nn >= 0, sim, 0.0))
    return d.cache[key]


def known_codes(d, exclude=()):
    """Known (smorf, code) pairs, not excluded, sorted by smORF"""
    t = (d.targets.filter(~pl.col("cat").is_in(list(exclude)))
         .select(pl.col("smorf").cast(pl.Int64), code=pl.col("cat").cast(pl.Enum(d.labels[:d.n_cats]))
                 .to_physical().cast(pl.Int64))
         .sort("smorf", maintain_order=True))
    return t["smorf"].to_numpy(), t["code"].to_numpy()


def similar_vote(d, engine="annoy", k=10, max_rank=3, dims=512, exclude=(), vectors="cog_ids", **opts):
    """The k smORFs with the most similar neighbourhoods vote for their known
    categories, weighted by similarity. Not the smORF itself, not smORFs on
    the same target gene, not smORFs with nothing in common."""
    defaults = {"annoy": {"n_trees": 20}, "voyager": {"M": 12, "ef_construction": 200, "query_ef": 100}}
    found = nearest(d, engine, k, max_rank, dims, vectors, **{**defaults.get(engine, {}), **opts})
    if found is None:
        return always(d)
    who, nn, sim = found
    other = (nn >= 0) & (nn != who[:, None])
    use = other & (np.cumsum(other, axis=1) <= k)  # first k others
    i = np.repeat(who, use.sum(axis=1)).astype(np.int64)
    j, s = nn[use], sim[use]

    known_smorf, known_code = known_codes(d, exclude)
    n_known = np.bincount(known_smorf, minlength=d.n)
    first_known = np.r_[0, np.cumsum(n_known)][:-1]
    ok = (n_known[j] > 0) & (s > SIM_TOLERANCE)
    if (targets := target_genes(d)) is not None:
        ok &= ~share_gene(i, j, targets)
    i, j, s = i[ok], j[ok], s[ok]
    if not len(i):
        return always(d)
    n = n_known[j]  # one vote per category of the neighbour
    pos = np.repeat(first_known[j], n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    return resolve(d, np.repeat(i, n) * d.n_cats + known_code[pos], np.repeat(s / n, n), d.n)


# ------------------------------------------------------ machine learning ----

FOREST_POSITIONS = 5  # genes per side used as features
FOREST_SEED = 0       # also the seed of the TRAIN/TEST halves in cascade.py
FOREST_DROP = {"strand": ("_same_strand",), "distance": ("_log_distance", "_overlaps")}
GROUP_SIM = 0.9       # annotated smORFs at least this similar go in the same half
ML_MODELS = {  # name used in the rules: text in the rule
    "boosting":     "gradient boosting (HistGradientBoosting)",
    "extra_trees":  "extra trees (more randomised random forest)",
    "logistic":     "logistic regression",
    "naive_bayes":  "naive Bayes (Gaussian)",
    "knn":          "k-nearest neighbours (25 most similar smORFs, closer ones count more)",
    "mlp":          "neural network (MLP, one hidden layer of 128)",
}
ENSEMBLE = ("forest", "boosting", "extra_trees", "logistic", "mlp")


def ml_features(d):
    """smORF x features: per side and position 1-5, the gene's COG letters
    (one-hot, multi-letter genes split), same strand, log distance, overlap
    with the gene before; averaged over occurrences; plus log occurrences. Cached."""
    if "features" not in d.cache:
        nb = d.nb
        n_slot, nc = 2 * FOREST_POSITIONS, d.n_cats
        use = nb["rank"] <= FOREST_POSITIONS
        side = (d.info["side"][nb["ent"]].astype(np.int64) if d.info is not None
                else np.zeros(len(nb["ent"]), dtype=np.int64))
        key = nb["smorf"].astype(np.int64) * n_slot + side * FOREST_POSITIONS + nb["rank"].astype(np.int64) - 1
        feats = [np.bincount(key[use] * nc + nb["code"][use], weights=1 / nb["ncat"][use],
                             minlength=d.n * n_slot * nc).reshape(d.n, n_slot * nc)]
        if d.info is not None:
            first = run_starts(nb["ent"])
            first = first[use[first]]  # one entry per gene
            gene = nb["ent"][first]
            for value in (d.info["same"][gene].astype(float), np.log1p(np.maximum(d.info["dist"][gene], 0)),
                          (d.info["gap"][gene] < 0).astype(float)):
                feats.append(np.bincount(key[first], weights=value, minlength=d.n * n_slot).reshape(d.n, n_slot))
        n_occ = np.maximum(np.bincount(d.occ_smorf, minlength=d.n), 1)
        x = np.hstack(feats) / n_occ[:, None]
        d.cache["features"] = np.hstack([x, np.log(n_occ)[:, None]]).astype(np.float32)
    return d.cache["features"]


def feature_names(d):
    slots = [f"{'up' if s < FOREST_POSITIONS else 'down'}{s % FOREST_POSITIONS + 1}"
             for s in range(2 * FOREST_POSITIONS)]
    names = [f"{s}_{c}" for s in slots for c in d.labels[:d.n_cats]]
    if d.info is not None:
        names += [f"{s}_{f}" for f in ("same_strand", "log_distance", "overlaps") for s in slots]
    return names + ["log_occurrences"]


def known_matrix(d, exclude=()):
    """smORF x category: known and not excluded? Cached."""
    key = ("known", tuple(exclude))
    if key not in d.cache:
        k = np.zeros((d.n, d.n_cats), dtype=bool)
        k[known_codes(d, exclude)] = True
        d.cache[key] = k
    return d.cache[key]


def halves(d):
    """TRAIN and TEST halves of the annotated smORFs, exactly as cascade.py"""
    annotated = [int(i) for i in np.flatnonzero(d.smorfs["annotated"].to_numpy())]
    random.Random(FOREST_SEED).shuffle(annotated)
    half = len(annotated) // 2
    return np.array(sorted(annotated[:half]), dtype=np.int64), np.array(sorted(annotated[half:]), dtype=np.int64)


def related_pairs(d, min_sim=GROUP_SIM):
    """(i, j) of related annotated smORFs: on the same target gene, or among
    each other's 25 most similar neighbourhoods (rules AL) with cosine >= min_sim"""
    key = ("related pairs", min_sim)
    if key in d.cache:
        return d.cache[key]
    annotated = d.smorfs["annotated"].to_numpy()
    i, j = [np.zeros(0, dtype=np.int64)], [np.zeros(0, dtype=np.int64)]
    if (found := nearest(d, "faiss", 25, 3, 512) or nearest(d, "exact", 25, 3, 512)) is not None:
        who, nn, sim = found
        a, b, s = np.repeat(who, nn.shape[1]), nn.ravel(), sim.ravel()
        ok = (b >= 0) & (a != b) & (s >= min_sim - SIM_TOLERANCE)
        a, b = a[ok], b[ok]
        ok = annotated[a] & annotated[b]
        i.append(a[ok])
        j.append(b[ok])
    if (targets := target_genes(d)) is not None:  # link each smORF to the next on the same gene
        smorf, gene, _ = targets
        order = np.argsort(gene, kind="stable")
        smorf, gene = smorf[order], gene[order]
        same = (gene[1:] == gene[:-1]) & (smorf[1:] != smorf[:-1]) & annotated[smorf[1:]] & annotated[smorf[:-1]]
        i.append(smorf[:-1][same])
        j.append(smorf[1:][same])
    d.cache[key] = np.concatenate(i), np.concatenate(j)
    return d.cache[key]


def related_groups(d, min_sim=GROUP_SIM):
    """Group per smORF, CD-HIT style: annotated smORFs in random order; each
    one not yet grouped starts a group with its ungrouped related smORFs (no chains)"""
    key = ("related groups", min_sim)
    if key in d.cache:
        return d.cache[key]
    i, j = related_pairs(d, min_sim)
    a, b = np.r_[i, j], np.r_[j, i]
    order = np.argsort(a, kind="stable")
    a, b = a[order], b[order]
    start = np.searchsorted(a, np.arange(d.n + 1))
    annotated = [int(x) for x in np.flatnonzero(d.smorfs["annotated"].to_numpy())]
    random.Random(FOREST_SEED).shuffle(annotated)
    group = np.full(d.n, -1, dtype=np.int64)
    g = 0
    for x in annotated:
        if group[x] < 0:
            group[x] = g
            mates = b[start[x]:start[x + 1]]
            group[mates[group[mates] < 0]] = g
            g += 1
    rest = group < 0
    group[rest] = g + np.arange(rest.sum())
    d.cache[key] = group
    return group


def grouped_halves(d, min_sim=GROUP_SIM):
    """As halves(), but related smORFs always go in the same half: whole
    groups, in random order, fill the first half up to half the annotated smORFs"""
    key = ("grouped halves", min_sim)
    if key not in d.cache:
        annotated = np.flatnonzero(d.smorfs["annotated"].to_numpy())
        g = related_groups(d, min_sim)[annotated]
        ids, size = np.unique(g, return_counts=True)
        order = list(range(len(ids)))
        random.Random(FOREST_SEED).shuffle(order)
        start = np.cumsum(size[order]) - size[order]
        first = np.zeros(len(ids), dtype=bool)
        first[np.array(order, dtype=np.int64)[start < len(annotated) // 2]] = True
        in_first = first[np.searchsorted(ids, g)]
        d.cache[key] = annotated[in_first], annotated[~in_first]

        i, j = related_pairs(d, min_sim)
        def across(part):
            half = np.zeros(d.n, dtype=np.int8)
            half[part] = 1
            return int((half[i] != half[j]).sum())
        print(f"    (grouped halves: {len(annotated):,} annotated smORFs in {len(ids):,} groups, largest "
              f"{size.max() if len(size) else 0:,}, {(size > 1).sum():,} groups with more than one; halves of "
              f"{in_first.sum():,} and {(~in_first).sum():,}; related pairs {len(i):,}, split across the "
              f"halves {across(annotated[in_first]):,} ({pct(across(annotated[in_first]), len(i)) or 0}%), "
              f"with random halves {pct(across(halves(d)[0]), len(i)) or 0}%)", flush=True)
    return d.cache[key]


def make_model(name):
    """Untrained model; logistic, kNN and MLP standardise the features first"""
    from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.naive_bayes import GaussianNB
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.neural_network import MLPClassifier
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    n, seed = cpus(), FOREST_SEED
    return {
        "forest": lambda: RandomForestClassifier(n_estimators=300, min_samples_leaf=3, n_jobs=n, random_state=seed),
        "boosting": lambda: HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True,
                                                           random_state=seed),
        "extra_trees": lambda: ExtraTreesClassifier(n_estimators=300, min_samples_leaf=3, n_jobs=n,
                                                    random_state=seed),
        "logistic": lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
        "naive_bayes": GaussianNB,
        "knn": lambda: make_pipeline(StandardScaler(),
                                     KNeighborsClassifier(n_neighbors=25, weights="distance", n_jobs=n)),
        "mlp": lambda: make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(128,), early_stopping=True,
                                                                     max_iter=200, random_state=seed)),
    }[name]()


def fit_model(name, x, k, rows):
    """One training row per known category, weight 1 / number of categories
    (kNN takes no weights)"""
    rows = rows[k[rows].any(axis=1) & x[rows].any(axis=1)]
    r, c = np.nonzero(k[rows])
    w = 1.0 / k[rows].sum(axis=1)[r]
    model = make_model(name)
    if name == "knn":
        model.fit(x[rows[r]], c)
    elif hasattr(model, "steps"):  # pipeline: the weight goes to its last step
        model.fit(x[rows[r]], c, **{f"{model.steps[-1][0]}__sample_weight": w})
    else:
        model.fit(x[rows[r]], c, sample_weight=w)
    return model


def predict_proba(model, x, n_cats, rows):
    """rows x categories; 0 for rows without features"""
    p = np.zeros((len(x), n_cats))
    go = rows[x[rows].any(axis=1)]
    for b in range(0, len(go), 200_000):
        g = go[b:b + 200_000]
        p[np.ix_(g, model.classes_)] = model.predict_proba(x[g])
    return p


def ml_probabilities(d, name, exclude=(), split="random", drop=()):
    """smORF x category probabilities: each annotated smORF by a model trained
    on the other half (split "random" or "grouped"), the rest by a model
    trained on all annotated smORFs. drop: "strand" and/or "distance" features. Cached."""
    key = ("ml", name, tuple(exclude), split, tuple(drop))
    if key not in d.cache:
        t = time.time()
        x = ml_features(d)
        if drop:
            gone = sum((FOREST_DROP[f] for f in drop), ())
            x = x[:, [i for i, f in enumerate(feature_names(d)) if not f.endswith(gone)]]
        k = known_matrix(d, exclude)
        a, b = grouped_halves(d) if split == "grouped" else halves(d)
        rest = np.flatnonzero(~d.smorfs["annotated"].to_numpy())
        prob = np.zeros((d.n, d.n_cats))
        for train, test in ((a, b), (b, a), (np.r_[a, b], rest)):
            if len(train) and len(test):
                prob += predict_proba(fit_model(name, x, k, train), x, d.n_cats, test)
        d.cache[key] = prob
        print(f"    ({name} trained 3 times in {time.time() - t:.1f}s)", flush=True)
    return d.cache[key]


def from_probabilities(d, prob, cutoff=None):
    """Most probable category; winner's % = its probability"""
    conf = prob.max(axis=1)
    pred = np.where(conf > 0, prob.argmax(axis=1), d.code("NO_EVIDENCE"))
    if cutoff:
        pred = low_confidence(d, pred, conf * 100 < cutoff - PCT_TOLERANCE)
    return pred, conf * 100


def have_sklearn():
    try:
        import sklearn  # noqa: F401
        return True
    except ImportError:
        print("  (scikit-learn not installed: the ML rules predict nothing)", flush=True)
        return False


def ml_model(d, model="boosting", exclude=(), cutoff=None, split="random", drop=()):
    """One model (forest or ML_MODELS); confidence = its probability"""
    if not have_sklearn():
        return always(d)
    return from_probabilities(d, ml_probabilities(d, model, exclude, split, drop), cutoff)


def random_forest(d, exclude=(), cutoff=None, drop=(), split="random"):
    return ml_model(d, "forest", exclude, cutoff, split, drop)


def ml_ensemble(d, models=ENSEMBLE, exclude=(), cutoff=None, split="random"):
    """Average of the models' probabilities"""
    if not have_sklearn():
        return always(d)
    return from_probabilities(d, np.mean([ml_probabilities(d, m, exclude, split) for m in models], axis=0), cutoff)


# ------------------------------------------------------------ the rules ----

CATS = (("all categories", ()), ("excluding S and R", SR))
CURVE = {"flat": "flat", "linear": "linear", "steeper": "steeper", "very_steep": "very steep"}
HALVES = "(each half predicted by a model trained on the other half)"
GROUPED = ("each half predicted by a model trained on the other half (related smORFs kept in the same half: "
           "same target gene or neighbourhood similarity >= 0.9)")
ENSEMBLE_TEXT = "average of random forest, gradient boosting, extra trees, logistic regression and neural network"
CUTS = (None, 50, 60, 70, 80, 90)


def need(cut):
    return f", predict only if confidence >= {cut}%" if cut else ""


def build_combinations(d):
    """All rules: id, group, rule, fn, settings, exclude. IDs follow the order."""
    combos = []

    def add(group, rule, fn, **settings):
        combos.append(dict(id=f"C{len(combos) + 1:03d}", group=group, rule=rule, fn=fn,
                           settings=settings, exclude=settings.get("exclude", ())))

    def both(group, rule, fn, tail="", **settings):
        """The rule for all categories, then excluding S and R"""
        for cname, ex in CATS:
            add(group, f"{rule}, {cname}{tail}", fn, exclude=ex, **settings)

    # A. weighted vote (34)
    for n in (1, 2, 3, 4, 5):
        for curve in ("linear", "steeper", "flat", "very_steep"):
            if n > 1 or curve == "linear":  # 1 neighbour: same answer for every weighting
                both("A. Weighted vote", f"{CURVE[curve]} weights, {n} neighbour(s)", weighted_vote,
                     tail=" (ORIGINAL METHOD)" if curve == "linear" and n == 5 else "", curve=curve, n=n)
    # B. majority rule (10)
    for n, k in ((3, 2), (5, 3), (5, 4), (7, 5), (9, 6)):
        both("B. Majority rule", f"at least {k} of the {n} nearest agree", majority, n=n, k=k)
    # C. confidence cutoff (20)
    for curve in ("linear", "steeper"):
        for cut in (30, 40, 50, 60, 70):
            both("C. Confidence cutoff", f"{curve} weights, 5 neighbours, predict only if winner has >= {cut}% "
                 f"of the vote", weighted_vote, curve=curve, n=5, cutoff=cut)
    # D. down-weight common categories (6)
    near = d.nb["rank"] <= 5
    background = np.bincount(d.nb["code"][near], weights=1 / d.nb["ncat"][near], minlength=d.n_cats)
    seen = np.bincount(d.nb["code"][near], minlength=d.n_cats) > 0
    background = background / (background[seen].mean() if seen.any() else 1)
    for curve, n in (("linear", 5), ("steeper", 5), ("steeper", 3)):
        both("D. Down-weight common categories", f"{curve} weights, {n} neighbours, votes divided by how common "
             f"the category is overall", weighted_vote, curve=curve, n=n, background=background)
    # E. consensus across genomes (8)
    for curve in ("linear", "steeper"):
        for n in (3, 5):
            both("E. Consensus across genomes", f"each genome predicts with {curve} weights and {n} neighbours, "
                 f"then genomes vote", consensus, curve=curve, n=n)
    # F. minimum evidence (6)
    for m in (2, 3, 5):
        both("F. Minimum evidence", f"steeper weights, 5 neighbours, predict only if >= {m} genomes contribute "
             f"votes", weighted_vote, curve="steeper", n=5, min_occ=m)
    # G. more neighbours (4)
    for n in (7, 10):
        both("G. More neighbours", f"steeper weights, {n} neighbours", weighted_vote, curve="steeper", n=n)
    # H. multi-letter categories (4)
    for multi, text in (("full", "count every letter fully"), ("ignore", "ignore multi-letter genes")):
        both("H. Multi-letter categories", f"steeper weights, 5 neighbours, {text}", weighted_vote,
             curve="steeper", n=5, multi=multi)
    # I. other exclusions (4)
    for curve in ("linear", "steeper"):
        for letter in SR:
            add("I. Other exclusions", f"{curve} weights, 5 neighbours, excluding {letter} only",
                weighted_vote, curve=curve, n=5, exclude=(letter,))
    # J. controls (4)
    freq = Counter(d.targets["cat"].to_list())
    freq_no_sr = Counter({c: v for c, v in freq.items() if c not in SR})
    for f, (cname, ex) in zip((freq, freq_no_sr), CATS):
        if f:
            top = f.most_common(1)[0][0]
            add("J. Control", f"always predict '{top}' (most common category), {cname}", control_always,
                category=top, **({"exclude": ex} if ex else {}))
    if freq:
        cats, counts = zip(*freq.most_common())
        add("J. Control", "random category, weighted by how common each is, all categories",
            control_random, categories=cats, freqs=counts, rng=random.Random(0))
    add("J. Control", "original method run on the neighbours of a different random smORF",
        control_shuffled, rng=random.Random(1))
    # K. conserved neighbourhood (24)
    for n in (1, 2, 3, 5):
        for frac in (None, 0.5, 0.7):
            both("K. Conserved neighbourhood", f"category found within {n} neighbour(s) in the most genomes"
                 + (f", only if in >= {frac:.0%} of its genomes" if frac else ""), conserved, n=n, min_frac=frac)
    # L. sandwiched (6)
    for genomes, frac, text in ((1, None, "in at least 1 genome"), (2, None, "in at least 2 genomes"),
                                (1, 0.5, "in at least 50% of its genomes")):
        both("L. Sandwiched", f"closest gene on both sides has the same category, {text}", sandwiched,
             min_genomes=genomes, min_frac=frac)
    # M. two methods agree (4)
    for text, first, second in (
            ("very steep weights, 3 neighbours, and conserved neighbourhood (3 neighbours)",
             (weighted_vote, dict(curve="very_steep", n=3)), (conserved, dict(n=3))),
            ("original method and 'at least 3 of the 5 nearest agree'",
             (weighted_vote, dict(curve="linear", n=5)), (majority, dict(n=5, k=3)))):
        for cname, ex in CATS:
            add("M. Two methods agree", f"{text} predict the same category, {cname}", agree, exclude=ex,
                first=(first[0], {**first[1], "exclude": ex}), second=(second[0], {**second[1], "exclude": ex}))
    # N. enrichment (4)
    for n in (3, 5):
        both("N. Enrichment over background", f"category most enriched within {n} neighbours compared with the "
             f"whole data, found in >= 2 genomes", enriched, n=n)
    # O. conservation-weighted vote (4)
    for curve, n in (("very_steep", 3), ("linear", 5)):
        both("O. Conservation-weighted vote", f"{CURVE[curve]} weights, {n} neighbours, each category's vote "
             f"times the share of genomes that have it", conserved_vote, curve=curve, n=n)
    # P. very steep cutoffs (8)
    for cut in (50, 60, 70, 80):
        both("P. Confidence cutoff (very steep)", f"very steep weights, 3 neighbours, predict only if winner has "
             f">= {cut}% of the vote", weighted_vote, curve="very_steep", n=3, cutoff=cut)
    # Q. operon vote (12)
    for gap in (50, 100, 200):
        for curve in ("flat", "linear"):
            both("Q. Operon vote", f"{curve} weights, only neighbours in the target gene's operon (same strand, "
                 f"<= {gap} bp between genes)", operon_vote, max_gap=gap, curve=curve, n=5)
    # R. sandwiched with strand/distance (12)
    for strand, dist, text in ((True, None, "both on the target's strand"),
                               (False, 100, "both within 100 bp of the target gene"),
                               (True, 100, "both on the target's strand and within 100 bp")):
        for genomes, frac, need_g in ((2, None, "in at least 2 genomes"), (1, 0.5, "in at least 50% of its genomes")):
            both("R. Sandwiched (strand/distance)", f"closest gene on both sides has the same category, {text}, "
                 f"{need_g}", sandwiched, same_strand=strand, max_dist=dist, min_genomes=genomes, min_frac=frac)
    # S. distance-weighted vote (6)
    for half in (250, 500, 1000):
        both("S. Distance-weighted vote", f"5 neighbours, vote halves every {half} bp from the target gene",
             distance_vote, half=half, n=5)
    # T. operon vote + cutoff (4)
    for cut in (50, 70):
        both("T. Operon vote + cutoff", f"linear weights, operon neighbours (<= 100 bp), predict only if winner "
             f"has >= {cut}% of the vote", operon_vote, max_gap=100, curve="linear", n=5, cutoff=cut)
    # U, V. similar neighbourhoods: Annoy (4), Voyager (4)
    for group, engine, name in (("U. Annoy nearest neighbours", "annoy", "Annoy"),
                                ("V. Voyager nearest neighbours", "voyager", "Voyager")):
        for k in (10, 25):
            both(group, f"the {k} smORFs with the most similar neighbourhoods (COG IDs within 3 genes, {name}) "
                 f"vote", similar_vote, engine=engine, k=k)
    # W. overlapping / tightly coupled genes (6)
    for gap, text in ((-1, "overlapping"), (0, "overlapping or touching"), (20, "<= 20 bp apart")):
        both("W. Coupled genes", f"linear weights, only neighbours in a run of same-strand genes {text} up to "
             f"the target gene", operon_vote, max_gap=gap, curve="linear", n=5)
    # X. one side only (4)
    for side in ("up", "down"):
        both("X. One side only", f"linear weights, 5 neighbours, only {side}stream genes vote", side_vote,
             where=side, curve="linear", n=5)
    # Y. same strand, any distance (2)
    both("Y. Same strand", "linear weights, 5 neighbours, only genes on the target's strand vote (any distance)",
         side_vote, where="same_strand", curve="linear", n=5)
    # Z. operon + distance weighting (4)
    for half in (250, 500):
        both("Z. Operon + distance-weighted vote", f"5 neighbours, operon neighbours (<= 100 bp), vote halves "
             f"every {half} bp from the target gene", distance_vote, half=half, n=5, max_gap=100)
    # AA. margin cutoff (12)
    for curve in ("linear", "steeper"):
        for m in (10, 20, 30):
            both("AA. Margin cutoff", f"{curve} weights, 5 neighbours, predict only if winner has >= {m}% more "
                 f"of the vote than the runner-up", weighted_vote, curve=curve, n=5, margin_cutoff=m)
    # AB. tie-breaking (8)
    for curve in ("linear", "steeper"):
        for tb, text in (("closest", "the closest gene's category"),
                         ("common", "the category most common among neighbours")):
            both("AB. Tie-breaking", f"{curve} weights, 5 neighbours, ties go to {text}", weighted_vote,
                 curve=curve, n=5, tiebreak=tb)
    # AC. exact nearest neighbours on the Annoy/Voyager vectors (4)
    for k in (10, 25):
        both("AC. Exact nearest neighbours", f"the {k} smORFs with the most similar neighbourhoods (COG IDs "
             f"within 3 genes, exact search) vote", similar_vote, engine="exact", k=k)
    # AD. exact nearest neighbours on COG letters (4)
    for k in (10, 25):
        both("AD. Similar letter profiles", f"the {k} smORFs with the most similar COG letters (5 genes per side, "
             f"linear weights, up/downstream kept apart, exact search) vote", similar_vote, engine="exact", k=k,
             vectors="letters")
    # AE. random forest (8)
    for cut in (None, 40, 50, 60):
        both("AE. Random forest", f"random forest on COG letters, strand and distance of 5 genes per side "
             f"{HALVES}{need(cut)}", random_forest, cutoff=cut)
    # AF. random forest and the strict sandwiched rule (C171/C172) agree (2)
    for cname, ex in CATS:
        add("AF. Random forest + sandwiched agree", f"random forest and 'closest gene on both sides has the same "
            f"category, both on the target's strand and within 100 bp, in at least 2 genomes' predict the same "
            f"category, {cname}", agree, exclude=ex, first=(random_forest, dict(exclude=ex)),
            second=(sandwiched, dict(exclude=ex, min_genomes=2, same_strand=True, max_dist=100)))
    # AG. random forest, stricter cutoffs (6)
    for cut in (70, 80, 90):
        both("AG. Random forest (strict)", f"random forest on COG letters, strand and distance of 5 genes per "
             f"side {HALVES}{need(cut)}", random_forest, cutoff=cut)
    # AH. sandwiched by the same COG ID (2)
    both("AH. Sandwiched (same COG ID)", "closest gene on both sides in the same COG (same COG ID), both on the "
         "target's strand and within 100 bp, in at least 2 genomes", sandwiched, min_genomes=2, same_strand=True,
         max_dist=100, same_cog_id=True)
    # AI. downstream genes only, same strand, within 100 bp (2)
    both("AI. Downstream, same strand, close", "linear weights, 5 neighbours, only downstream genes on the "
         "target's strand within 100 bp vote", close_vote, side="down", max_dist=100, curve="linear", n=5)
    # AJ. sandwiched with a tie-break (2)
    both("AJ. Sandwiched (tie-break)", "closest gene on both sides has the same category, both on the target's "
         "strand and within 100 bp, in at least 2 genomes; ties go to the category whose genes are closer",
         sandwiched, min_genomes=2, same_strand=True, max_dist=100, tiebreak="closest")
    # AK. sandwiched, very close (2)
    both("AK. Sandwiched (very close)", "closest gene on both sides has the same category, both on the target's "
         "strand and within 50 bp, in at least 2 genomes", sandwiched, min_genomes=2, same_strand=True, max_dist=50)
    # AL. the exact searches of AC/AD with FAISS (8)
    for vectors, text in (("cog_ids", "neighbourhoods (COG IDs within 3 genes, exact search with FAISS)"),
                          ("letters", "COG letters (5 genes per side, linear weights, up/downstream kept apart, "
                                      "exact search with FAISS)")):
        for k in (10, 25):
            both("AL. FAISS nearest neighbours", f"the {k} smORFs with the most similar {text} vote", similar_vote,
                 engine="faiss", k=k, vectors=vectors)
    # AM. random forest without strand and/or distance (42)
    for drop, used in ((("strand",), "COG letters and distance (no strand)"),
                       (("distance",), "COG letters and strand (no distance)"),
                       (("strand", "distance"), "COG letters only (no strand, no distance)")):
        for cut in (None, 40, 50, 60, 70, 80, 90):
            both("AM. Random forest (fewer features)", f"random forest on {used} of 5 genes per side "
                 f"{HALVES}{need(cut)}", random_forest, cutoff=cut, drop=drop)
    # AN. other ML models, same features and halves (72)
    for model, text in ML_MODELS.items():
        for cut in CUTS:
            both("AN. Other ML models", f"{text} on COG letters, strand and distance of 5 genes per side "
                 f"{HALVES}{need(cut)}", ml_model, model=model, cutoff=cut)
    # AO. ensemble (12)
    for cut in CUTS:
        both("AO. ML ensemble", f"{ENSEMBLE_TEXT} (each half predicted by models trained on the other half)"
             f"{need(cut)}", ml_ensemble, cutoff=cut)
    # AP. the ML rules again with related smORFs in the same half; a big drop
    # from AE/AG/AN/AO means those scores were inflated (84)
    for model, text in (("forest", "random forest"), ("extra_trees", ML_MODELS["extra_trees"]),
                        ("boosting", ML_MODELS["boosting"]), ("logistic", ML_MODELS["logistic"]),
                        ("knn", ML_MODELS["knn"]), ("mlp", ML_MODELS["mlp"]), ("ensemble", ENSEMBLE_TEXT)):
        fn, s = ((ml_ensemble, {}) if model == "ensemble" else
                 (random_forest, {}) if model == "forest" else (ml_model, {"model": model}))
        for cut in CUTS:
            both("AP. Grouped halves", f"{text} on COG letters, strand and distance of 5 genes per side, "
                 f"{GROUPED}{need(cut)}", fn, **s, cutoff=cut, split="grouped")
    # AQ. random forest on COG letters only, grouped halves (12)
    for cut in CUTS:
        both("AQ. Grouped halves, COG letters only", f"random forest on COG letters only (no strand, no distance) "
             f"of 5 genes per side, {GROUPED}{need(cut)}", random_forest, cutoff=cut, drop=("strand", "distance"),
             split="grouped")
    return combos


# --------------------------------------------------- predict, evaluate ----

def pct(a, b):
    return f"{a / b * 100:.1f}" if b else ""


def predict(d, combos, save=True):
    preds, pcts = {}, {}
    for c in combos:
        t = time.time()
        preds[c["id"]], pcts[c["id"]] = p, _ = c["fn"](d, **c["settings"])
        print(f"  {c['id']}  predicted {(p < d.n_cats).mean() * 100:5.1f}% of smORFs  "
              f"({time.time() - t:.1f}s)  {c['rule']}", flush=True)
    if save:
        labels = pl.Series(d.labels)
        (d.smorfs.select("bin", "smorf_id", annotated=pl.when("annotated").then(pl.lit("yes")).otherwise(pl.lit("no")),
                         known_category="known_category")
         .with_columns(labels.gather(preds[c["id"]]).alias(c["id"]) for c in combos)
         .write_csv("predictions_all_smorfs.tsv", separator="\t", quote_style="never"))
    return preds, pcts


def evaluate(d, combos, preds, pcts, save=True):
    """Score each rule on the smORFs with a known (not excluded) category:
    correct if any known category. Coverage = % of annotated smORFs evaluated."""
    bins = sorted(set(d.smorfs["bin"].to_list()))
    bin_of = d.smorfs["bin"].cast(pl.Enum(bins)).to_physical().to_numpy()
    annotated = d.smorfs["annotated"].to_numpy()
    annot_by_bin = np.bincount(bin_of[annotated], minlength=len(bins))
    rows, by_bin = [], []
    for c in combos:
        k = known_matrix(d, c["exclude"])
        p, s = preds[c["id"]], pcts[c["id"]]
        made = p < d.n_cats
        ev = made & k.any(axis=1)
        ok = ev & k[np.arange(len(p)), np.minimum(p, d.n_cats - 1)]
        tier = np.select([~(s < low - PCT_TOLERANCE) for _, low in TIERS], range(len(TIERS)), -1)  # NaN: top tier
        row = {"ID": c["id"], "Group": c["group"], "Rule": c["rule"], "Accuracy (%)": pct(ok.sum(), ev.sum()),
               "Coverage (% of annotated smORFs evaluated)": pct(ev.sum(), annotated.sum()),
               "smORFs predicted (% of all)": pct(made.sum(), len(p)),
               "Evaluated (n)": int(ev.sum()), "Correct (n)": int(ok.sum())}
        for i, (t, _) in enumerate(TIERS):
            row[f"{t} n"] = int((ev & (tier == i)).sum())
            row[f"{t} accuracy (%)"] = pct((ok & (tier == i)).sum(), row[f"{t} n"])
        rows.append(row)
        n_ev, n_ok = np.bincount(bin_of[ev], minlength=len(bins)), np.bincount(bin_of[ok], minlength=len(bins))
        by_bin += [{"ID": c["id"], "Rule": c["rule"], "Bin": b, "Coverage (%)": pct(n_ev[i], annot_by_bin[i]),
                    "Evaluated (n)": int(n_ev[i]), "Correct (n)": int(n_ok[i]), "Accuracy (%)": pct(n_ok[i], n_ev[i])}
                   for i, b in enumerate(bins)]

    ranked = sorted(rows, key=lambda r: float(r["Accuracy (%)"] or 0), reverse=True)
    for rank, r in enumerate(ranked, 1):
        r["Accuracy rank"] = rank
    if save:
        cols = ["ID", "Accuracy rank", "Group", "Rule", "Accuracy (%)", "Coverage (% of annotated smORFs evaluated)",
                "smORFs predicted (% of all)", "Evaluated (n)", "Correct (n)"] + \
               [f"{t} {x}" for t, _ in TIERS for x in ("n", "accuracy (%)")]
        pl.DataFrame(rows).select(cols).write_csv("combinations_results.csv")
        pl.DataFrame(by_bin).write_csv("results_by_bin.csv")

    cov = "Coverage (% of annotated smORFs evaluated)"
    print(f"\nAnnotated smORFs used for evaluation: {annotated.sum():,} of {d.n:,}")
    print("\nAll rules, most accurate first:")
    print(f"  {'rank':>4}  {'ID':4}  {'accuracy':>8}  {'coverage':>8}  {'evaluated':>9}  group: rule")
    for r in ranked:
        print(f"  {r['Accuracy rank']:>4}  {r['ID']}  {r['Accuracy (%)']:>7}%  {r[cov]:>7}%  "
              f"{r['Evaluated (n)']:>9,}  {r['Group']}: {r['Rule']}")
    orig = next(r for r in rows if "ORIGINAL" in r["Rule"] and "all categories" in r["Rule"])
    print(f"\nOriginal method ({orig['ID']}): accuracy {orig['Accuracy (%)']}%  coverage {orig[cov]}%")
    return rows


# -------------------------------------------------------------- report ----

CONFIDENCE_STEPS = (0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
CALIBRATION_BINS = (0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01)


def forest_report(d, out):
    """Random forest (rules AE) on the TEST half: accuracy by confidence
    cutoff and per confidence bin, top features; random_forest_predictions.tsv"""
    train, test = halves(d)
    out.append("1. RANDOM FOREST on the TEST half (predicted by the model trained on the TRAIN half)")
    for cname, ex in CATS:
        prob, k = ml_probabilities(d, "forest", ex), known_matrix(d, ex)
        conf = prob.max(axis=1)
        pred = np.where(conf > 0, prob.argmax(axis=1), d.code("NO_EVIDENCE"))
        made = (pred[test] < d.n_cats) & k[test].any(axis=1)
        right = made & k[test, np.minimum(pred[test], d.n_cats - 1)]
        c = conf[test]
        out.append(f"\n  {cname}\n  predict only if confidence >=\tcoverage %\taccuracy %")
        for step in CONFIDENCE_STEPS:
            keep = made & (c >= step)
            out.append(f"  {step:.1f}\t{pct(keep.sum(), len(test))}\t{pct((right & keep).sum(), keep.sum())}")
        out.append("  confidence bin\tTEST smORFs\taccuracy %")
        for lo, hi in zip(CALIBRATION_BINS[:-1], CALIBRATION_BINS[1:]):
            b = made & (c >= lo) & (c < hi)
            out.append(f"  {lo:.1f}-{min(hi, 1):.1f}\t{b.sum():,}\t{pct((right & b).sum(), b.sum())}")
        if not ex:
            imp = fit_model("forest", ml_features(d), k, train).feature_importances_
            names = feature_names(d)
            out.append("  most important features (model trained on TRAIN): "
                       + ", ".join(f"{names[i]} ({imp[i]:.3f})" for i in np.argsort(-imp)[:15]))
            part = np.full(d.n, "", dtype=object)
            part[train], part[test] = "TRAIN", "TEST"
            (d.smorfs.select("bin", "smorf_id", "known_category")
             .with_columns(half=pl.Series(part.tolist()), prediction=pl.Series(d.labels).gather(pred),
                           confidence=pl.Series(np.round(conf, 3)))
             .write_csv("random_forest_predictions.tsv", separator="\t"))


def shown(v):
    if isinstance(v, np.ndarray):
        return "divided by how common"
    if isinstance(v, tuple):
        return ",".join(v) if all(isinstance(x, str) for x in v) and v else ("none" if not v else str(v))
    return str(v)


def one_factor_report(combos, rows, out):
    """Every pair of rules (same method) differing in exactly one setting, with
    the change in accuracy and coverage, summarised per setting"""
    cov = "Coverage (% of annotated smORFs evaluated)"
    res = {r["ID"]: r for r in rows if r["Accuracy (%)"]}

    def full(c):  # settings with the function's defaults filled in
        return {**{k: p.default for k, p in inspect.signature(c["fn"]).parameters.items()
                   if p.default is not inspect.Parameter.empty}, **c["settings"]}

    def same(a, b):
        if isinstance(a, (np.ndarray, random.Random)) or isinstance(b, np.ndarray):
            return a is b
        try:
            return bool(a == b)
        except Exception:
            return a is b

    by_fn = {}
    for c in combos:
        if c["id"] in res and not c["group"].startswith("J."):
            by_fn.setdefault(c["fn"].__name__, []).append((c, full(c)))
    pairs = []
    for fn, cs in by_fn.items():
        for i, (ci, fi) in enumerate(cs):
            for cj, fj in cs[i + 1:]:
                diff = [k for k in fi.keys() | fj.keys() if not same(fi.get(k), fj.get(k))]
                if len(diff) != 1:
                    continue
                f = diff[0]
                (ca, va), (cb, vb) = (ci, fi.get(f)), (cj, fj.get(f))
                if (isinstance(va, (int, float)) and isinstance(vb, (int, float)) and va > vb) or \
                        (vb is None and va is not None):
                    (ca, va), (cb, vb) = (cb, vb), (ca, va)
                ra, rb = res[ca["id"]], res[cb["id"]]
                pairs.append({"method": fn, "setting": f, "from": shown(va), "to": shown(vb),
                              "rule A": ca["id"], "rule B": cb["id"],
                              "accuracy A": float(ra["Accuracy (%)"]), "accuracy B": float(rb["Accuracy (%)"]),
                              "coverage A": float(ra[cov]), "coverage B": float(rb[cov])})
    out.append("\n2. ONE FACTOR AT A TIME: rule pairs that differ in exactly one setting (all annotated smORFs)")
    if not pairs:
        out.append("  no pairs")
        return
    df = pl.DataFrame(pairs).with_columns(
        (pl.col("accuracy B") - pl.col("accuracy A")).round(1).alias("accuracy change"),
        (pl.col("coverage B") - pl.col("coverage A")).round(1).alias("coverage change"))
    df.write_csv("one_factor_comparisons.tsv", separator="\t")
    out.append(f"  {len(pairs)} pairs (all in one_factor_comparisons.tsv); median change from A to B:")
    out.append("method\tsetting\tfrom -> to\tpairs\taccuracy change\tcoverage change")
    for r in (df.group_by("method", "setting", "from", "to")
              .agg(pl.len(), pl.col("accuracy change").median(), pl.col("coverage change").median())
              .sort("method", "setting", "from", "to").iter_rows(named=True)):
        out.append(f"{r['method']}\t{r['setting']}\t{r['from']} -> {r['to']}\t{r['len']}\t"
                   f"{r['accuracy change']:+.1f}\t{r['coverage change']:+.1f}")
    out.append("  (setting 'exclude': the two rules are scored on different smORFs, so compare with care)")


def report(d, combos, rows):
    out = []
    forest_report(d, out)
    one_factor_report(combos, rows, out)
    with open("report.txt", "w") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n" + "\n".join(out[:40]) + "\n  ... (all in report.txt)", flush=True)


def main():
    args = [a for a in sys.argv[1:] if a not in ("--save", "--report")]
    save, with_report = "--save" in sys.argv, "--report" in sys.argv
    if len(args) != 1:
        sys.exit(__doc__)
    t0 = time.time()
    d = load_db(args[0])
    print(f"Loaded {d.n:,} smORFs in {time.time() - t0:.1f}s", flush=True)
    combos = build_combinations(d)
    print(f"\nPART 1 - PREDICT: {len(combos)} rules on all smORFs", flush=True)
    preds, pcts = predict(d, combos, save)
    print("\nPART 2 - EVALUATE: annotated smORFs only", flush=True)
    rows = evaluate(d, combos, preds, pcts, save)
    if with_report:
        print("\nPART 3 - REPORT", flush=True)
        report(d, combos, rows)
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min."
          + (" Saved predictions_all_smorfs.tsv, combinations_results.csv, results_by_bin.csv" if save
             else " No files saved (add --save to save them).")
          + (" Saved report.txt, random_forest_predictions.tsv, one_factor_comparisons.tsv" if with_report else ""))


if __name__ == "__main__":
    main()
