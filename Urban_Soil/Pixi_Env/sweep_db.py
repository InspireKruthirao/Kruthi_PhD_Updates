#!/usr/bin/env python3
"""
Predict the COG category of every smORF from its neighbour genes with 192
rules, and score each rule on the smORFs whose category is known.

    python3 sweep_db.py smorf_neighbourhoods.sqlite           show the results
    python3 sweep_db.py smorf_neighbourhoods.sqlite --save    also save:
        predictions_all_smorfs.tsv   every smORF x every rule
        combinations_results.csv     one row per rule: accuracy and coverage
        results_by_bin.csv           the same, per bin

The input (from smorf_neighbourhoods_sqlite.py) has one row per smORF with
a JSON list of its occurrences: neighbour genes as [rank, [COG letters]],
neighbour_info (side, same strand, distance, gap, COG ID) and the smORF's
own COG letters (target_categories). C027/C028 are the original method.
To add a rule, add it at the end of build_combinations().
"""

import csv
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
TIERS = (("VERY_HIGH", 50), ("HIGH", 40), ("MEDIUM", 30), ("LOW", 0))
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

# JSON patterns, rewritten as objects so that polars can parse them
NEIGHBOUR = r"\[(\d+), \[([^\]]*)\]\]"
NEIGHBOUR_AS_OBJECT = r'{"rank":${1},"cats":[${2}]}'
OCCS_TYPE = pl.List(pl.Struct({
    "neighbours": pl.List(pl.Struct({"rank": pl.Int64, "cats": pl.List(pl.String)})),
    "target_categories": pl.List(pl.String),
}))
INFO = r'\["(up|down)", (\d+), (-?\d+), (-?\d+), ("[^"]*")\]'
INFO_AS_OBJECT = r'{"side":"${1}","same":${2},"dist":${3},"gap":${4},"cog":${5}}'
INFO_TYPE = pl.List(pl.Struct({"neighbour_info": pl.List(pl.Struct(
    {"side": pl.String, "same": pl.Int8, "dist": pl.Int64, "gap": pl.Int64, "cog": pl.String}))}))


# ---------------------------------------------------------------- data ----

@dataclass
class Data:
    """Everything loaded from the database. Categories are integer codes
    (indices into labels). nb has one entry per (neighbour gene, category):
    smorf, occ, ent (gene), rank, code, key (smorf * n_cats + code), ncat
    (categories of the gene) and ord (position of the gene by rank)."""
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

    def code(self, label):
        return self.labels.index(label)

    def excluded(self, codes, exclude):
        """Is the category of each code in exclude?"""
        return np.array([c in exclude for c in self.labels[:self.n_cats]])[codes]

    def left(self, exclude):
        """Per nb entry, with exclude dropped: kept?, categories left on the
        gene (k), and its index in vote_table (vote)."""
        if exclude not in self._left:
            nb = self.nb
            keep = ~self.excluded(nb["code"], exclude)
            starts = run_starts(nb["ent"])
            k = np.repeat(np.add.reduceat(keep.astype(np.int32), starts), np.diff(np.r_[starts, len(keep)]))
            vote = (nb["rank"].astype(np.int32) * (nb["ncat"].max() + 1) + k) * self.n_cats + nb["code"]
            self._left[exclude] = keep, k, vote
        return self._left[exclude]

    def voters(self, n, exclude, multi, columns):
        """The given columns (plus vote) of the nb entries that vote."""
        sel = (n, exclude, multi == "ignore", columns)
        if sel not in self._voters:
            keep, k, vote = self.left(exclude)
            keep = keep & (self.nb["rank"] <= n)
            if multi == "ignore":
                keep &= k == 1
            if len(self._voters) == 2:  # keep only the last two
                del self._voters[next(iter(self._voters))]
            self._voters[sel] = select({**self.nb, "vote": vote}, keep, columns + ("vote",))
        return dict(self._voters[sel])

    def counts(self, exclude):
        """The kept nb entries counted by (key, rank, k), so votes can be
        added up from counts."""
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


def parse(raw, pattern, as_object, dtype):
    """Parse the JSON of raw in 64 chunks, in parallel."""
    step = -(-raw.height // 64)
    return pl.concat([raw.slice(i, step).lazy().select(
        "smorf", pl.col("json").str.replace_all(pattern, as_object).str.json_decode(dtype).alias("o"))
        .explode("o", empty_as_null=True).drop_nulls("o")
        for i in range(0, raw.height, step)], parallel=True).collect()


def load_db(db_path):
    """Read the database once into a Data."""
    conn = sqlite3.connect(db_path)
    raw = pl.DataFrame(conn.execute(
        "SELECT bin_name, smorf_id, occurrences_json FROM smorfs ORDER BY bin_name, smorf_id"
    ).fetchall(), schema=["bin", "smorf_id", "json"], orient="row").with_row_index("smorf")
    conn.close()

    occs = parse(raw, NEIGHBOUR, NEIGHBOUR_AS_OBJECT, OCCS_TYPE).with_row_index("occ").unnest("o")
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

    # gene strands, distances and COG IDs (rules C151-C192)
    if raw.height and len(d.nb["ent"]):
        info = (parse(raw, INFO, INFO_AS_OBJECT, INFO_TYPE).with_row_index("occ").unnest("o")
                .explode("neighbour_info", empty_as_null=True).drop_nulls("neighbour_info")
                .unnest("neighbour_info"))
        if info.height > d.nb["ent"].max():
            x = {c: info[c].to_numpy() for c in ("occ", "same", "dist", "gap")}
            if np.array_equal(x["occ"][d.nb["ent"]], d.nb["occ"]):  # same genes, same order
                x["side"] = (info["side"] == "down").to_numpy().astype(np.int8)
                x["cog"] = info["cog"]
                d.info = x
    if d.info is None:
        print("  (no gene distances in this database: rules C151-C192 predict nothing)", flush=True)
    return d


# -------------------------------------------------------- vote helpers ----

def run_starts(a):
    """Where each run of equal values starts"""
    return np.flatnonzero(np.r_[True, a[1:] != a[:-1]]) if len(a) else np.zeros(0, dtype=int)


def resolve(d, key, v, n_groups):
    """Winner of each group from votes v at key (group * n_cats + code):
    (code, winner's % of the votes); UNRESOLVED on a tie, NO_EVIDENCE if no votes."""
    votes = np.bincount(key, weights=v, minlength=n_groups * d.n_cats)
    present = np.zeros(n_groups * d.n_cats, dtype=bool)
    present[key] = True
    keys = np.flatnonzero(present)
    return pick(d, keys, votes[keys], n_groups)


def pick(d, keys, x, n_groups):
    """As resolve, from totals x at sorted unique keys"""
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
    """columns of arrays where keep"""
    if keep.all():
        return {c: arrays[c] for c in columns}
    keep = np.flatnonzero(keep)
    return {c: arrays[c].take(keep) for c in columns}


def vote_table(d, weights, multi, background):
    """Vote by (rank, k, code): rank weight, / k if "split", / background."""
    v = np.array([np.nan] + [weights[r] for r in range(1, MAX_RANK + 1)])[:, None, None]
    ks = d.nb["ncat"].max() + 1
    if multi == "split":
        v = v / np.maximum(np.arange(ks), 1)[None, :, None]
    if background is not None:
        v = v / background[None, None, :]
    return np.broadcast_to(v, (MAX_RANK + 1, ks, d.n_cats)).ravel()


def tally(d, weights, n, exclude, multi, background, columns=("key",)):
    """The voting nb entries (columns) with their vote v"""
    x = d.voters(n, tuple(exclude), multi, columns)
    x["v"] = vote_table(d, weights, multi, background).take(x.pop("vote"))
    return x


def nothing(d):
    """NO_EVIDENCE for every smORF"""
    n = d.smorfs.height
    return np.full(n, d.code("NO_EVIDENCE")), np.zeros(n)


def low_confidence(d, pred, low):
    return np.where((pred < d.n_cats) & low, d.code("LOW_CONFIDENCE"), pred)


# --------------------------------------------------------- the methods ----
# each returns (predicted code, winner's %) per smORF

def weighted_vote(d, curve="linear", n=5, exclude=(), multi="split",
                  background=None, cutoff=None, min_occ=None):
    """Neighbours (rank <= n) vote for their categories, weighted by curve,
    over all occurrences. multi: "split", "full" or "ignore" multi-letter
    genes. cutoff: winner needs >= cutoff %. min_occ: votes from >= min_occ
    occurrences."""
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
        pred = low_confidence(d, pred, pct < cutoff - PCT_TOLERANCE)
    if min_occ:
        x = d.voters(n, tuple(exclude), multi, ("smorf", "occ"))
        few = np.bincount(x["smorf"][run_starts(x["occ"])], minlength=n_smorfs) < min_occ
        pred = np.where(few, d.code("TOO_FEW_OCCURRENCES"), pred)
        pct = np.where(few, 0.0, pct)
    return pred, pct


def majority(d, n=5, k=4, exclude=()):
    """Most common category among the n closest neighbours of all
    occurrences, if at least k/n of them have it (ties: first counted)."""
    nb = d.nb_by_rank
    keep = (nb["ord"] < n) & ~d.excluded(nb["code"], exclude)
    key, ent = select(nb, keep, ("key", "ent")).values()
    n_smorfs = d.smorfs.height
    total = np.bincount(key[run_starts(ent)] // d.n_cats, minlength=n_smorfs)
    counts = np.bincount(key, minlength=n_smorfs * d.n_cats).reshape(n_smorfs, d.n_cats)
    first = np.full(n_smorfs * d.n_cats, len(key))
    np.minimum.at(first, key, np.arange(len(key)))
    top = counts.max(axis=1)
    top_cat = np.where(counts == top[:, None], first.reshape(n_smorfs, d.n_cats), len(key)).argmin(axis=1)
    share = np.divide(top, total, out=np.zeros(n_smorfs), where=total > 0)
    agree = (total > 0) & (share >= k / n)
    return (np.where(agree, top_cat, d.code("NO_EVIDENCE")),
            np.where(agree, share * 100, 0.0))


def consensus(d, curve="linear", n=5, exclude=()):
    """Each occurrence predicts by weighted vote, then occurrences vote."""
    x = tally(d, WEIGHTS[curve], n, exclude, "split", None, ("occ", "code"))
    new = np.r_[True, x["occ"][1:] != x["occ"][:-1]][:len(x["occ"])]
    occ_pred, _ = resolve(d, (np.cumsum(new) - 1) * d.n_cats + x["code"], x["v"], new.sum())
    made = occ_pred < d.n_cats
    smorf = d.occ_smorf[x["occ"][new][made]]
    return resolve(d, smorf.astype(np.int64) * d.n_cats + occ_pred[made], np.ones(made.sum()),
                   d.smorfs.height)


def control_always(d, category, exclude=()):
    """Control: always the same category"""
    n = d.smorfs.height
    return np.full(n, d.code(category)), np.zeros(n)


def control_random(d, categories, freqs, rng, exclude=()):
    """Control: a random category, weighted by how common it is"""
    n = d.smorfs.height
    return np.array([d.code(c) for c in rng.choices(categories, freqs, k=n)]), np.zeros(n)


def control_shuffled(d, rng, exclude=()):
    """Control: the original method on another random smORF's neighbours"""
    n = d.smorfs.height
    other = np.array([rng.randrange(n) for _ in range(n)])
    pred, pct = weighted_vote(d)
    return pred[other], pct[other]


def presence(d, n, exclude):
    """smORF x category: in how many occurrences the category is within
    rank n; and the number of occurrences per smORF (cached, read-only)."""
    made = ("presence", n, tuple(exclude))
    if made not in d.cache:
        nb = d.nb
        keep = (nb["rank"] <= n) & ~d.excluded(nb["code"], exclude)
        pair = np.unique(nb["occ"][keep].astype(np.int64) * d.n_cats + nb["code"][keep])
        key = d.occ_smorf[pair // d.n_cats].astype(np.int64) * d.n_cats + pair % d.n_cats
        n_smorfs = d.smorfs.height
        count = np.bincount(key, minlength=n_smorfs * d.n_cats).reshape(n_smorfs, d.n_cats)
        d.cache[made] = count, np.bincount(d.occ_smorf, minlength=n_smorfs)
    return d.cache[made]


def top_of(d, score, eligible):
    """Highest-scoring eligible category per row: (code, column)"""
    s = np.where(eligible, score, -np.inf)
    top = np.where(eligible.any(axis=1), s.max(axis=1), 0.0)
    win = eligible & (np.abs(s - top[:, None]) < TIE_TOLERANCE)
    pred = np.where(win.sum(axis=1) == 1, s.argmax(axis=1), d.code("UNRESOLVED"))
    return np.where(eligible.any(axis=1), pred, d.code("NO_EVIDENCE")), s.argmax(axis=1)


def conserved(d, n=3, exclude=(), min_frac=None):
    """Category within rank n in the most occurrences (if in >= min_frac)."""
    count, n_occ = presence(d, n, exclude)
    pred, col = top_of(d, count, count > 0)
    frac = count[np.arange(len(pred)), col] / np.maximum(n_occ, 1)
    pct = np.where(pred < d.n_cats, frac * 100, 0.0)
    if min_frac:
        pred = low_confidence(d, pred, frac < min_frac - PCT_TOLERANCE)
    return pred, pct


def sandwiched(d, exclude=(), min_genomes=1, min_frac=None, same_strand=False, max_dist=None):
    """Closest gene on both sides shares a category, in >= min_genomes
    occurrences (and >= min_frac of them). Optionally both genes on the
    target's strand and/or within max_dist bp."""
    nb = d.nb
    keep = (nb["rank"] == 1) & ~d.excluded(nb["code"], exclude)
    if same_strand or max_dist is not None:
        if d.info is None:
            return nothing(d)
        ok = np.ones(len(d.info["occ"]), dtype=bool)
        if same_strand:
            ok &= d.info["same"] == 1
        if max_dist is not None:
            ok &= d.info["dist"] <= max_dist
        keep &= ok[nb["ent"]]
    pair, k = np.unique(nb["occ"][keep].astype(np.int64) * d.n_cats + nb["code"][keep],
                        return_counts=True)
    both = pair[k >= 2]  # on both rank-1 genes
    key = d.occ_smorf[both // d.n_cats].astype(np.int64) * d.n_cats + both % d.n_cats
    n_smorfs = d.smorfs.height
    count = np.bincount(key, minlength=n_smorfs * d.n_cats).reshape(n_smorfs, d.n_cats)
    n_occ = np.bincount(d.occ_smorf, minlength=n_smorfs)
    pred, col = top_of(d, count, count >= max(min_genomes, 1))
    frac = count[np.arange(n_smorfs), col] / np.maximum(n_occ, 1)
    if min_frac:
        pred = np.where((pred < d.n_cats) & (frac < min_frac - PCT_TOLERANCE), d.code("NO_EVIDENCE"), pred)
    return pred, np.where(pred < d.n_cats, frac * 100, 0.0)


def enriched(d, n=3, exclude=(), min_genomes=2):
    """Category most enriched within rank n compared with all occurrences."""
    count, n_occ = presence(d, n, exclude)
    frac = count / np.maximum(n_occ, 1)[:, None]
    background = count.sum(axis=0) / max(n_occ.sum(), 1)
    enrichment = np.divide(frac, background, out=np.zeros_like(frac), where=background > 0)
    pred, col = top_of(d, enrichment, count >= min_genomes)
    return pred, np.where(pred < d.n_cats, frac[np.arange(len(pred)), col] * 100, 0.0)


def conserved_vote(d, curve="very_steep", n=3, exclude=()):
    """Weighted vote times the share of occurrences that have the category."""
    x = tally(d, WEIGHTS[curve], n, exclude, "split", None)
    n_smorfs = d.smorfs.height
    votes = np.bincount(x["key"], weights=x["v"], minlength=n_smorfs * d.n_cats)
    count, n_occ = presence(d, n, exclude)
    scaled = votes * (count / np.maximum(n_occ, 1)[:, None]).ravel()
    keys = np.flatnonzero(count.ravel() > 0)
    return pick(d, keys, scaled[keys], n_smorfs)


def agree(d, first, second, exclude=()):
    """Predict only where two (method, settings) agree."""
    p1, s1 = first[0](d, **first[1])
    p2, _ = second[0](d, **second[1])
    ok = (p1 < d.n_cats) & (p1 == p2)
    return np.where(ok, p1, d.code("NO_EVIDENCE")), np.where(ok, s1, 0.0)


def in_operon(d, max_gap):
    """Per gene: it and every gene up to the target are on the target's
    strand with <= max_gap bp between them."""
    info = d.info
    fail = ~((info["same"] == 1) & (info["gap"] <= max_gap))
    run = np.r_[True, (info["occ"][1:] != info["occ"][:-1]) | (info["side"][1:] != info["side"][:-1])]
    failed = np.cumsum(fail)
    starts = np.flatnonzero(run)
    before = np.repeat(failed[starts] - fail[starts], np.diff(np.r_[starts, len(fail)]))
    return failed - before == 0


def vote_where(d, gene_ok, weight, n, exclude, cutoff=None):
    """Weighted vote (split) of the genes with gene_ok, rank <= n."""
    nb = d.nb
    keep, k, _ = d.left(tuple(exclude))
    use = keep & (nb["rank"] <= n) & gene_ok[nb["ent"]]
    pred, pct = resolve(d, nb["key"][use], weight[use] / k[use], d.smorfs.height)
    if cutoff:
        pred = low_confidence(d, pred, pct < cutoff - PCT_TOLERANCE)
    return pred, pct


def operon_vote(d, max_gap=100, curve="linear", n=5, exclude=(), cutoff=None):
    """Only genes in the target gene's operon vote."""
    if d.info is None:
        return nothing(d)
    weight = np.array([0.0] + [WEIGHTS[curve][r] for r in range(1, MAX_RANK + 1)])[d.nb["rank"]]
    return vote_where(d, in_operon(d, max_gap), weight, n, exclude, cutoff)


def distance_vote(d, half=500, n=5, exclude=()):
    """Vote halves every `half` bp from the target gene."""
    if d.info is None:
        return nothing(d)
    weight = 0.5 ** (np.maximum(d.info["dist"][d.nb["ent"]], 0) / half)
    return vote_where(d, np.ones(len(d.info["occ"]), dtype=bool), weight, n, exclude)


# ------------------------------------- similar neighbourhoods (Annoy, Voyager) ----

def context_vectors(d, max_rank=3, dims=512):
    """smORF x COG ID: share of occurrences with the COG ID within
    max_rank, for the dims most common COG IDs."""
    key = ("ctx", max_rank, dims)
    if key not in d.cache:
        if d.info is None:
            return None
        nb = d.nb
        first = run_starts(nb["ent"])  # one entry per gene
        use = first[nb["rank"][first] <= max_rank]
        cog = d.info["cog"].gather(nb["ent"][use]).to_numpy()
        occ = nb["occ"][use]
        ok = cog != ""
        ids, code = np.unique(cog[ok], return_inverse=True)
        pair = np.unique(occ[ok].astype(np.int64) * len(ids) + code)
        sm = d.occ_smorf[pair // len(ids)].astype(np.int64)
        c = pair % len(ids)
        top = np.argsort(-np.bincount(np.unique(sm * len(ids) + c) % len(ids), minlength=len(ids)))[:dims]
        col = np.full(len(ids), -1)
        col[top] = np.arange(len(top))
        keep = col[c] >= 0
        print(f"  C185-C192: using {len(top)} of {len(ids)} COG IDs, {100 * keep.mean():.1f}% of "
              f"(occurrence, COG ID) pairs", flush=True)
        n_smorfs = d.smorfs.height
        x = np.bincount(sm[keep] * len(top) + col[c[keep]],
                        minlength=n_smorfs * len(top)).reshape(n_smorfs, len(top)).astype(np.float32)
        x /= np.maximum(np.bincount(d.occ_smorf, minlength=n_smorfs), 1)[:, None]
        d.cache[key] = x
    return d.cache[key]


def target_genes(d):
    """(smorf, gene) pairs of target genes (sample|Original_ID), sorted;
    None if the database has no neighbours table."""
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
                         gene=pl.col("parts").list.get(1) + "|" + pl.col("original_id")))
            index = d.smorfs.select("smorf_id", smorf=pl.col("smorf").cast(pl.Int64))
            t = t.join(index, on="smorf_id", how="inner")
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


def nearest(d, engine, k, max_rank, dims, **opts):
    """For each smORF with a context vector (who): its k + 1 most similar
    smORFs (nn, -1 if none) and their cosine similarity (sim). Cached."""
    x = context_vectors(d, max_rank, dims)
    if x is None:
        return None
    key = (engine, k, max_rank, dims) + tuple(sorted(opts.items()))
    if key in d.cache:
        return d.cache[key]
    who = np.flatnonzero(x.any(axis=1))
    m = k + 1
    nn = np.full((len(who), m), -1, dtype=np.int64)
    dist = np.zeros((len(who), m))
    if engine == "annoy":
        try:
            from annoy import AnnoyIndex
        except ImportError:
            print("  (python-annoy not installed: rules C185-C188 predict nothing)", flush=True)
            return None
        ikey = ("annoy index", max_rank, dims, opts["n_trees"])
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
                dist[row, :len(got)] = dj

        with ThreadPoolExecutor(cpus()) as pool:  # same results as one thread
            list(pool.map(search, np.array_split(np.arange(len(who)), 16 * cpus())))
        sim = 1.0 - dist * dist / 2  # angular distance -> cosine
    elif engine == "voyager":
        try:
            from voyager import Index, Space
        except ImportError:
            print("  (voyager not installed: rules C189-C192 predict nothing)", flush=True)
            return None
        ikey = ("voyager index", max_rank, dims, opts["M"], opts["ef_construction"])
        if ikey not in d.cache and len(who):
            index = Index(Space.Cosine, num_dimensions=x.shape[1], M=opts["M"],
                          ef_construction=opts["ef_construction"], random_seed=1, max_elements=len(who))
            index.add_items(x[who], ids=who.tolist(), num_threads=1)  # one thread: same index every run
            d.cache[ikey] = index
        if len(who):
            m = min(m, len(who))
            got, dj = d.cache[ikey].query(x[who], k=m, query_ef=max(opts["query_ef"], m))
            nn[:, :m] = got.astype(np.int64)
            dist[:, :m] = dj.astype(np.float64)
        sim = 1.0 - dist  # cosine distance -> cosine
    else:
        raise ValueError(f"unknown search engine: {engine}")
    d.cache[key] = (who, nn, np.where(nn >= 0, sim, 0.0))
    return d.cache[key]


def similar_vote(d, engine="annoy", k=10, max_rank=3, dims=512, exclude=(), **opts):
    """The k smORFs with the most similar neighbourhoods vote for their known
    categories, weighted by similarity. Not the smORF itself, not smORFs on
    the same target gene, not smORFs with nothing in common."""
    defaults = {"annoy": {"n_trees": 20}, "voyager": {"M": 12, "ef_construction": 200, "query_ef": 100}}
    found = nearest(d, engine, k, max_rank, dims, **{**defaults.get(engine, {}), **opts})
    if found is None:
        return nothing(d)
    who, nn, sim = found

    other = (nn >= 0) & (nn != who[:, None])
    use = other & (np.cumsum(other, axis=1) <= k)  # first k others
    i = np.repeat(who, use.sum(axis=1)).astype(np.int64)
    j, s = nn[use], sim[use]

    known = d.targets.filter(~pl.col("cat").is_in(list(exclude)))
    known_smorf = known["smorf"].to_numpy().astype(np.int64)
    order = np.argsort(known_smorf, kind="stable")
    known_smorf = known_smorf[order]
    known_code = known["cat"].cast(pl.Enum(d.labels[:d.n_cats])).to_physical().to_numpy().astype(np.int64)[order]
    n_known = np.bincount(known_smorf, minlength=d.smorfs.height)
    first_known = np.r_[0, np.cumsum(n_known)][:-1]

    ok = (n_known[j] > 0) & (s > SIM_TOLERANCE)
    targets = target_genes(d)
    if targets is not None:
        ok &= ~share_gene(i, j, targets)
    i, j, s = i[ok], j[ok], s[ok]
    if not len(i):
        return nothing(d)

    n = n_known[j]  # one vote per category of the neighbour
    pos = np.repeat(first_known[j], n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    w = np.repeat(s / n, n)
    keys = np.repeat(i, n) * d.n_cats + known_code[pos]
    return resolve(d, keys, w, d.smorfs.height)


# ------------------------------------------------------------- the rules ----

def build_combinations(d):
    """All rules: id, group, rule (text), fn, settings, exclude. IDs follow
    the order, so add new rules at the end."""
    combos = []

    def add(group, rule, fn, **settings):
        combos.append(dict(id=f"C{len(combos) + 1:03d}", group=group, rule=rule, fn=fn,
                           settings=settings, exclude=settings.get("exclude", ())))

    cats_opts = (("all categories", ()), ("excluding S and R", SR))
    curve_name = {"flat": "flat", "linear": "linear", "steeper": "steeper", "very_steep": "very steep"}

    # A. weighted vote (34)
    for n in (1, 2, 3, 4, 5):
        for curve in ("linear", "steeper", "flat", "very_steep"):
            if n == 1 and curve != "linear":
                continue  # same answer for every weighting
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

    # D. down-weight common categories (6)
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

    # H. multi-letter categories (4)
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

    # K. conserved neighbourhood (24)
    for n in (1, 2, 3, 5):
        for frac in (None, 0.5, 0.7):
            for cname, ex in cats_opts:
                need = f", only if in >= {frac:.0%} of its genomes" if frac else ""
                add("K. Conserved neighbourhood", f"category found within {n} neighbour(s) in the most "
                    f"genomes{need}, {cname}", conserved, n=n, exclude=ex, min_frac=frac)

    # L. sandwiched (6)
    for genomes, frac, text in ((1, None, "in at least 1 genome"), (2, None, "in at least 2 genomes"),
                                (1, 0.5, "in at least 50% of its genomes")):
        for cname, ex in cats_opts:
            add("L. Sandwiched", f"closest gene on both sides has the same category, {text}, {cname}",
                sandwiched, exclude=ex, min_genomes=genomes, min_frac=frac)

    # M. two methods agree (4)
    for cname, ex in cats_opts:
        add("M. Two methods agree", f"very steep weights, 3 neighbours, and conserved neighbourhood "
            f"(3 neighbours) predict the same category, {cname}", agree, exclude=ex,
            first=(weighted_vote, dict(curve="very_steep", n=3, exclude=ex)),
            second=(conserved, dict(n=3, exclude=ex)))
    for cname, ex in cats_opts:
        add("M. Two methods agree", f"original method and 'at least 3 of the 5 nearest agree' predict "
            f"the same category, {cname}", agree, exclude=ex,
            first=(weighted_vote, dict(curve="linear", n=5, exclude=ex)),
            second=(majority, dict(n=5, k=3, exclude=ex)))

    # N. enrichment (4)
    for n in (3, 5):
        for cname, ex in cats_opts:
            add("N. Enrichment over background", f"category most enriched within {n} neighbours compared "
                f"with the whole data, found in >= 2 genomes, {cname}", enriched, n=n, exclude=ex)

    # O. conservation-weighted vote (4)
    for curve, n in (("very_steep", 3), ("linear", 5)):
        for cname, ex in cats_opts:
            add("O. Conservation-weighted vote", f"{curve_name[curve]} weights, {n} neighbours, each "
                f"category's vote times the share of genomes that have it, {cname}",
                conserved_vote, curve=curve, n=n, exclude=ex)

    # P. very steep cutoffs (8)
    for cut in (50, 60, 70, 80):
        for cname, ex in cats_opts:
            add("P. Confidence cutoff (very steep)", f"very steep weights, 3 neighbours, predict only if "
                f"winner has >= {cut}% of the vote, {cname}", weighted_vote,
                curve="very_steep", n=3, exclude=ex, cutoff=cut)

    # Q. operon vote (12)
    for gap in (50, 100, 200):
        for curve in ("flat", "linear"):
            for cname, ex in cats_opts:
                add("Q. Operon vote", f"{curve} weights, only neighbours in the target gene's operon (same "
                    f"strand, <= {gap} bp between genes), {cname}", operon_vote,
                    max_gap=gap, curve=curve, n=5, exclude=ex)

    # R. sandwiched with strand/distance (12)
    for strand, dist, text in ((True, None, "both on the target's strand"),
                               (False, 100, "both within 100 bp of the target gene"),
                               (True, 100, "both on the target's strand and within 100 bp")):
        for genomes, frac, need in ((2, None, "in at least 2 genomes"), (1, 0.5, "in at least 50% of its genomes")):
            for cname, ex in cats_opts:
                add("R. Sandwiched (strand/distance)", f"closest gene on both sides has the same category, "
                    f"{text}, {need}, {cname}", sandwiched, same_strand=strand, max_dist=dist,
                    exclude=ex, min_genomes=genomes, min_frac=frac)

    # S. distance-weighted vote (6)
    for half in (250, 500, 1000):
        for cname, ex in cats_opts:
            add("S. Distance-weighted vote", f"5 neighbours, vote halves every {half} bp from the target "
                f"gene, {cname}", distance_vote, half=half, n=5, exclude=ex)

    # T. operon vote + cutoff (4)
    for cut in (50, 70):
        for cname, ex in cats_opts:
            add("T. Operon vote + cutoff", f"linear weights, operon neighbours (<= 100 bp), predict only if "
                f"winner has >= {cut}% of the vote, {cname}", operon_vote,
                max_gap=100, curve="linear", n=5, exclude=ex, cutoff=cut)

    # U, V. similar neighbourhoods: Annoy (4), Voyager (4)
    for group, engine in (("U. Annoy nearest neighbours", "annoy"), ("V. Voyager nearest neighbours", "voyager")):
        name = "Annoy" if engine == "annoy" else "Voyager"
        for k in (10, 25):
            for cname, ex in cats_opts:
                add(group, f"the {k} smORFs with the most similar neighbourhoods "
                    f"(COG IDs within 3 genes, {name}) vote, {cname}", similar_vote, engine=engine, k=k,
                    exclude=ex)
    return combos


# ------------------------------------------------------------- predict ----

def predict(d, combos, save=True):
    """Every rule predicts for every smORF."""
    preds, pcts = {}, {}
    for c in combos:
        t = time.time()
        preds[c["id"]], pcts[c["id"]] = p, s = c["fn"](d, **c["settings"])
        made = (p < d.n_cats).sum()
        print(f"  {c['id']}  predicted {made / d.smorfs.height * 100:5.1f}% of smORFs  "
              f"({time.time() - t:.1f}s)  {c['rule']}", flush=True)
    if save:
        labels = pl.Series(d.labels)
        (d.smorfs.select("bin", "smorf_id",
                         annotated=pl.when("annotated").then(pl.lit("yes")).otherwise(pl.lit("no")),
                         known_category="known_category")
         .with_columns(labels.gather(preds[c["id"]]).alias(c["id"]) for c in combos)
         .write_csv("predictions_all_smorfs.tsv", separator="\t", line_terminator="\r\n",
                    quote_style="never"))
    return preds, pcts


# ------------------------------------------------------------ evaluate ----

def pct(a, b):
    return f"{a / b * 100:.1f}" if b else ""


def evaluate(d, combos, preds, pcts, save=True):
    """Score each rule on the smORFs with a known category (not excluded by
    the rule): correct if any known category. Coverage = % of annotated
    smORFs evaluated."""
    bins = sorted(set(d.smorfs["bin"].to_list()))
    bin_of = d.smorfs["bin"].cast(pl.Enum(bins)).to_physical().to_numpy()
    annotated = d.smorfs["annotated"].to_numpy()
    annot_by_bin = np.bincount(bin_of[annotated], minlength=len(bins))
    smorf = d.targets["smorf"].to_numpy()
    cat = d.targets["cat"].cast(pl.Enum(d.labels[:d.n_cats])).to_physical().to_numpy()
    rows, by_bin = [], []
    known = {}  # exclude -> (smorf x category known, smorf has one)

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
            tier = np.where(~(s < low - PCT_TOLERANCE), i, tier)

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
    if save:
        for name, out, fields in (("combinations_results.csv", rows, cols),
                                  ("results_by_bin.csv", by_bin, list(by_bin[0]))):
            with open(name, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields)
                w.writeheader()
                w.writerows(out)

    print(f"\nAnnotated smORFs used for evaluation: {d.smorfs['annotated'].sum():,} of {d.smorfs.height:,}")
    print("\nAll rules, most accurate first:")
    print(f"  {'rank':>4}  {'ID':4}  {'accuracy':>8}  {'coverage':>8}  {'evaluated':>9}  group: rule")
    for r in ranked:
        print(f"  {r['Accuracy rank']:>4}  {r['ID']}  {r['Accuracy (%)']:>7}%  "
              f"{r['Coverage (% of annotated smORFs evaluated)']:>7}%  {r['Evaluated (n)']:>9,}  "
              f"{r['Group']}: {r['Rule']}")
    orig = next(r for r in rows if "ORIGINAL" in r["Rule"] and "all categories" in r["Rule"])
    print(f"\nOriginal method ({orig['ID']}): accuracy {orig['Accuracy (%)']}%  "
          f"coverage {orig['Coverage (% of annotated smORFs evaluated)']}%")


def main():
    args = [a for a in sys.argv[1:] if a != "--save"]
    save = "--save" in sys.argv[1:]
    if len(args) != 1:
        sys.exit(__doc__)
    t0 = time.time()
    d = load_db(args[0])
    print(f"Loaded {d.smorfs.height:,} smORFs in {time.time() - t0:.1f}s", flush=True)
    combos = build_combinations(d)

    print(f"\nPART 1 - PREDICT: {len(combos)} rules on all smORFs", flush=True)
    preds, pcts = predict(d, combos, save)
    print("\nPART 2 - EVALUATE: annotated smORFs only", flush=True)
    evaluate(d, combos, preds, pcts, save)

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min."
          + (" Saved predictions_all_smorfs.tsv, combinations_results.csv, results_by_bin.csv" if save
             else " No files saved (add --save to save them)."))


if __name__ == "__main__":
    main()
