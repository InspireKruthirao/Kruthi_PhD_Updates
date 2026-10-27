#!/usr/bin/env python3
"""
COG function prediction for smORFs: 100 rule combinations, one file.

PART 1 - PREDICT   every combination predicts a COG category for ALL smORFs
                   (annotated or not)            -> predictions_all_smorfs.tsv
PART 2 - EVALUATE  each combination is scored ONLY on smORFs that already
                   have a known COG category     -> combinations_results.csv
                                                    results_by_bin.csv

    python3 sweep_db.py [db]          show the results on screen, save nothing
    python3 sweep_db.py [db] --save   also save the files below
                                      (default db: smorfs_all_bins.sqlite)

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
PCT_TOLERANCE = 1e-9  # scores are not rounded: exactly 50% may come out as 49.99999999999999
SIM_TOLERANCE = 1e-6  # similarity 0 may come out as 3e-08 (Annoy's distances are rounded)

MAX_RANK = 10
WEIGHTS = {  # weight of the neighbour at each distance (rank 1 = closest)
    "flat":       {r: 1.0 for r in range(1, MAX_RANK + 1)},
    "linear":     {r: max(0.0, 1.0 - 0.2 * (r - 1)) for r in range(1, MAX_RANK + 1)},  # 1, .8, .6, .4, .2
    "steeper":    {r: 0.5 ** (r - 1) for r in range(1, MAX_RANK + 1)},                # 1, .5, .25, ...
    "very_steep": {r: 0.25 ** (r - 1) for r in range(1, MAX_RANK + 1)},               # 1, .25, .06, ...
}


# =========================================================== load once ==

def read_smorfs(db_path):
    """The smorfs table: bin, smorf_id and the JSON of the occurrences, in
    (bin, smorf_id) order"""
    conn = sqlite3.connect(db_path)
    raw = pl.DataFrame(conn.execute(
        "SELECT bin_name, smorf_id, occurrences_json FROM smorfs ORDER BY bin_name, smorf_id"
    ).fetchall(), schema=["bin", "smorf_id", "json"], orient="row")
    conn.close()
    return raw


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


def load_db(db_path, raw=None):
    """Read every smORF from the database at db_path into a Data.

    smORFs are numbered (smorf) in (bin, smorf_id) order, occurrences (occ)
    and neighbour genes (ent) in file order. Neighbours without a category
    are dropped from nb, but still count for ord. A smORF is annotated if
    any of its occurrences has target categories (known_category lists
    them all, sorted and comma-separated). raw is the smorfs table if it
    was already read (see read_smorfs), so that it is read only once."""
    if raw is None:
        raw = read_smorfs(db_path)
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


# =========================================== methods added for C101-C150 ==

def presence(d, n, exclude):
    """For each smORF and category: in how many of the smORF's occurrences
    (genomes) the category is found among the neighbours with rank <= n,
    ignoring the categories in exclude. Returns (count, n_occ): count is a
    smORF x category array, n_occ the number of occurrences of each smORF.
    Made once for each (n, exclude) and shared by the rules (do not modify).
    """
    if not hasattr(d, "_presence"):
        d._presence = {}
    made = (n, tuple(exclude))
    if made in d._presence:
        return d._presence[made]
    nb = d.nb
    keep = (nb["rank"] <= n) & ~d.excluded(nb["code"], exclude)
    pair = np.unique(nb["occ"][keep].astype(np.int64) * d.n_cats + nb["code"][keep])
    key = d.occ_smorf[pair // d.n_cats].astype(np.int64) * d.n_cats + pair % d.n_cats
    n_smorfs = d.smorfs.height
    count = np.bincount(key, minlength=n_smorfs * d.n_cats).reshape(n_smorfs, d.n_cats)
    d._presence[made] = count, np.bincount(d.occ_smorf, minlength=n_smorfs)
    return d._presence[made]


def top_of(d, score, eligible):
    """The category with the highest score in each row of score (a smORF x
    category array), among the eligible entries: UNRESOLVED on a tie,
    NO_EVIDENCE if no entry is eligible. Returns the codes and the row
    index of the winner (for looking up other values)."""
    s = np.where(eligible, score, -np.inf)
    top = np.where(eligible.any(axis=1), s.max(axis=1), 0.0)
    win = eligible & (np.abs(s - top[:, None]) < TIE_TOLERANCE)
    pred = np.where(win.sum(axis=1) == 1, s.argmax(axis=1), d.code("UNRESOLVED"))
    return np.where(eligible.any(axis=1), pred, d.code("NO_EVIDENCE")), s.argmax(axis=1)


def conserved(d, n=3, exclude=(), min_frac=None):
    """Conserved neighbourhood: predict the category found near the smORF
    (rank <= n) in the most genomes. pct is the share of the smORF's
    genomes that have it. With min_frac, predict LOW_CONFIDENCE if fewer
    than that share of genomes have the winner."""
    count, n_occ = presence(d, n, exclude)
    pred, col = top_of(d, count, count > 0)
    frac = count[np.arange(len(pred)), col] / np.maximum(n_occ, 1)
    pct = np.where(pred < d.n_cats, frac * 100, 0.0)
    if min_frac:
        pred = np.where((pred < d.n_cats) & (frac < min_frac - PCT_TOLERANCE), d.code("LOW_CONFIDENCE"), pred)
    return pred, pct


def sandwiched(d, exclude=(), min_genomes=1, min_frac=None, same_strand=False, max_dist=None):
    """The smORF sits between two genes of the same category: in a genome,
    the closest gene on each side (both rank 1) share a category. Predict
    the category sandwiching the smORF in the most genomes, if in at least
    min_genomes genomes (and, with min_frac, that share of its genomes).
    pct is the share of the smORF's genomes where it is sandwiched.

    With same_strand and/or max_dist (rules C163-C174), the closest gene on
    each side only counts if it is on the target's strand and/or at most
    max_dist bp from the target gene (needs the gene distances, see
    gene_info)."""
    nb = d.nb
    keep = (nb["rank"] == 1) & ~d.excluded(nb["code"], exclude)
    if same_strand or max_dist is not None:
        info = gene_info(d)
        if info is None:
            return nothing(d)
        ok = np.ones(len(info["occ"]), dtype=bool)
        if same_strand:
            ok &= info["same"] == 1
        if max_dist is not None:
            ok &= info["dist"] <= max_dist
        keep &= ok[nb["ent"]]
    pair, k = np.unique(nb["occ"][keep].astype(np.int64) * d.n_cats + nb["code"][keep],
                        return_counts=True)
    both = pair[k >= 2]  # the category is on both rank-1 genes of this occurrence
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
    """Enrichment over background: a category's share of the smORF's
    genomes that have it near the smORF (rank <= n), divided by its share
    over all occurrences in the data. Predict the most enriched category
    among those found in at least min_genomes genomes. pct is the share of
    the smORF's genomes that have the winner."""
    count, n_occ = presence(d, n, exclude)
    frac = count / np.maximum(n_occ, 1)[:, None]
    background = count.sum(axis=0) / max(n_occ.sum(), 1)
    enrichment = np.divide(frac, background, out=np.zeros_like(frac), where=background > 0)
    pred, col = top_of(d, enrichment, count >= min_genomes)
    return pred, np.where(pred < d.n_cats, frac[np.arange(len(pred)), col] * 100, 0.0)


def conserved_vote(d, curve="very_steep", n=3, exclude=()):
    """Weighted vote (as weighted_vote) where each category's total vote
    is multiplied by the share of the smORF's genomes that have it near
    the smORF (as conserved): categories seen in only a few genomes count
    less. pct is the winner's share of these votes."""
    x = tally(d, WEIGHTS[curve], n, exclude, "split", None)
    n_smorfs = d.smorfs.height
    votes = np.bincount(x["key"], weights=x["v"], minlength=n_smorfs * d.n_cats)
    count, n_occ = presence(d, n, exclude)
    scaled = votes * (count / np.maximum(n_occ, 1)[:, None]).ravel()
    keys = np.flatnonzero(count.ravel() > 0)
    return pick(d, keys, scaled[keys], n_smorfs)


def agree(d, first, second, exclude=()):
    """Two methods must agree: predict the category of first only where
    second predicts the same one, else NO_EVIDENCE. first and second are
    (method, settings) pairs. pct is that of first."""
    p1, s1 = first[0](d, **first[1])
    p2, _ = second[0](d, **second[1])
    ok = (p1 < d.n_cats) & (p1 == p2)
    return np.where(ok, p1, d.code("NO_EVIDENCE")), np.where(ok, s1, 0.0)


# =========================================== methods added for C151-C184 ==
# These use the gene strands and distances stored in neighbour_info by
# smorf_neighbourhoods_sqlite.py (smorf_neighbourhoods.sqlite). For each
# neighbour gene: side (up/down), same_strand (as the target gene), dist_bp
# (bp between it and the target gene), gap_bp (bp between it and the
# previous gene towards the target) and its COG ID (used by C185-C192). The
# target gene is the gene the smORF overlaps or, for an intergenic smORF,
# the next gene.

INFO = r'\["(up|down)", (\d+), (-?\d+), (-?\d+), ("[^"]*")\]'
INFO_AS_OBJECT = r'{"side":"${1}","same":${2},"dist":${3},"gap":${4},"cog":${5}}'
INFO_TYPE = pl.List(pl.Struct({"neighbour_info": pl.List(pl.Struct(
    {"side": pl.String, "same": pl.Int8, "dist": pl.Int64, "gap": pl.Int64, "cog": pl.String}))}))

_load_db_without_info = load_db


def load_db(db_path):
    """load_db, reading the database only once: the same JSON also gives
    the gene strands, distances and COG IDs (see gene_info). Remembers
    db_path (for target_genes)."""
    raw = read_smorfs(db_path)
    d = _load_db_without_info(db_path, raw)
    d.db_path = db_path
    d._info = read_info(d, raw)
    if d._info is None:
        print("  (no gene distances in this database: rules C151-C192 predict nothing; "
              "use smorf_neighbourhoods.sqlite)", flush=True)
    return d


def read_info(d, raw):
    """Arrays indexed by neighbour gene (ent, as in d.nb): occ, side (0 up,
    1 down), same (1 if on the target's strand), dist and gap (bp), and cog
    (the COG ID, a polars Series, "" if none). None if raw has no
    neighbour_info (smorfs_all_bins.sqlite)."""
    if not raw.height:
        return None
    step = -(-raw.height // 64)
    info = (pl.concat([raw.slice(i, step).lazy().select(
        pl.col("json").str.replace_all(INFO, INFO_AS_OBJECT).str.json_decode(INFO_TYPE).alias("o"))
        .explode("o", empty_as_null=True).drop_nulls("o")
        for i in range(0, raw.height, step)], parallel=True)
            .collect().with_row_index("occ").unnest("o")
            .explode("neighbour_info", empty_as_null=True).drop_nulls("neighbour_info")
            .unnest("neighbour_info"))
    if not len(d.nb["ent"]) or info.height <= d.nb["ent"].max():
        return None
    x = {c: info[c].to_numpy() for c in ("occ", "same", "dist", "gap")}
    if not np.array_equal(x["occ"][d.nb["ent"]], d.nb["occ"]):  # same genes, same order
        return None
    x["side"] = (info["side"] == "down").to_numpy().astype(np.int8)
    x["cog"] = info["cog"]
    return x


def gene_info(d):
    """The gene strands, distances and COG IDs read by load_db (see
    read_info), or None if the database has none"""
    return getattr(d, "_info", None)


def nothing(d):
    """NO_EVIDENCE for every smORF (for databases without distances)"""
    n = d.smorfs.height
    return np.full(n, d.code("NO_EVIDENCE")), np.zeros(n)


def in_operon(d, max_gap):
    """For each neighbour gene: is it in the same operon as the target
    gene? That is, it and every gene between it and the target are on the
    target's strand, with at most max_gap bp between consecutive genes."""
    info = gene_info(d)
    fail = ~((info["same"] == 1) & (info["gap"] <= max_gap))
    run = np.r_[True, (info["occ"][1:] != info["occ"][:-1]) | (info["side"][1:] != info["side"][:-1])]
    failed = np.cumsum(fail)  # genes are in rank order within each (occurrence, side)
    starts = np.flatnonzero(run)
    before = np.repeat(failed[starts] - fail[starts], np.diff(np.r_[starts, len(fail)]))
    return failed - before == 0


def vote_where(d, gene_ok, weight, n, exclude, cutoff=None):
    """Weighted vote (as weighted_vote, "split") over the neighbour genes
    with gene_ok[ent], rank <= n; weight is the vote of each nb entry
    before splitting between the gene's categories."""
    nb = d.nb
    keep, k, _ = d.left(tuple(exclude))
    use = keep & (nb["rank"] <= n) & gene_ok[nb["ent"]]
    pred, pct = resolve(d, nb["key"][use], weight[use] / k[use], d.smorfs.height)
    if cutoff:
        pred = np.where((pred < d.n_cats) & (pct < cutoff - PCT_TOLERANCE), d.code("LOW_CONFIDENCE"), pred)
    return pred, pct


def operon_vote(d, max_gap=100, curve="linear", n=5, exclude=(), cutoff=None):
    """Only the neighbours in the target gene's operon vote (see
    in_operon), closer ones more (curve). With cutoff, predict
    LOW_CONFIDENCE if the winner has < cutoff % of the vote."""
    if gene_info(d) is None:
        return nothing(d)
    weight = np.array([0.0] + [WEIGHTS[curve][r] for r in range(1, MAX_RANK + 1)])[d.nb["rank"]]
    return vote_where(d, in_operon(d, max_gap), weight, n, exclude, cutoff)


def distance_vote(d, half=500, n=5, exclude=()):
    """Neighbours vote with a weight that halves every `half` bp of
    distance from the target gene (by bp, not by rank)."""
    info = gene_info(d)
    if info is None:
        return nothing(d)
    weight = 0.5 ** (np.maximum(info["dist"][d.nb["ent"]], 0) / half)
    return vote_where(d, np.ones(len(info["occ"]), dtype=bool), weight, n, exclude)


# =========================================== methods added for C185-C192 ==
# Nearest-neighbour vote: each smORF is described by which neighbour genes
# (COG IDs, e.g. COG0146) sit near it, as the share of its occurrences that
# have each one. A nearest-neighbour search finds the smORFs with the most
# similar neighbourhoods; their known categories vote. Two search engines:
# Annoy (github.com/spotify/annoy, trees; C185-C188) and Voyager
# (github.com/spotify/voyager, Spotify's successor to Annoy, a graph (HNSW),
# usually more accurate; new smORFs can be added to an index later; C189-
# C192). Needs smorf_neighbourhoods.sqlite (COG IDs) and the python-annoy
# and voyager packages.

def context_vectors(d, max_rank=3, dims=512):
    """smORF x COG ID matrix: the share of the smORF's occurrences with that
    COG ID among the neighbours with rank <= max_rank. Uses the dims COG IDs
    found near the most smORFs. None if the database has no COG IDs."""
    key = (max_rank, dims)
    if not hasattr(d, "_ctx"):
        d._ctx = {}
    if key not in d._ctx:
        info = gene_info(d)
        if info is None:
            return None
        nb = d.nb
        first = run_starts(nb["ent"])  # one entry per gene (genes with several letters repeat)
        use = first[nb["rank"][first] <= max_rank]
        cog = info["cog"].gather(nb["ent"][use]).to_numpy()
        occ = nb["occ"][use]
        ok = cog != ""
        ids, code = np.unique(cog[ok], return_inverse=True)
        pair = np.unique(occ[ok].astype(np.int64) * len(ids) + code)  # (occurrence, COG ID) once
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
        d._ctx[key] = x
    return d._ctx[key]


def target_genes(d):
    """The target genes (sample, Original_ID) of each smORF's occurrences,
    as (smorf, gene) pairs: two sorted int64 arrays, without duplicates.
    Read from the neighbours table; None if the database has none (a
    filtered copy), and then smORFs on the same gene are not excluded."""
    if not hasattr(d, "_targets"):
        d._targets = None
        conn = sqlite3.connect(d.db_path)
        try:
            rows = conn.execute("SELECT occ_id, original_id FROM neighbours WHERE is_target = 1").fetchall()
        except sqlite3.Error:
            rows = []
            print("  (no gene tables: smORFs sharing a target gene are not excluded in C185-C192)", flush=True)
        conn.close()
        if rows:
            # occ_id is smorf_id|sample|contig|start (see smorf_neighbourhoods_sqlite.py)
            t = (pl.DataFrame(rows, schema=["occ_id", "original_id"], orient="row")
                 .with_columns(parts=pl.col("occ_id").str.split("|"))
                 .select(smorf_id=pl.col("parts").list.get(0),
                         gene=pl.col("parts").list.get(1) + "|" + pl.col("original_id")))
            index = d.smorfs.select("smorf_id", smorf=pl.col("smorf").cast(pl.Int64))
            t = t.join(index, on="smorf_id", how="inner")
            gene = t["gene"].cast(pl.Categorical).to_physical().to_numpy().astype(np.int64)
            pair = np.unique(t["smorf"].to_numpy() * (gene.max() + 1) + gene)
            d._targets = (pair // (gene.max() + 1), pair % (gene.max() + 1), gene.max() + 1)
    return d._targets


def cpus():
    """The number of CPUs this job may use"""
    try:
        return max(len(os.sched_getaffinity(0)), 1)
    except AttributeError:
        return os.cpu_count() or 1


def nearest(d, engine, k, max_rank, dims, **opts):
    """For each smORF with a context vector (who, ascending): the k + 1
    smORFs with the most similar neighbourhoods (nn, -1 where fewer were
    found) and their cosine similarity (sim), in order, from the search
    engine "annoy" or "voyager". The index and the search are made once and
    shared by the rules that differ only in the categories they exclude.
    None if the engine is not installed or there are no COG IDs.

    Each engine is searched again for each k: its answer for the first k + 1
    depends on how many are asked for. Annoy searches on all CPUs (with
    the same results as on one)."""
    x = context_vectors(d, max_rank, dims)
    if x is None:
        return None
    key = (engine, k, max_rank, dims) + tuple(sorted(opts.items()))
    if key in d._ctx:
        return d._ctx[key]
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
        if ikey not in d._ctx:
            index = AnnoyIndex(x.shape[1], "angular")
            index.set_seed(0)
            for i in who:
                index.add_item(int(i), x[i])
            index.build(opts["n_trees"], n_jobs=1)
            d._ctx[ikey] = index
        index = d._ctx[ikey]

        def search(rows):  # (Annoy lets several searches run at the same time)
            for row in rows:
                got, dj = index.get_nns_by_item(int(who[row]), m, include_distances=True)
                nn[row, :len(got)] = got
                dist[row, :len(got)] = dj

        with ThreadPoolExecutor(cpus()) as pool:
            list(pool.map(search, np.array_split(np.arange(len(who)), 16 * cpus())))
        sim = 1.0 - dist * dist / 2  # angular distance -> cosine similarity
    elif engine == "voyager":
        try:
            from voyager import Index, Space
        except ImportError:
            print("  (voyager not installed: rules C189-C192 predict nothing)", flush=True)
            return None
        ikey = ("voyager index", max_rank, dims, opts["M"], opts["ef_construction"])
        if ikey not in d._ctx and len(who):
            index = Index(Space.Cosine, num_dimensions=x.shape[1], M=opts["M"],
                          ef_construction=opts["ef_construction"], random_seed=1, max_elements=len(who))
            index.add_items(x[who], ids=who.tolist(), num_threads=1)  # (one thread: the same index every run)
            d._ctx[ikey] = index
        if len(who):
            m = min(m, len(who))
            got, dj = d._ctx[ikey].query(x[who], k=m, query_ef=max(opts["query_ef"], m))
            nn[:, :m] = got.astype(np.int64)
            dist[:, :m] = dj.astype(np.float64)
        sim = 1.0 - dist  # Voyager's cosine distance is 1 - cosine similarity
    else:
        raise ValueError(f"unknown search engine: {engine}")
    d._ctx[key] = (who, nn, np.where(nn >= 0, sim, 0.0))
    return d._ctx[key]


def similar_vote(d, engine="annoy", k=10, max_rank=3, dims=512, exclude=(), **opts):
    """Predict the category that the k smORFs with the most similar
    neighbourhoods (see nearest; cosine similarity of context_vectors) most
    often have; each one's vote is its cosine similarity, split between its
    known categories. A smORF with similarity 0 (nothing in common; up to
    SIM_TOLERANCE, as Annoy's distances are rounded) does not vote, so a smORF whose neighbours share nothing with it gets
    NO_EVIDENCE. The smORF itself is never its own neighbour, and smORFs
    that share a target gene (see target_genes) do not vote for each other.
    pct is the winner's share of the votes.

    opts are the settings of the engine: n_trees for Annoy; M,
    ef_construction and query_ef (how many candidates each search looks
    through) for Voyager."""
    defaults = {"annoy": {"n_trees": 20}, "voyager": {"M": 12, "ef_construction": 200, "query_ef": 100}}
    found = nearest(d, engine, k, max_rank, dims, **{**defaults.get(engine, {}), **opts})
    if found is None:
        return nothing(d)
    who, nn, sim = found

    # the first k neighbours other than the smORF itself, in order (as the search returns them)
    other = (nn >= 0) & (nn != who[:, None])
    use = other & (np.cumsum(other, axis=1) <= k)
    i = np.repeat(who, use.sum(axis=1)).astype(np.int64)
    j, s = nn[use], sim[use]

    # the known categories of each smORF, in the order of d.targets
    known = d.targets.filter(~pl.col("cat").is_in(list(exclude)))
    known_smorf = known["smorf"].to_numpy().astype(np.int64)
    order = np.argsort(known_smorf, kind="stable")
    known_smorf = known_smorf[order]
    known_code = known["cat"].cast(pl.Enum(d.labels[:d.n_cats])).to_physical().to_numpy().astype(np.int64)[order]
    n_known = np.bincount(known_smorf, minlength=d.smorfs.height)
    first_known = np.r_[0, np.cumsum(n_known)][:-1]

    ok = (n_known[j] > 0) & (s > SIM_TOLERANCE)  # a neighbour with nothing in common does not vote
    targets = target_genes(d)
    if targets is not None:
        ok &= ~share_gene(i, j, targets)
    i, j, s = i[ok], j[ok], s[ok]
    if not len(i):
        return nothing(d)

    # one vote per (pair, category of the neighbour), each worth its share of the similarity
    n = n_known[j]
    pos = np.repeat(first_known[j], n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    w = np.repeat(s / n, n)
    keys = np.repeat(i, n) * d.n_cats + known_code[pos]
    return resolve(d, keys, w, d.smorfs.height)


def share_gene(i, j, targets, block=1_000_000):
    """For each pair (i[p], j[p]) of smORFs: do they have a target gene in
    common? targets is as returned by target_genes."""
    smorf, gene, n_genes = targets
    key = smorf * n_genes + gene  # sorted
    count = np.bincount(smorf, minlength=max(int(smorf.max()) + 1, int(i.max()) + 1, int(j.max()) + 1))
    start = np.r_[0, np.cumsum(count)][:-1]
    out = np.zeros(len(i), dtype=bool)
    for b in range(0, len(i), block):  # in blocks, to limit memory
        ii, jj = i[b:b + block], j[b:b + block]
        n = count[ii]
        pair = np.repeat(np.arange(len(ii)), n)
        g = gene[np.repeat(start[ii], n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))]
        want = jj[pair] * n_genes + g
        hit = key[np.minimum(np.searchsorted(key, want), len(key) - 1)] == want
        out[b:b + block] = np.bincount(pair[hit], minlength=len(ii)) > 0
    return out


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

    # ---- added rules C101-C150 (the rules above keep their IDs) ----

    # K. conserved neighbourhood (24)
    for n in (1, 2, 3, 5):
        for frac in (None, 0.5, 0.7):
            for cname, ex in cats_opts:
                need = f", only if in >= {frac:.0%} of its genomes" if frac else ""
                add("K. Conserved neighbourhood", f"category found within {n} neighbour(s) in the most "
                    f"genomes{need}, {cname}", conserved, n=n, exclude=ex, min_frac=frac)

    # L. sandwiched between genes of the same category (6)
    for genomes, frac, text in ((1, None, "in at least 1 genome"), (2, None, "in at least 2 genomes"),
                                (1, 0.5, "in at least 50% of its genomes")):
        for cname, ex in cats_opts:
            add("L. Sandwiched", f"closest gene on both sides has the same category, {text}, {cname}",
                sandwiched, exclude=ex, min_genomes=genomes, min_frac=frac)

    # M. two methods must agree (4)
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

    # N. enrichment over background (4)
    for n in (3, 5):
        for cname, ex in cats_opts:
            add("N. Enrichment over background", f"category most enriched within {n} neighbours compared "
                f"with the whole data, found in >= 2 genomes, {cname}", enriched, n=n, exclude=ex)

    # O. weighted vote scaled by conservation (4)
    for curve, n in (("very_steep", 3), ("linear", 5)):
        for cname, ex in cats_opts:
            add("O. Conservation-weighted vote", f"{curve_name[curve]} weights, {n} neighbours, each "
                f"category's vote times the share of genomes that have it, {cname}",
                conserved_vote, curve=curve, n=n, exclude=ex)

    # P. confidence cutoffs on very steep weights, 3 neighbours (8)
    for cut in (50, 60, 70, 80):
        for cname, ex in cats_opts:
            add("P. Confidence cutoff (very steep)", f"very steep weights, 3 neighbours, predict only if "
                f"winner has >= {cut}% of the vote, {cname}", weighted_vote,
                curve="very_steep", n=3, exclude=ex, cutoff=cut)

    # ---- added rules C151-C184: gene strands and distances (need smorf_neighbourhoods.sqlite) ----

    # Q. operon vote (12)
    for gap in (50, 100, 200):
        for curve in ("flat", "linear"):
            for cname, ex in cats_opts:
                add("Q. Operon vote", f"{curve} weights, only neighbours in the target gene's operon (same "
                    f"strand, <= {gap} bp between genes), {cname}", operon_vote,
                    max_gap=gap, curve=curve, n=5, exclude=ex)

    # R. sandwiched, with strand and distance (12)
    for strand, dist, text in ((True, None, "both on the target's strand"),
                               (False, 100, "both within 100 bp of the target gene"),
                               (True, 100, "both on the target's strand and within 100 bp")):
        for genomes, frac, need in ((2, None, "in at least 2 genomes"), (1, 0.5, "in at least 50% of its genomes")):
            for cname, ex in cats_opts:
                add("R. Sandwiched (strand/distance)", f"closest gene on both sides has the same category, "
                    f"{text}, {need}, {cname}", sandwiched, same_strand=strand, max_dist=dist,
                    exclude=ex, min_genomes=genomes, min_frac=frac)

    # S. vote weighted by distance in bp (6)
    for half in (250, 500, 1000):
        for cname, ex in cats_opts:
            add("S. Distance-weighted vote", f"5 neighbours, vote halves every {half} bp from the target "
                f"gene, {cname}", distance_vote, half=half, n=5, exclude=ex)

    # T. operon vote with confidence cutoff (4)
    for cut in (50, 70):
        for cname, ex in cats_opts:
            add("T. Operon vote + cutoff", f"linear weights, operon neighbours (<= 100 bp), predict only if "
                f"winner has >= {cut}% of the vote, {cname}", operon_vote,
                max_gap=100, curve="linear", n=5, exclude=ex, cutoff=cut)

    # ---- added rules C185-C188: nearest-neighbour vote with Annoy ----

    # U. Annoy nearest neighbours (4)
    for k in (10, 25):
        for cname, ex in cats_opts:
            add("U. Annoy nearest neighbours", f"the {k} smORFs with the most similar neighbourhoods "
                f"(COG IDs within 3 genes, Annoy) vote, {cname}", similar_vote, engine="annoy", k=k, exclude=ex)

    # ---- added rules C189-C192: the same with Voyager ----

    # V. Voyager nearest neighbours (4)
    for k in (10, 25):
        for cname, ex in cats_opts:
            add("V. Voyager nearest neighbours", f"the {k} smORFs with the most similar neighbourhoods "
                f"(COG IDs within 3 genes, Voyager) vote, {cname}", similar_vote, engine="voyager", k=k,
                exclude=ex)
    return combos


# ===================================================== PART 1: PREDICT ==

def predict(d, combos, save=True):
    """Every combination predicts for every smORF, annotated or not.

    Writes predictions_all_smorfs.tsv (one row per smORF, one column per
    rule) if save, and returns the predictions and their pct, as dicts by
    rule ID."""
    preds, pcts = {}, {}
    for c in combos:
        t = time.time()
        preds[c["id"]], pcts[c["id"]] = p, s = c["fn"](d, **c["settings"])
        made = (p < d.n_cats).sum()
        print(f"  {c['id']}  predicted {made / d.smorfs.height * 100:5.1f}% of smORFs  "
              f"({time.time() - t:.1f}s)  {c['rule']}", flush=True)

    if not save:
        return preds, pcts
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


def evaluate(d, combos, preds, pcts, save=True):
    """Score each combination ONLY on smORFs with a known COG category.

    A smORF is evaluated if the rule predicts a category for it and it has
    a known category not excluded by the rule; the prediction is correct if
    it is any of its known categories. Coverage is the share of annotated
    smORFs that are evaluated. Accuracy is also given per confidence tier
    (TIERS, by pct). Shows every rule, most accurate first, and, if save,
    writes combinations_results.csv and results_by_bin.csv."""
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
        if not save:
            break
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
    db = args[0] if args else "smorfs_all_bins.sqlite"
    t0 = time.time()
    d = load_db(db)
    print(f"Loaded {d.smorfs.height:,} smORFs in {time.time() - t0:.1f}s", flush=True)
    combos = build_combinations(d)

    print(f"\nPART 1 - PREDICT: {len(combos)} combinations on all smORFs", flush=True)
    preds, pcts = predict(d, combos, save)

    print("\nPART 2 - EVALUATE: annotated smORFs only", flush=True)
    evaluate(d, combos, preds, pcts, save)

    if not save:
        print(f"\nDone in {(time.time() - t0) / 60:.1f} min. No files saved (add --save to save them).")
        return
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Saved:")
    print("  combinations_results.csv   the sheet: one row per rule, accuracy in %")
    print("  results_by_bin.csv         the same, per bin")
    print("  predictions_all_smorfs.tsv every smORF x every rule (C001...)")


if __name__ == "__main__":
    main()
