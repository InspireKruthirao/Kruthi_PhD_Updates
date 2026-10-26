#!/usr/bin/env python3
"""
Does the neighbourhood signal fade with distance from the smORF?

For the annotated smORFs, every neighbour gene with a COG category is a
"hit" if it shares a category with the smORF. The hit rate (% of neighbour
genes that are hits) is reported:

  - by rank (1 = closest gene, 5 = farthest), for genes on the same strand
    as the target gene and on the opposite strand
  - by distance from the target gene, in bp

The background is the hit rate when every occurrence is given the
neighbourhood of a random other occurrence (as "swap" in
permutation_test.py): what the hit rate is with no real signal.

Also, how good a prediction is from the genes at one rank only: the genes
at rank r on both sides vote (all categories, over all occurrences, a tie
is no prediction), r = 1...5, and for the context size: the genes at ranks
1...r vote (weight 1 each).

    python3 distance_decay.py DB [SEED=0]
        -> distance_decay.tsv and distance_decay.png/.pdf in the current folder
"""

import sys
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from permutation_test import WINDOW, Neighbourhoods, predict, score

DIST_BINS = (0, 100, 250, 500, 1000, 2000, 5000, 10**9)  # bp


def hits(nb, gene_in_slot):
    """For every filled slot whose gene has a category and whose smORF is
    annotated: the slot and whether the gene shares a category with the
    smORF of the slot"""
    slot = np.flatnonzero(gene_in_slot >= 0)
    gene = gene_in_slot[slot]
    slot_of_gene = np.full(len(nb.occ), -1, dtype=np.int64)
    slot_of_gene[gene] = slot
    rows = slot_of_gene[nb.c_gene] >= 0
    row_slot = slot_of_gene[nb.c_gene[rows]]
    smorf = nb.occ_smorf[nb.occ[row_slot]]
    row_hit = nb.known[smorf, nb.c_cat[rows]]
    has_cat = np.bincount(row_slot, minlength=len(nb.occ)) > 0
    hit = np.bincount(row_slot, weights=row_hit, minlength=len(nb.occ)) > 0
    keep = has_cat & nb.annotated[nb.occ_smorf[nb.occ]]
    return np.flatnonzero(keep), hit[keep]


def rate(hit, mask):
    """(hit rate %, n); the rate is NaN if n is 0"""
    n = int(mask.sum())
    return (100 * hit[mask].sum() / n if n else float("nan")), n


def f(x):
    return "" if x != x else f"{x:.1f}"


def main(db_path, seed=0):
    t0 = time.time()
    nb = Neighbourhoods(db_path)
    print(f"{nb.n_smorfs:,} smORFs ({nb.annotated.sum():,} annotated), {len(nb.occ):,} neighbour genes "
          f"({time.time() - t0:.0f}s)", flush=True)
    slot, hit = hits(nb, nb.real())
    b_slot, b_hit = hits(nb, nb.swap(np.random.default_rng(seed)))
    lines = []
    def say(text=""):
        print(text, flush=True)
        lines.append(text)

    say("1. Hit rate by rank (% of neighbour genes sharing a category with the smORF)")
    say("rank\tsame strand (%)\tn\topposite strand (%)\tn\tall (%)\tn\tbackground (%)\tn")
    by_rank = []
    for r in range(1, WINDOW + 1):
        at, b_at = nb.rank[slot] == r, nb.rank[b_slot] == r
        same, opp, a, bg = (rate(hit, at & nb.same[slot]), rate(hit, at & ~nb.same[slot]),
                            rate(hit, at), rate(b_hit, b_at))
        by_rank.append((r, same[0], opp[0], a[0], bg[0]))
        say(f"{r}\t{f(same[0])}\t{same[1]:,}\t{f(opp[0])}\t{opp[1]:,}\t{f(a[0])}\t{a[1]:,}\t{f(bg[0])}\t{bg[1]:,}")

    say("\n2. Hit rate by distance from the target gene")
    say("distance (bp)\tsame strand (%)\tn\topposite strand (%)\tn\tall (%)\tn\tbackground (%)\tn")
    by_dist = []
    for lo, hi in zip(DIST_BINS, DIST_BINS[1:]):
        at, b_at = ((nb.dist[s] >= lo) & (nb.dist[s] < hi) for s in (slot, b_slot))
        same, opp, a, bg = (rate(hit, at & nb.same[slot]), rate(hit, at & ~nb.same[slot]),
                            rate(hit, at), rate(b_hit, b_at))
        label = f"{lo}-{hi}" if hi < 10**9 else f">={lo}"
        by_dist.append((label, same[0], opp[0], a[0], bg[0]))
        say(f"{label}\t{f(same[0])}\t{same[1]:,}\t{f(opp[0])}\t{opp[1]:,}\t{f(a[0])}\t{a[1]:,}\t{f(bg[0])}\t{bg[1]:,}")

    # votes of the genes at one rank, and at ranks 1...r
    cat_slot = nb.c_gene  # real data: gene i in slot i
    smorf, rank = nb.occ_smorf[nb.occ[cat_slot]], nb.rank[cat_slot]
    say("\n3. Prediction from the genes at rank r only, and from ranks 1...r (context size)")
    say("r\trank r only: accuracy (%)\tcoverage (%)\tranks 1...r: accuracy (%)\tcoverage (%)")
    by_vote = []
    for r in range(1, WINDOW + 1):
        one, upto = rank == r, rank <= r
        acc1, cov1 = score(nb, predict(nb, smorf[one], nb.c_cat[one], np.ones(one.sum())))
        acc2, cov2 = score(nb, predict(nb, smorf[upto], nb.c_cat[upto], np.ones(upto.sum())))
        by_vote.append((r, acc1, cov1, acc2, cov2))
        say(f"{r}\t{acc1:.1f}\t{cov1:.1f}\t{acc2:.1f}\t{cov2:.1f}")

    with open("distance_decay.tsv", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    ranks = [x[0] for x in by_rank]
    for i, (name, style) in enumerate((("same strand", "o-"), ("opposite strand", "s-"), ("all", "^-"),
                                       ("background", "k--")), 1):
        ax.plot(ranks, [x[i] for x in by_rank], style, label=name)
    ax.set_xlabel("Rank of the neighbour gene (1 = closest)")
    ax.set_ylabel("Hit rate (% sharing a category with the smORF)")
    ax.set_xticks(ranks)
    ax.set_title("By rank")
    ax.legend(fontsize=8)
    ax = axes[1]
    xs = range(len(by_dist))
    for i, (name, style) in enumerate((("same strand", "o-"), ("opposite strand", "s-"), ("all", "^-"),
                                       ("background", "k--")), 1):
        ax.plot(xs, [x[i] for x in by_dist], style, label=name)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([x[0] for x in by_dist], rotation=30, fontsize=8)
    ax.set_xlabel("Distance from the target gene (bp)")
    ax.set_title("By distance")
    ax = axes[2]
    ax.plot(ranks, [x[1] for x in by_vote], "o-", label="rank r only: accuracy")
    ax.plot(ranks, [x[3] for x in by_vote], "s-", label="ranks 1...r: accuracy")
    ax.plot(ranks, [x[4] for x in by_vote], "s:", color="grey", label="ranks 1...r: coverage")
    ax.set_xticks(ranks)
    ax.set_xlabel("r")
    ax.set_ylabel("%")
    ax.set_title("Vote of the genes at rank r / ranks 1...r")
    ax.legend(fontsize=8)
    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"distance_decay.{ext}", dpi=200)
    print(f"\nSaved distance_decay.tsv and distance_decay.png/.pdf ({(time.time() - t0) / 60:.1f} min)")


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], *(int(x) for x in sys.argv[2:]))
