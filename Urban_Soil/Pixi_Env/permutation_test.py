#!/usr/bin/env python3
"""
Is the neighbourhood signal real? A permutation test (the statistical
validation of Aim 1 in the confirmation document).

Three simple rules are scored on the real neighbourhoods, then again on
neighbourhoods randomised in two ways, many times, to get what accuracy
looks like by chance:

  order   Inside each occurrence's window, the neighbour genes are put in a
          random order (same genes, same categories, different positions).
          Tests: does it matter which genes are CLOSEST to the smORF?
  swap    Every occurrence gets the window of another occurrence, picked at
          random from all smORFs. Tests: does a smORF's neighbourhood say
          anything about THIS smORF at all?

The rules (all categories, votes added up over all occurrences of a smORF,
a tie is no prediction):

  nearest   the closest gene on each side votes for its categories
  weighted  the 5 genes on each side vote, the gene at rank r with weight 6 - r
  sandwich  only occurrences where the closest gene on both sides shares a
            category vote, for that shared category

For each rule and randomisation: the real accuracy, the mean and spread of
the random accuracies, z = (real - mean) / sd, and an empirical p-value,
(1 + randomisations at least as accurate as the real data) / (1 + N).
Coverage is given too: the % of annotated smORFs evaluated.

    python3 permutation_test.py DB [N=100] [SEED=0]
        -> permutation_test.tsv (summary) and permutation_null.tsv (every
           randomisation), in the current folder
"""

import csv
import json
import sqlite3
import sys
import time

import numpy as np

WINDOW = 5
RULES = ("nearest", "weighted", "sandwich")


class Neighbourhoods:
    """Every neighbour gene of every occurrence, as numpy arrays.

    Genes: occ (occurrence), up (1 upstream, 0 downstream), rank (1 is the
    closest), dist (bp from the target gene). Category rows: gene, cat (a
    gene with several categories has several rows, a gene with none has
    none). Each occurrence's genes are its "slots": randomising changes
    which gene sits in which slot."""

    def __init__(self, db_path):
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        rows = conn.execute("SELECT bin_name, smorf_id, occurrences_json FROM smorfs "
                            "ORDER BY bin_name, smorf_id")
        cats, known = {}, []
        occ_smorf, g_occ, g_up, g_rank, g_dist, c_gene, c_cat = [], [], [], [], [], [], []
        for s, (_, _, text) in enumerate(rows):
            k = set()
            for o in json.loads(text):
                occ = len(occ_smorf)
                occ_smorf.append(s)
                k.update(o["target_categories"])
                for c in o["target_categories"]:
                    cats.setdefault(c, len(cats))
                for (rank, gene_cats), info in zip(o["neighbours"], o["neighbour_info"]):
                    for c in gene_cats:
                        c_gene.append(len(g_occ))
                        c_cat.append(cats.setdefault(c, len(cats)))
                    g_occ.append(occ)
                    g_up.append(info[0] == "up")
                    g_rank.append(rank)
                    g_dist.append(info[2])
            known.append(k)
        conn.close()
        self.cats = sorted(cats, key=cats.get)
        self.n_cats = len(self.cats)
        self.n_smorfs = len(known)
        self.occ_smorf = np.array(occ_smorf, dtype=np.int64)
        self.occ = np.array(g_occ, dtype=np.int64)
        self.up = np.array(g_up, dtype=bool)
        self.rank = np.array(g_rank, dtype=np.int64)
        self.dist = np.array(g_dist, dtype=np.int64)
        self.c_gene = np.array(c_gene, dtype=np.int64)
        self.c_cat = np.array(c_cat, dtype=np.int64)
        # known[s, c]: c is a known category of smORF s
        self.known = np.zeros((self.n_smorfs, self.n_cats), dtype=bool)
        for s, k in enumerate(known):
            self.known[s, [cats[c] for c in k]] = True
        self.annotated = self.known.any(axis=1)

    def real(self):
        """Gene i sits in its own slot"""
        return np.arange(len(self.occ))

    def shuffle_order(self, rng):
        """gene_in_slot: inside each occurrence, the genes in a random order"""
        order = np.lexsort((rng.random(len(self.occ)), self.occ))  # genes by occ, random within
        slots = np.argsort(self.occ, kind="stable")  # slots by occ, in their order
        gene_in_slot = np.empty(len(self.occ), dtype=np.int64)
        gene_in_slot[slots] = order
        return gene_in_slot

    def swap(self, rng):
        """Each occurrence's slots take the genes of another, random,
        occurrence (as far as it has genes for them: a slot whose side and
        rank the other window lacks stays empty)"""
        n_occ = len(self.occ_smorf)
        other = rng.permutation(n_occ)
        key = lambda occ, up, rank: (occ * 2 + up) * (WINDOW + 1) + np.minimum(rank, WINDOW)
        index = np.full(n_occ * 2 * (WINDOW + 1), -1, dtype=np.int64)
        index[key(self.occ, self.up, self.rank)[::-1]] = np.arange(len(self.occ))[::-1]  # first wins
        return index[key(other[self.occ], self.up, self.rank)]


def predict(nb, smorf, cat, weight):
    """The category with the most votes for each smORF (-1: no votes or a tie)"""
    votes = np.bincount(smorf * nb.n_cats + cat, weights=weight,
                        minlength=nb.n_smorfs * nb.n_cats).reshape(nb.n_smorfs, nb.n_cats)
    top = votes.max(axis=1)
    best = votes.argmax(axis=1)
    ties = (np.abs(votes - top[:, None]) < 1e-9).sum(axis=1) > 1
    return np.where((top > 0) & ~ties, best, -1)


def score(nb, pred):
    """(accuracy %, coverage %) over the annotated smORFs"""
    ev = (pred >= 0) & nb.annotated
    ok = ev & nb.known[np.arange(nb.n_smorfs), np.maximum(pred, 0)]
    return 100 * ok.sum() / max(ev.sum(), 1), 100 * ev.sum() / max(nb.annotated.sum(), 1)


def run_rules(nb, gene_in_slot):
    """{rule: (accuracy, coverage)} with gene gene_in_slot[i] in slot i
    (-1: empty slot)"""
    slot_of = np.full(len(nb.occ), -1, dtype=np.int64)  # the slot each gene now sits in
    filled = gene_in_slot >= 0
    slot_of[gene_in_slot[filled]] = np.flatnonzero(filled)
    rows = slot_of[nb.c_gene] >= 0
    slot = slot_of[nb.c_gene[rows]]
    cat = nb.c_cat[rows]
    smorf = nb.occ_smorf[nb.occ[slot]]
    rank = nb.rank[slot]
    out = {}

    near = rank == 1
    out["nearest"] = score(nb, predict(nb, smorf[near], cat[near], np.ones(near.sum())))
    w = np.clip(WINDOW + 1 - rank, 0, None).astype(float)
    out["weighted"] = score(nb, predict(nb, smorf, cat, w))

    # sandwich: a category on the rank-1 gene of both sides of the same occurrence
    occ, up = nb.occ[slot[near]], nb.up[slot[near]]
    key = np.unique((occ * nb.n_cats + cat[near]) * 2 + up)  # (occ, cat, side), once each
    pair, n_sides = np.unique(key // 2, return_counts=True)
    both = pair[n_sides == 2]
    out["sandwich"] = score(nb, predict(nb, nb.occ_smorf[both // nb.n_cats], both % nb.n_cats,
                                        np.ones(len(both))))
    return out


def main(db_path, n=100, seed=0):
    t0 = time.time()
    nb = Neighbourhoods(db_path)
    print(f"{nb.n_smorfs:,} smORFs ({nb.annotated.sum():,} annotated), {len(nb.occ_smorf):,} occurrences, "
          f"{len(nb.occ):,} neighbour genes ({time.time() - t0:.0f}s)", flush=True)
    real = run_rules(nb, nb.real())
    rng = np.random.default_rng(seed)
    null = {"order": [], "swap": []}
    with open("permutation_null.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["randomisation", "i", "rule", "accuracy (%)", "coverage (%)"])
        for i in range(n):
            for how, f in (("order", nb.shuffle_order), ("swap", nb.swap)):
                res = run_rules(nb, f(rng))
                null[how].append(res)
                for r in RULES:
                    w.writerow([how, i, r, f"{res[r][0]:.2f}", f"{res[r][1]:.2f}"])
            if (i + 1) % 10 == 0:
                print(f"  {i + 1}/{n} randomisations ({time.time() - t0:.0f}s)", flush=True)

    head = ["rule", "randomisation", "real accuracy (%)", "random accuracy mean (%)", "random accuracy sd",
            "z", "p", "real coverage (%)", "random coverage mean (%)", "N"]
    lines = ["\t".join(head)]
    for r in RULES:
        for how in null:
            acc = np.array([x[r][0] for x in null[how]])
            cov = np.array([x[r][1] for x in null[how]])
            sd = acc.std(ddof=1) if n > 1 else 0.0
            z = (real[r][0] - acc.mean()) / sd if sd > 0 else float("inf")
            p = (1 + (acc >= real[r][0]).sum()) / (1 + n)
            lines.append("\t".join([r, how, f"{real[r][0]:.1f}", f"{acc.mean():.1f}", f"{sd:.2f}", f"{z:.1f}",
                                    f"{p:.3g}", f"{real[r][1]:.1f}", f"{cov.mean():.1f}", str(n)]))
    with open("permutation_test.tsv", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nSaved permutation_test.tsv and permutation_null.tsv ({(time.time() - t0) / 60:.1f} min)")


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], *(int(x) for x in sys.argv[2:]))
