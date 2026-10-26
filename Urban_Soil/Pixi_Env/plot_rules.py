#!/usr/bin/env python3
"""
Accuracy against coverage for every rule in combinations_results.csv (the
output of sweep_db.py): one dot per rule, coloured by rule group.

A rule is on the best-trade-off line (Pareto front) if no other rule is both
more accurate and covers more smORFs. These rules are labelled, and listed
in rules_pareto.tsv, from highest coverage to highest accuracy. The control
rules (group J) are drawn as crosses and never on the line: they show what
accuracy looks like without real neighbourhood signal.

    python3 plot_rules.py [folder with combinations_results.csv]
        -> rules_accuracy_coverage.png, rules_accuracy_coverage.pdf, rules_pareto.tsv
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ACC = "Accuracy (%)"
COV = "Coverage (% of annotated smORFs evaluated)"


def pareto(rules):
    """The rules that no other rule beats on both accuracy and coverage,
    from highest coverage down"""
    front, best_acc = [], -1.0
    for r in sorted(rules, key=lambda r: (-r["cov"], -r["acc"])):
        if r["acc"] > best_acc:
            front.append(r)
            best_acc = r["acc"]
    return front


def main():
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    with open(folder / "combinations_results.csv", newline="") as fh:
        rules = [r for r in csv.DictReader(fh) if r[ACC] and r[COV]]
    for r in rules:
        r["acc"], r["cov"] = float(r[ACC]), float(r[COV])
    controls = [r for r in rules if r["Group"].startswith("J.")]
    real = [r for r in rules if not r["Group"].startswith("J.")]
    front = pareto(real)

    groups = sorted({r["Group"] for r in real})
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(11, 7))
    for i, g in enumerate(groups):
        rs = [r for r in real if r["Group"] == g]
        ax.scatter([r["cov"] for r in rs], [r["acc"] for r in rs], s=28, alpha=0.8,
                   color=cmap(i % 20), label=f"{g} ({len(rs)})", edgecolor="none")
    ax.scatter([r["cov"] for r in controls], [r["acc"] for r in controls], marker="x", s=50,
               color="black", label=f"J. Control ({len(controls)})")
    ax.plot([r["cov"] for r in front], [r["acc"] for r in front], color="grey", lw=1, ls="--",
            zorder=0, label="best trade-off (Pareto front)")
    for r in front:
        ax.annotate(r["ID"], (r["cov"], r["acc"]), textcoords="offset points", xytext=(4, 3), fontsize=7)

    ax.set_xlabel("Coverage (% of annotated smORFs that get a prediction)")
    ax.set_ylabel("Accuracy (% of predictions that are a known COG category)")
    ax.set_title(f"{len(rules)} rules: accuracy against coverage")
    ax.set_xlim(0, 102)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(folder / f"rules_accuracy_coverage.{ext}", dpi=200)

    with open(folder / "rules_pareto.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ID", "Group", "Rule", ACC, COV, "Evaluated (n)"])
        for r in front:
            w.writerow([r["ID"], r["Group"], r["Rule"], r[ACC], r[COV], r["Evaluated (n)"]])
    print(f"{len(rules)} rules, {len(front)} on the Pareto front:")
    for r in front:
        print(f"  {r['ID']}\t{r['acc']:5.1f}% accuracy\t{r['cov']:5.1f}% coverage\t{r['Rule']}")
    print(f"Saved rules_accuracy_coverage.png/.pdf and rules_pareto.tsv in {folder}")


if __name__ == "__main__":
    main()
