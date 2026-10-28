#!/usr/bin/env python3
"""
Plots of the rule results for a supervisor meeting.

    python3 meeting_plots.py [sweep_results folder] [output folder]

Reads combinations_results.csv, predictions_all_smorfs.tsv and
cascade_predictions.tsv (from run_sweep.sh) and saves, as PNG and PDF:

  1_accuracy_vs_coverage   every rule as a dot; best trade-offs joined; key rules named
  2_key_rules              accuracy and coverage of the key rules, side by side
  3_new_predictions        per key rule: smORFs predicted that are already annotated
                           vs smORFs with no known function (new predictions)
  4_cascade                the cascade step by step: smORFs predicted (annotated / new)
                           and the expected accuracy of each step

Also prints the numbers behind each plot.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl

RESULTS = "/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL/sweep_results"
NAVY, ORANGE, GREY = "#1f3a5f", "#c0581f", "#b8bec6"
KEY = {  # rule -> short name on the plots
    "C172": "Sandwich, same strand, <=100 bp",
    "C166": "Sandwich, same strand, >=50% genomes",
    "C126": "Sandwich, >=1 genome",
    "C184": "Operon vote, >=70%",
    "C189": "Similar neighbourhoods (Voyager)",
    "C185": "Similar neighbourhoods (Annoy)",
    "C175": "Distance-weighted vote",
    "C027": "Original method",
    "C099": "Control: random",
}
OFFSET = {  # label position of each key rule (points), so close rules do not overlap
    "C172": (8, 6), "C166": (8, -12), "C126": (8, 2), "C184": (8, 2), "C189": (-38, 6),
    "C185": (8, -12), "C175": (-38, 6), "C027": (-38, -12), "C099": (-38, 4)}
MIN_EVALUATED = 500  # a rule must be checked on at least this many smORFs to be on the best trade-off line
NO_CALL = {"NO_EVIDENCE", "UNRESOLVED", "LOW_CONFIDENCE", "TOO_FEW_OCCURRENCES", "NO_PREDICTION", ""}

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlesize": 13, "savefig.dpi": 200,
                     "savefig.bbox": "tight"})


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}")
    plt.close(fig)
    print(f"  saved {out / name}.png/.pdf\n")


def pareto(df):
    """Rules no other rule beats on both accuracy and coverage"""
    front, best = [], -1.0
    for r in df.sort(["cov", "acc"], descending=True).iter_rows(named=True):
        if r["acc"] > best:
            front.append(r)
            best = r["acc"]
    return pl.DataFrame(front)


def plot_accuracy_coverage(rules, out):
    real = rules.filter(~pl.col("Group").str.starts_with("J."))
    ctrl = rules.filter(pl.col("Group").str.starts_with("J."))
    front = pareto(real.filter(pl.col("Evaluated (n)").cast(pl.Int64) >= MIN_EVALUATED))
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(real["cov"], real["acc"], s=22, color=GREY, label=f"{real.height} rules")
    ax.scatter(ctrl["cov"], ctrl["acc"], s=40, marker="x", color=ORANGE, label="Controls (random)")
    ax.plot(front["cov"], front["acc"], color=NAVY, lw=2, label="Best trade-off")
    for r in rules.filter(pl.col("ID").is_in(list(KEY))).iter_rows(named=True):
        ax.scatter(r["cov"], r["acc"], s=60, color=NAVY, zorder=3)
        ax.annotate(r["ID"], (r["cov"], r["acc"]), xytext=OFFSET[r["ID"]], textcoords="offset points",
                    fontsize=9, fontweight="bold", color=NAVY)
    ax.set(xlim=(0, 105), ylim=(0, 100), xlabel="Coverage (% of annotated smORFs predicted)",
           ylabel="Accuracy (% of predictions correct)",
           title="Strict rules are accurate, loose rules cover more")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, loc="upper right")
    print("1. Accuracy vs coverage: best trade-offs")
    for r in front.iter_rows(named=True):
        print(f"   {r['ID']}  accuracy {r['acc']:5.1f}%  coverage {r['cov']:5.1f}%  {r['Rule']}")
    save(fig, out, "1_accuracy_vs_coverage")


def plot_key_rules(rules, out):
    key = rules.filter(pl.col("ID").is_in(list(KEY))).sort("acc")
    names = [f"{i}  {KEY[i]}" for i in key["ID"]]
    y = range(key.height)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.barh([i + 0.2 for i in y], key["acc"], height=0.4, color=NAVY, label="Accuracy")
    ax.barh([i - 0.2 for i in y], key["cov"], height=0.4, color=ORANGE, label="Coverage")
    for i, (a, c) in enumerate(zip(key["acc"], key["cov"])):
        ax.text(a + 1, i + 0.2, f"{a:.1f}%", va="center", fontsize=8.5, color=NAVY)
        ax.text(c + 1, i - 0.2, f"{c:.1f}%", va="center", fontsize=8.5, color=ORANGE)
    ax.set_yticks(list(y), names)
    ax.set(xlim=(0, 112), xlabel="%")
    ax.set_title("Key rules: accuracy and coverage", pad=28)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    print("2. Key rules")
    for r in key.sort("acc", descending=True).iter_rows(named=True):
        print(f"   {r['ID']}  accuracy {r['acc']:5.1f}%  coverage {r['cov']:5.1f}%  {KEY[r['ID']]}")
    save(fig, out, "2_key_rules")


def plot_new_predictions(preds, out):
    rows = []
    for rid in KEY:
        if rid not in preds.columns or rid == "C099":  # (the random control predicts nothing real)
            continue
        made = preds.filter(~pl.col(rid).is_in(list(NO_CALL)))
        known = made.filter(pl.col("annotated") == "yes").height
        rows.append((rid, known, made.height - known))
    rows.sort(key=lambda r: r[1] + r[2])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    names = [f"{r[0]}  {KEY[r[0]]}" for r in rows]
    known = [r[1] for r in rows]
    new = [r[2] for r in rows]
    ax.barh(names, known, color=NAVY, label="Already annotated (checkable)")
    ax.barh(names, new, left=known, color=ORANGE, label="No known function (new prediction)")
    for i, (k, n) in enumerate(zip(known, new)):
        ax.text(k + n, i, f"  {n:,} new", va="center", fontsize=8.5, color=ORANGE)
    ax.set(xlabel="smORFs predicted", title=f"New functional predictions per rule (of {preds.height:,} smORFs)")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlim(0, max(k + n for k, n in zip(known, new)) * 1.2)
    ax.legend(frameon=False, loc="lower right")
    print(f"3. smORFs predicted per rule (all {preds.height:,} smORFs)")
    for r in sorted(rows, key=lambda r: -(r[1] + r[2])):
        print(f"   {r[0]}  predicted {r[1] + r[2]:>7,}  already annotated {r[1]:>7,}  new {r[2]:>7,}  {KEY[r[0]]}")
    save(fig, out, "3_new_predictions")


def plot_cascade(cas, out):
    made = cas.filter(pl.col("step") != "")
    steps = (made.group_by("step", "rule")
             .agg(known=(pl.col("annotated") == "yes").sum(), new=(pl.col("annotated") == "no").sum(),
                  acc=pl.col("expected accuracy (%)").first().cast(pl.Float64))
             .with_columns(pl.col("step").cast(pl.Int32)).sort("step"))
    labels = [f"{s}\n{r}" for s, r in zip(steps["step"], steps["rule"])]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(labels, steps["known"], color=NAVY, label="Already annotated")
    ax.bar(labels, steps["new"], bottom=steps["known"], color=ORANGE, label="No known function (new)")
    ax.set(xlabel="Cascade step and rule", ylabel="smORFs predicted",
           title=f"Cascade: {made.height:,} of {cas.height:,} smORFs predicted "
                 f"({(made['annotated'] == 'no').sum():,} new)")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax2 = ax.twinx()
    ax2.plot(labels, steps["acc"], color=NAVY, marker="o", lw=2, ls="--", label="Expected accuracy")
    for x, a in zip(labels, steps["acc"]):
        ax2.annotate(f"{a:.0f}%", (x, a), xytext=(0, 7), textcoords="offset points", ha="center",
                     fontsize=8.5, color=NAVY)
    ax2.set(ylim=(0, 105), ylabel="Expected accuracy of the step (%)")
    ax2.spines["right"].set_visible(True)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, loc="upper right")
    print("4. Cascade (all smORFs)")
    for r in steps.iter_rows(named=True):
        print(f"   step {r['step']}  {r['rule']}  expected accuracy {r['acc']:5.1f}%  "
              f"annotated {r['known']:>7,}  new {r['new']:>7,}")
    print(f"   not predicted: {cas.height - made.height:,}")
    save(fig, out, "4_cascade")


def main():
    folder = Path(sys.argv[1] if len(sys.argv) > 1 else RESULTS)
    out = Path(sys.argv[2] if len(sys.argv) > 2 else folder / "meeting_plots")
    out.mkdir(parents=True, exist_ok=True)
    rules = (pl.read_csv(folder / "combinations_results.csv", infer_schema_length=0)
             .with_columns(acc=pl.col("Accuracy (%)").cast(pl.Float64, strict=False),
                           cov=pl.col("Coverage (% of annotated smORFs evaluated)").cast(pl.Float64, strict=False))
             .drop_nulls(["acc", "cov"]))
    plot_accuracy_coverage(rules, out)
    plot_key_rules(rules, out)

    preds = pl.read_csv(folder / "predictions_all_smorfs.tsv", separator="\t", infer_schema_length=0)
    preds = preds.with_columns(pl.all().str.strip_chars("\r"))
    plot_new_predictions(preds, out)

    cascade = folder / "cascade_predictions.tsv"
    if cascade.exists():
        cas = (pl.read_csv(cascade, separator="\t", infer_schema_length=0)
               .with_columns(pl.all().str.strip_chars("\r").fill_null("")))
        plot_cascade(cas, out)
    print(f"Plots in {out}")


if __name__ == "__main__":
    main()
