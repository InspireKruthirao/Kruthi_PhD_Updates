#!/usr/bin/env python3
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import numpy as np

BASE = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_26_30_no_min_overlap"]

MIN_SUPPORT = 3

RULE_KEYS = [
    "R1_all3",
    "R2_both_min3",
    "R3_down",
    "R4_up",
    "R5_either"
]

RULE_LABELS = {
    "R1_all3":      "All three match (upstream + downstream + smORF)",
    "R2_both_min3": f"Both upstream & downstream same (min {MIN_SUPPORT})",
    "R3_down":      "Downstream neighbour matches",
    "R4_up":        "Upstream neighbour matches",
    "R5_either":    "Either neighbour matches",
}

# ─────────────────────────────────────────────────────────────────────────────
# Your existing helper functions (parse_attrs, cog_set, is_annotated, majority, read_gff,
# aggregate_smorf_evidence, match_category) go here unchanged
# ─────────────────────────────────────────────────────────────────────────────

def evaluate_rules(ev):
    tc = ev["target"]
    uc = ev["upstream"]
    dc = ev["downstream"]
    is_ann = ev["n_annotated"] > ev["n_unannotated"]

    r = {}

    el1 = is_annotated(uc) and is_annotated(dc) and uc == dc and is_annotated(tc) and tc == uc
    r["R1_all3"] = {"eligible": el1, "correct": el1}

    el2 = is_annotated(uc) and is_annotated(dc) and uc == dc
    ok2 = el2 and is_ann and tc == uc
    r["R2_both_min3"] = {"eligible": el2, "correct": ok2}

    el3 = is_annotated(dc)
    ok3 = el3 and is_ann and tc == dc
    r["R3_down"] = {"eligible": el3, "correct": ok3}

    el4 = is_annotated(uc)
    ok4 = el4 and is_ann and tc == uc
    r["R4_up"] = {"eligible": el4, "correct": ok4}

    el5 = is_annotated(uc) or is_annotated(dc)
    ok5 = el5 and is_ann and (tc == uc or tc == dc)
    r["R5_either"] = {"eligible": el5, "correct": ok5}

    return r, is_ann

def plot_line_cov_acc(counters, total_ann, title):
    covs, accs, lbls = [], [], []
    for k in RULE_KEYS:
        el = counters[k]["eligible"]
        co = counters[k]["correct"]
        cov = el / total_ann * 100 if total_ann else 0
        acc = co / el * 100 if el else 0
        covs.append(cov)
        accs.append(acc)
        lbls.append(RULE_LABELS[k])

    fig, ax1 = plt.subplots(figsize=(10, 5.5))

    ax1.plot(range(len(RULE_KEYS)), covs, marker='o', color="#1f77b4", label="Coverage %", linewidth=2.2)
    ax1.set_ylabel("Coverage (%)", color="#1f77b4")
    ax1.tick_params(axis='y', labelcolor="#1f77b4")
    ax1.set_ylim(0, 110)

    ax2 = ax1.twinx()
    ax2.plot(range(len(RULE_KEYS)), accs, marker='s', color="#ff7f0e", label="Accuracy %", linewidth=2.2)
    ax2.set_ylabel("Accuracy (%)", color="#ff7f0e")
    ax2.tick_params(axis='y', labelcolor="#ff7f0e")
    ax2.set_ylim(0, 110)

    ax1.set_xticks(range(len(RULE_KEYS)))
    ax1.set_xticklabels(lbls, rotation=25, ha="right", fontsize=9)
    ax1.set_title(f"Coverage vs Accuracy Trade-off — {title}")
    ax1.grid(True, ls="--", alpha=0.35)

    lines1, lbls1 = ax1.get_legend_handles_labels()
    lines2, lbls2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, lbls1 + lbls2, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)

    plt.tight_layout()
    plt.show()

def main():
    for folder in RANGE_DIRS:
        root = BASE / folder
        print(f"{'═'*69}")
        print(f"Rule performance — {folder}")
        print(f"{'═'*69}")

        counters = {k: {"eligible": 0, "correct": 0} for k in RULE_KEYS}
        match_cnt = defaultdict(int)
        total_sm = ann_sm = 0

        unann_any = unann_up = unann_down = unann_both_same = unann_both_diff = 0

        for d in root.iterdir():
            if not d.is_dir() or not d.name.startswith("SHD1_SM.100AA"):
                continue
            total_sm += 1

            ev = aggregate_smorf_evidence(d)
            rules, is_ann = evaluate_rules(ev)

            if is_ann:
                ann_sm += 1
                for k in RULE_KEYS:
                    if rules[k]["eligible"]:
                        counters[k]["eligible"] += 1
                        counters[k]["correct"] += int(rules[k]["correct"])

                cat = match_category(ev, is_ann)
                if cat:
                    match_cnt[cat] += 1
            else:
                has_up = bool(ev["upstream"])
                has_down = bool(ev["downstream"])
                if has_up or has_down:
                    unann_any += 1
                    if has_up: unann_up += 1
                    if has_down: unann_down += 1
                    if has_up and has_down:
                        if ev["upstream"] == ev["downstream"]:
                            unann_both_same += 1
                        else:
                            unann_both_diff += 1

        # Print rule performance table
        print(f"{'Rule':<45} {'Eligible':>9} {'Correct':>9} {'Cov%':>7} {'Acc%':>7}")
        print("─"*69)
        for k in RULE_KEYS:
            el = counters[k]["eligible"]
            co = counters[k]["correct"]
            cov = el / ann_sm * 100 if ann_sm else 0
            acc = co / el * 100 if el else 0
            flag = "✓" if acc >= ACCURACY_TARGET else ""
            print(f"{RULE_LABELS[k]:<45} {el:>9} {co:>9} {cov:>7.1f} {acc:>7.1f}{flag}")

        print(f"{'═'*69}")
        print(f"{folder}")
        print(f"{'═'*69}")

        print(f"Total unique smORFs:      {total_sm:,}")
        print(f"Annotated:                {ann_sm:,} ({ann_sm/total_sm*100:.1f}% if total_sm else 0)")
        print(f"Unannotated:              {total_sm - ann_sm:,}\n")

        print("Unannotated smORFs with annotated neighbours:")
        print(f"  • Has any annotated neighbour:   {unann_any:,}")
        print(f"  • Upstream neighbour annotated:  {unann_up:,}")
        print(f"  • Downstream neighbour annotated: {unann_down:,}")
        print(f"  • Both neighbours same function: {unann_both_same:,}")
        print(f"  • Both annotated but different:  {unann_both_diff:,}\n")

        # Show the line plot
        plot_line_cov_acc(counters, ann_sm, folder)

if __name__ == "__main__":
    main()
