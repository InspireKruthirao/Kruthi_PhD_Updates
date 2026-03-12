#!/usr/bin/env python3
"""
smORF neighbourhood rule evaluator — v6
Structure:
  1. Rule table  (eligible / correct / coverage / accuracy)
  2. COG analysis per rule:
       - which COG categories appear most in CORRECT predictions
       - which COG categories have the highest accuracy per rule
       - upstream / downstream / target breakdown
  3. Unannotated smORF summary
"""

from pathlib import Path
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── paths ──────────────────────────────────────────────────────────────────────
BASE       = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_21_25_no_min_overlap"]

# ── parameters ────────────────────────────────────────────────────────────────
N_NEIGHBOURS    = 5
ACCURACY_TARGET = 80.0
MIN_COG_COUNT   = 3    # minimum eligible count to report accuracy for a COG letter
# ──────────────────────────────────────────────────────────────────────────────

FEATURE_TYPE = "CDS"
UNKNOWN      = "Unknown"

RULE_KEYS = ["R1_both", "R2_down", "R3_up", "R4_either"]

RULE_LABELS = {
    "R1_both":   "Both neighbours agree",
    "R2_down":   "Downstream neighbour matches",
    "R3_up":     "Upstream neighbour matches",
    "R4_either": "Either neighbour matches",
}

# COG functional category descriptions + supergroup colours
COG_META = {
    # Information storage & processing
    "J": ("Translation, ribosomal structure",        "#e07b54"),
    "A": ("RNA processing and modification",          "#e07b54"),
    "K": ("Transcription",                            "#e07b54"),
    "L": ("Replication, recombination, repair",       "#e07b54"),
    "B": ("Chromatin structure and dynamics",         "#e07b54"),
    # Cellular processes & signalling
    "D": ("Cell cycle control, division",             "#5b8db8"),
    "V": ("Defense mechanisms",                       "#5b8db8"),
    "T": ("Signal transduction",                      "#5b8db8"),
    "M": ("Cell wall/membrane biogenesis",            "#5b8db8"),
    "N": ("Cell motility",                            "#5b8db8"),
    "U": ("Intracellular trafficking, secretion",     "#5b8db8"),
    "O": ("Post-translational modification",          "#5b8db8"),
    "W": ("Extracellular structures",                 "#5b8db8"),
    "X": ("Mobilome: prophages, transposons",         "#5b8db8"),
    # Metabolism
    "C": ("Energy production and conversion",         "#6aaa64"),
    "G": ("Carbohydrate transport and metabolism",    "#6aaa64"),
    "E": ("Amino acid transport and metabolism",      "#6aaa64"),
    "F": ("Nucleotide transport and metabolism",      "#6aaa64"),
    "H": ("Coenzyme transport and metabolism",        "#6aaa64"),
    "I": ("Lipid transport and metabolism",           "#6aaa64"),
    "P": ("Inorganic ion transport and metabolism",   "#6aaa64"),
    "Q": ("Secondary metabolites biosynthesis",       "#6aaa64"),
    # Poorly characterised
    "R": ("General function prediction only",         "#aaaaaa"),
    "S": ("Function unknown",                         "#aaaaaa"),
    "Z": ("Cytoskeleton",                             "#aaaaaa"),
    "Y": ("Nuclear structure",                        "#aaaaaa"),
}


# ══════════════════════════════════════════════════════════════════════════════
#  GFF / COG helpers
# ══════════════════════════════════════════════════════════════════════════════

def parse_attrs(s: str) -> dict:
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def cog_set(name: str):
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return frozenset()
    i = str(name).rfind("-")
    if i == -1:
        return frozenset()
    return frozenset(ch for ch in str(name)[i + 1:] if ch.isalpha())


def is_annotated(cog) -> bool:
    return bool(cog)


def majority(counter: Counter):
    if not counter:
        return None
    return counter.most_common(1)[0][0]


def read_gff_ordered(gff_file: Path):
    features = []
    with gff_file.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue
            a = parse_attrs(parts[8])
            features.append({
                "start":     int(parts[3]),
                "end":       int(parts[4]),
                "strand":    parts[6],
                "is_target": a.get("target", "0") == "1",
                "cog":       cog_set(a.get("Name", UNKNOWN)),
            })

    features.sort(key=lambda x: x["start"])
    target_strands = [f["strand"] for f in features if f["is_target"]]
    flipped = bool(target_strands and target_strands[0] == "-")

    if flipped:
        features = list(reversed(features))
        flip_map = {"+": "-", "-": "+"}
        for f in features:
            f["norm_strand"] = flip_map.get(f["strand"], f["strand"])
    else:
        for f in features:
            f["norm_strand"] = f["strand"]

    return features, flipped


# ══════════════════════════════════════════════════════════════════════════════
#  Evidence aggregation
# ══════════════════════════════════════════════════════════════════════════════

def aggregate_smorf_evidence(smorf_dir: Path) -> dict:
    target_votes     = Counter()
    up_votes         = Counter()
    down_votes       = Counter()
    multi_up_votes   = defaultdict(Counter)
    multi_down_votes = defaultdict(Counter)
    n_annotated = n_unannotated = n_total = n_flipped = 0

    for contig in smorf_dir.iterdir():
        if not contig.is_dir():
            continue
        for gff in contig.glob("*.gff"):
            seq, flipped = read_gff_ordered(gff)
            if flipped:
                n_flipped += 1
            for i, f in enumerate(seq):
                if not f["is_target"]:
                    continue
                n_total += 1
                tc = f["cog"]
                left_cog  = seq[i - 1]["cog"] if i - 1 >= 0       else frozenset()
                right_cog = seq[i + 1]["cog"] if i + 1 < len(seq) else frozenset()

                if is_annotated(tc):
                    n_annotated += 1
                    target_votes[tc] += 1
                else:
                    n_unannotated += 1

                if is_annotated(left_cog):
                    up_votes[left_cog] += 1
                if is_annotated(right_cog):
                    down_votes[right_cog] += 1

                for dist in range(1, N_NEIGHBOURS + 1):
                    ui, di = i - dist, i + dist
                    if ui >= 0 and is_annotated(seq[ui]["cog"]):
                        multi_up_votes[dist][seq[ui]["cog"]] += 1
                    if di < len(seq) and is_annotated(seq[di]["cog"]):
                        multi_down_votes[dist][seq[di]["cog"]] += 1

    dom_up   = majority(up_votes)   or frozenset()
    dom_down = majority(down_votes) or frozenset()

    return {
        "target":              majority(target_votes) or frozenset(),
        "upstream":            dom_up,
        "downstream":          dom_down,
        "multi_up_votes":      multi_up_votes,
        "multi_down_votes":    multi_down_votes,
        "n_annotated":         n_annotated,
        "n_unannotated":       n_unannotated,
        "n_total":             n_total,
        "n_flipped":           n_flipped,
    }


def is_majority_annotated(ev) -> bool:
    return ev["n_annotated"] > ev["n_unannotated"]


# ══════════════════════════════════════════════════════════════════════════════
#  Rule evaluation
# ══════════════════════════════════════════════════════════════════════════════

def evaluate_rules(ev) -> dict:
    tc, uc, dc = ev["target"], ev["upstream"], ev["downstream"]
    is_ann = is_majority_annotated(ev)
    both   = is_annotated(uc) and is_annotated(dc) and uc == dc

    return {
        "R1_both": {
            "eligible": both,
            "correct":  both   and is_ann and tc == uc,
            "pred_cog": uc     if both   else frozenset(),
        },
        "R2_down": {
            "eligible": is_annotated(dc),
            "correct":  is_annotated(dc) and is_ann and tc == dc,
            "pred_cog": dc,
        },
        "R3_up": {
            "eligible": is_annotated(uc),
            "correct":  is_annotated(uc) and is_ann and tc == uc,
            "pred_cog": uc,
        },
        "R4_either": {
            "eligible": is_annotated(uc) or is_annotated(dc),
            "correct":  (is_annotated(uc) or is_annotated(dc)) and is_ann and (tc == uc or tc == dc),
            "pred_cog": uc if is_annotated(uc) else dc,
        },
    }


# ══════════════════════════════════════════════════════════════════════════════
#  COG analysis helpers
# ══════════════════════════════════════════════════════════════════════════════

def letters(cog_fs) -> list:
    return list(cog_fs) if cog_fs else []


def build_cog_rule_stats(ann_evidence: list) -> dict:
    """
    For each rule, for each COG letter that appeared as the predicted COG:
      - how many times was it eligible
      - how many times was the prediction correct
    Also separately track target / upstream / downstream position frequencies.
    Returns:
      {
        rule_key: {
          "cog_eligible":  Counter(letter -> count),
          "cog_correct":   Counter(letter -> count),
        },
        "positions": {
          "target":     Counter(letter -> count),
          "upstream":   Counter(letter -> count),
          "downstream": Counter(letter -> count),
        }
      }
    """
    stats = {k: {"cog_eligible": Counter(), "cog_correct": Counter()}
             for k in RULE_KEYS}
    positions = {"target": Counter(), "upstream": Counter(), "downstream": Counter()}

    for ev in ann_evidence:
        rules = evaluate_rules(ev)

        # position frequencies
        for letter in letters(ev["target"]):
            positions["target"][letter] += 1
        for letter in letters(ev["upstream"]):
            positions["upstream"][letter] += 1
        for letter in letters(ev["downstream"]):
            positions["downstream"][letter] += 1

        for k in RULE_KEYS:
            r = rules[k]
            if not r["eligible"]:
                continue
            for letter in letters(r["pred_cog"]):
                stats[k]["cog_eligible"][letter] += 1
                if r["correct"]:
                    stats[k]["cog_correct"][letter] += 1

    return stats, positions


def cog_accuracy_table(cog_eligible: Counter, cog_correct: Counter,
                        min_count: int = MIN_COG_COUNT) -> list:
    """
    Returns list of (letter, eligible, correct, accuracy%) sorted by accuracy desc.
    Only includes letters with eligible >= min_count.
    """
    rows = []
    for letter, el in cog_eligible.items():
        if el < min_count:
            continue
        co  = cog_correct.get(letter, 0)
        acc = co / el * 100
        rows.append((letter, el, co, acc))
    rows.sort(key=lambda x: -x[3])
    return rows


# ══════════════════════════════════════════════════════════════════════════════
#  Plotting
# ══════════════════════════════════════════════════════════════════════════════

def _cog_bar(ax, letters_list, counts, title, colour_by_group=True):
    colours = [COG_META.get(l, ("", "#cccccc"))[1]
               for l in letters_list] if colour_by_group else ["#4c72b0"] * len(letters_list)
    bars = ax.bar(range(len(letters_list)), counts, color=colours,
                  edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(letters_list)))
    ax.set_xticklabels(letters_list, fontsize=10, fontweight="bold")
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("Count")
    ax.grid(True, axis="y", ls="--", alpha=0.3)
    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
                str(cnt), ha="center", va="bottom", fontsize=7)
    return bars


def plot_rule_overview(counters, ann_sm, title):
    """Coverage vs Accuracy for all 4 rules."""
    labels = [RULE_LABELS[k] for k in RULE_KEYS]
    covs   = [counters[k]["eligible"] / ann_sm * 100 if ann_sm else 0 for k in RULE_KEYS]
    accs   = [counters[k]["correct"] / counters[k]["eligible"] * 100
              if counters[k]["eligible"] else 0 for k in RULE_KEYS]

    fig, ax1 = plt.subplots(figsize=(9, 5))
    x = list(range(len(RULE_KEYS)))

    ax1.bar(x, covs, color="#1f77b4", alpha=0.6, label="Coverage %")
    ax1.set_ylabel("Coverage (%)", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax1.set_ylim(0, 120)

    ax2 = ax1.twinx()
    ax2.plot(x, accs, marker="o", color="#ff7f0e", linewidth=2.2,
             label="Accuracy %", zorder=5)
    ax2.axhline(ACCURACY_TARGET, color="#ff7f0e", linestyle=":", alpha=0.6,
                label=f"{ACCURACY_TARGET}% target")
    ax2.set_ylabel("Accuracy (%)", color="#ff7f0e")
    ax2.tick_params(axis="y", labelcolor="#ff7f0e")
    ax2.set_ylim(0, 110)

    for xi, (c, a) in enumerate(zip(covs, accs)):
        ax1.text(xi, c + 1, f"{c:.1f}%", ha="center", fontsize=8, color="#1f77b4")
        ax2.text(xi, a + 2, f"{a:.1f}%", ha="center", fontsize=8, color="#ff7f0e")

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=9, rotation=10, ha="right")
    ax1.set_title(f"Rule Coverage & Accuracy — {title}")
    ax1.grid(True, axis="y", ls="--", alpha=0.3)

    lines1, l1 = ax1.get_legend_handles_labels()
    lines2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, l1 + l2, loc="upper right", fontsize=8)

    plt.tight_layout()
    plt.show()


def plot_position_cog_frequencies(positions: dict, title: str):
    """
    3-panel figure: COG letter frequencies for target / upstream / downstream.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"COG Category Frequencies by Position — {title}", fontsize=11)

    panel_titles = {
        "target":     "smORF Target",
        "upstream":   "Upstream Neighbour",
        "downstream": "Downstream Neighbour",
    }

    for ax, pos in zip(axes, ["target", "upstream", "downstream"]):
        top = positions[pos].most_common(15)
        if not top:
            ax.set_title(panel_titles[pos])
            continue
        ltrs  = [item[0] for item in top]
        cnts  = [item[1] for item in top]
        _cog_bar(ax, ltrs, cnts, panel_titles[pos])

    # shared colour-group legend
    legend_handles = [
        mpatches.Patch(color="#e07b54", label="Information storage & processing"),
        mpatches.Patch(color="#5b8db8", label="Cellular processes & signalling"),
        mpatches.Patch(color="#6aaa64", label="Metabolism"),
        mpatches.Patch(color="#aaaaaa", label="Poorly characterised"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=4,
               fontsize=8, bbox_to_anchor=(0.5, -0.04))

    plt.tight_layout()
    plt.show()


def plot_cog_accuracy_per_rule(cog_stats: dict, title: str):
    """
    4-panel figure — one panel per rule.
    Each panel: bars = eligible count per COG letter,
                line = accuracy % per COG letter.
    Sorted by accuracy descending.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"COG Category: Eligible Count & Accuracy per Rule — {title}",
                 fontsize=11)

    for ax, k in zip(axes.flat, RULE_KEYS):
        rows = cog_accuracy_table(
            cog_stats[k]["cog_eligible"],
            cog_stats[k]["cog_correct"],
        )
        if not rows:
            ax.set_title(RULE_LABELS[k])
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        ltrs  = [r[0] for r in rows]
        els   = [r[1] for r in rows]
        accs  = [r[3] for r in rows]
        cols  = [COG_META.get(l, ("", "#cccccc"))[1] for l in ltrs]

        ax2 = ax.twinx()
        ax.bar(range(len(ltrs)), els, color=cols, alpha=0.7,
               edgecolor="white", linewidth=0.5, label="Eligible")
        ax2.plot(range(len(ltrs)), accs, marker="o", color="#333333",
                 linewidth=1.8, label="Accuracy %", zorder=5)
        ax2.axhline(ACCURACY_TARGET, color="#333333", linestyle=":",
                    alpha=0.5, linewidth=1)
        ax2.set_ylim(0, 110)
        ax2.set_ylabel("Accuracy (%)", fontsize=8)

        ax.set_xticks(range(len(ltrs)))
        ax.set_xticklabels(ltrs, fontsize=9, fontweight="bold")
        ax.set_ylabel("Eligible count", fontsize=8)
        ax.set_title(RULE_LABELS[k], fontsize=9)
        ax.grid(True, axis="y", ls="--", alpha=0.3)

        # accuracy labels on line points
        for xi, acc in enumerate(accs):
            ax2.text(xi, acc + 2, f"{acc:.0f}%", ha="center",
                     fontsize=7, color="#333333")

        lines1, l1 = ax.get_legend_handles_labels()
        lines2, l2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, l1 + l2, fontsize=7, loc="upper right")

    plt.tight_layout()
    plt.show()


def plot_correct_cog_frequencies(cog_stats: dict, title: str):
    """
    4-panel figure — for each rule, bar chart of COG letters in CORRECT predictions.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(f"COG Categories in Correct Predictions per Rule — {title}",
                 fontsize=11)

    for ax, k in zip(axes.flat, RULE_KEYS):
        top = cog_stats[k]["cog_correct"].most_common(15)
        if not top:
            ax.set_title(RULE_LABELS[k])
            ax.text(0.5, 0.5, "No correct predictions", ha="center",
                    va="center", transform=ax.transAxes)
            continue
        ltrs = [item[0] for item in top]
        cnts = [item[1] for item in top]
        _cog_bar(ax, ltrs, cnts, RULE_LABELS[k])

    # shared legend
    legend_handles = [
        mpatches.Patch(color="#e07b54", label="Information storage & processing"),
        mpatches.Patch(color="#5b8db8", label="Cellular processes & signalling"),
        mpatches.Patch(color="#6aaa64", label="Metabolism"),
        mpatches.Patch(color="#aaaaaa", label="Poorly characterised"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=4,
               fontsize=8, bbox_to_anchor=(0.5, -0.03))

    plt.tight_layout()
    plt.show()


# ══════════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    for folder in RANGE_DIRS:
        root = BASE / folder

        print(f"\n{'═'*72}")
        print(f"  smORF Rule Evaluator v6  —  {folder}")
        print(f"{'═'*72}\n")

        # ── load evidence ─────────────────────────────────────────────────
        all_evidence = []
        for d in root.iterdir():
            if d.is_dir() and d.name.startswith("SHD1_SM.100AA"):
                all_evidence.append(aggregate_smorf_evidence(d))

        total_sm     = len(all_evidence)
        ann_evidence = [ev for ev in all_evidence if is_majority_annotated(ev)]
        unann_evid   = [ev for ev in all_evidence if not is_majority_annotated(ev)]
        ann_sm       = len(ann_evidence)

        print(f"  Total smORF clusters : {total_sm:,}")
        print(f"  Annotated            : {ann_sm:,}  ({ann_sm/total_sm*100:.1f}%)")
        print(f"  Unannotated          : {len(unann_evid):,}\n")

        # ══════════════════════════════════════════════════════════════════
        #  1. Rule table
        # ══════════════════════════════════════════════════════════════════
        counters = {k: {"eligible": 0, "correct": 0} for k in RULE_KEYS}
        for ev in ann_evidence:
            r = evaluate_rules(ev)
            for k in RULE_KEYS:
                if r[k]["eligible"]:
                    counters[k]["eligible"] += 1
                    counters[k]["correct"]  += int(r[k]["correct"])

        print(f"{'Rule':<35} {'Eligible':>9} {'Correct':>9} {'Cov%':>7} {'Acc%':>7}")
        print("─" * 63)
        for k in RULE_KEYS:
            el, co = counters[k]["eligible"], counters[k]["correct"]
            cov = el / ann_sm * 100 if ann_sm else 0
            acc = co / el * 100     if el     else 0
            flag = " ✓" if acc >= ACCURACY_TARGET else ""
            print(f"{RULE_LABELS[k]:<35} {el:>9} {co:>9} {cov:>7.1f} {acc:>7.1f}{flag}")

        # ══════════════════════════════════════════════════════════════════
        #  2. COG analysis
        # ══════════════════════════════════════════════════════════════════
        cog_stats, positions = build_cog_rule_stats(ann_evidence)

        print(f"\n{'═'*72}")
        print("  COG Category Analysis")
        print(f"{'═'*72}")

        # ── 2a. Position frequencies ──────────────────────────────────────
        for pos in ["target", "upstream", "downstream"]:
            label = {"target": "smORF target",
                     "upstream": "Upstream neighbour",
                     "downstream": "Downstream neighbour"}[pos]
            top = positions[pos].most_common(10)
            print(f"\n  {label} — top COG categories:")
            print(f"  {'Letter':<8} {'Count':>6}   Description")
            print("  " + "─" * 52)
            for letter, cnt in top:
                desc = COG_META.get(letter, ("Unknown", ""))[0]
                print(f"  {letter:<8} {cnt:>6}   {desc}")

        # ── 2b. Per-rule COG accuracy table ───────────────────────────────
        print(f"\n{'═'*72}")
        print(f"  COG Accuracy per Rule  (min {MIN_COG_COUNT} eligible)")
        print(f"{'═'*72}")

        for k in RULE_KEYS:
            rows = cog_accuracy_table(
                cog_stats[k]["cog_eligible"],
                cog_stats[k]["cog_correct"],
            )
            print(f"\n  {RULE_LABELS[k]}")
            print(f"  {'Letter':<8} {'Eligible':>9} {'Correct':>9} {'Acc%':>7}   Description")
            print("  " + "─" * 65)
            if not rows:
                print("  (no COG letters with enough data)")
                continue
            for letter, el, co, acc in rows:
                desc = COG_META.get(letter, ("Unknown", ""))[0]
                flag = " ✓" if acc >= ACCURACY_TARGET else ""
                print(f"  {letter:<8} {el:>9} {co:>9} {acc:>7.1f}{flag}   {desc}")

        # ── 2c. Most common COGs in correct predictions ───────────────────
        print(f"\n{'═'*72}")
        print("  Most common COG categories in CORRECT predictions")
        print(f"{'═'*72}")
        for k in RULE_KEYS:
            top3 = cog_stats[k]["cog_correct"].most_common(5)
            top3_str = "  ".join(
                f"{l}={n}" for l, n in top3
            ) if top3 else "—"
            print(f"  {RULE_LABELS[k]:<35}: {top3_str}")

        # ══════════════════════════════════════════════════════════════════
        #  3. Unannotated summary
        # ══════════════════════════════════════════════════════════════════
        unann_any = unann_up = unann_down = unann_both_same = unann_both_diff = 0
        pred_methods: Counter = Counter()

        for ev in unann_evid:
            has_up, has_down = bool(ev["upstream"]), bool(ev["downstream"])
            if has_up or has_down:
                unann_any += 1
                if has_up:   unann_up   += 1
                if has_down: unann_down += 1
                if has_up and has_down:
                    if ev["upstream"] == ev["downstream"]: unann_both_same += 1
                    else:                                   unann_both_diff += 1

            score_map: Counter = Counter()
            for dist in range(1, N_NEIGHBOURS + 1):
                w = 1.0 / dist
                for cog, cnt in ev["multi_up_votes"].get(dist, {}).items():
                    score_map[cog] += cnt * w
                for cog, cnt in ev["multi_down_votes"].get(dist, {}).items():
                    score_map[cog] += cnt * w
            if score_map:
                uc, dc = ev["upstream"], ev["downstream"]
                if   is_annotated(uc) and is_annotated(dc) and uc == dc: method = "both_agree"
                elif is_annotated(uc) and not is_annotated(dc):           method = "upstream_only"
                elif is_annotated(dc) and not is_annotated(uc):           method = "downstream_only"
                else:                                                      method = "weighted_vote"
                pred_methods[method] += 1

        print(f"\n{'═'*72}")
        print("  Unannotated smORFs")
        print(f"{'═'*72}")
        print(f"  Total:                          {len(unann_evid):,}")
        print(f"  Has any annotated neighbour:    {unann_any:,}")
        print(f"  Upstream annotated:             {unann_up:,}")
        print(f"  Downstream annotated:           {unann_down:,}")
        print(f"  Both same function:             {unann_both_same:,}")
        print(f"  Both annotated but different:   {unann_both_diff:,}  ← operon-boundary candidates")
        print(f"\n  Prediction methods (frequency-weighted):")
        for method, cnt in pred_methods.most_common():
            print(f"    • {method:<25}: {cnt:,}")

        # ══════════════════════════════════════════════════════════════════
        #  4. Plots
        # ══════════════════════════════════════════════════════════════════
        print(f"\n{'─'*72}")
        print("  Generating plots …\n")
        plot_rule_overview(counters, ann_sm, folder)
        plot_position_cog_frequencies(positions, folder)
        plot_cog_accuracy_per_rule(cog_stats, folder)
        plot_correct_cog_frequencies(cog_stats, folder)


if __name__ == "__main__":
    main()
