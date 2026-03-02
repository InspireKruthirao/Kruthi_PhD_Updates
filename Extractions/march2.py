#!/usr/bin/env python3
"""
smORF Neighbourhood Analysis  -  v4
====================================
Folder analysed: SmORF_neighbourhoods_26_30_no_min_overlap

New in v4 (supervisor feedback - Luis Pedro):
  1. Right vs Left neighbour match comparison
       Does the smORF share function more with RIGHT (downstream) or
       LEFT (upstream) neighbour?  Confirms/investigates R3 > R2 observation.

  2. Strand consistency check
       After strand-normalisation (flip so upstream=left), what fraction of
       smORFs have all three genes on the SAME strand?

  3. Gene orientation / synteny context
       Classify each neighbourhood by arrow pattern of (left, target, right):
         ">>>"  all same strand as target
         ">>>"  variations: "->>"  "<>>"  ">><"  etc.
       Measure Coverage & Accuracy per orientation pattern to see if
       co-directional genes predict function better.

Core principles (unchanged from v3):
  - One unique smORF directory = one smORF (no occurrence inflation)
  - Majority-vote across all occurrences for representative COG/strand
  - Coverage = eligible / total annotated smORFs
  - Accuracy  = predicted_correct / eligible
"""

from pathlib import Path
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
#  Configuration
# ─────────────────────────────────────────────────────────────────────────────
BASE       = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = [
    "SmORF_neighbourhoods_16_20_no_min_overlap",
    "SmORF_neighbourhoods_21_25_no_min_overlap",
    "SmORF_neighbourhoods_26_30_no_min_overlap",
]

FEATURE_TYPE    = "CDS"
UNKNOWN         = "Unknown"
ACCURACY_TARGET = 80.0

RULE_KEYS = ["R1_both_same", "R2_left_match", "R3_right_match", "R4_either_match"]

RULE_LABELS = {
    "R1_both_same":    "R1: Both neighbours\nsame function",
    "R2_left_match":   "R2: Left (upstream)\nneighbour",
    "R3_right_match":  "R3: Right (downstream)\nneighbour",
    "R4_either_match": "R4: Either neighbour\nmatches",
}

COLORS = {
    "R1_both_same":    "#2196F3",
    "R2_left_match":   "#4CAF50",
    "R3_right_match":  "#FF9800",
    "R4_either_match": "#9C27B0",
}

# Arrow symbols for orientation patterns
ARROW = {"+": "→", "-": "←", None: "?"}

# ─────────────────────────────────────────────────────────────────────────────
#  Low-level utilities
# ─────────────────────────────────────────────────────────────────────────────

def parse_attrs(s: str) -> dict:
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def cog_set(name: str) -> frozenset:
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return frozenset()
    i = str(name).rfind("-")
    if i == -1:
        return frozenset()
    return frozenset(ch for ch in str(name)[i + 1:] if ch.isalpha())


def is_annotated(cog: frozenset) -> bool:
    return bool(cog)


def majority(counter: Counter):
    """Return most common element or None if counter is empty."""
    if not counter:
        return None
    return counter.most_common(1)[0][0]


def orientation_label(left_strand, target_strand, right_strand):
    """
    Build a 3-character orientation string using arrows.
    All strands are expressed RELATIVE to the target after normalisation,
    so target is always '+' after flipping.
    Example: left='+', target='+', right='-'  →  '→→←'
    """
    def sym(s):
        return ARROW.get(s, "?")
    return f"{sym(left_strand)}{sym(target_strand)}{sym(right_strand)}"


# ─────────────────────────────────────────────────────────────────────────────
#  GFF parsing  (strand-normalised)
# ─────────────────────────────────────────────────────────────────────────────

def read_gff_ordered(gff_file: Path):
    """
    Parse GFF, sort by genomic position, then flip the whole list if the
    target is on the minus strand so that:
      - index i-1 = upstream  (left)
      - index i+1 = downstream (right)

    Returns list of feature dicts:
      start, end, strand (original), norm_strand (after flip), is_target, cog
    """
    features = []
    with gff_file.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue
            a      = parse_attrs(parts[8])
            strand = parts[6]
            start  = int(parts[3])
            end    = int(parts[4])
            features.append(dict(
                start     = start,
                end       = end,
                strand    = strand,
                is_target = a.get("target", "0") == "1",
                cog       = cog_set(a.get("Name", UNKNOWN)),
            ))

    features.sort(key=lambda x: x["start"])

    # Strand-normalise: if target is on minus strand, flip list AND invert strands
    target_strands = [f["strand"] for f in features if f["is_target"]]
    flipped = bool(target_strands and target_strands[0] == "-")

    if flipped:
        features = list(reversed(features))
        # Invert strands so target reads as '+'
        flip_map = {"+": "-", "-": "+"}
        for f in features:
            f["norm_strand"] = flip_map.get(f["strand"], f["strand"])
    else:
        for f in features:
            f["norm_strand"] = f["strand"]

    return features, flipped


# ─────────────────────────────────────────────────────────────────────────────
#  Per-smORF evidence aggregation  (across ALL occurrences)
# ─────────────────────────────────────────────────────────────────────────────

def aggregate_smorf_evidence(smorf_dir: Path):
    """
    Aggregate majority-vote evidence across all GFF occurrences.

    Returns dict with:
      target_cog_votes   Counter  COG frozenset -> count
      left_cog_votes     Counter  COG frozenset -> count
      right_cog_votes    Counter  COG frozenset -> count
      left_strand_votes  Counter  strand char   -> count  (norm_strand)
      right_strand_votes Counter  strand char   -> count  (norm_strand)
      target_strand_votes Counter strand char   -> count  (norm_strand, always '+' after flip)
      n_annotated        int
      n_unannotated      int
      n_total            int
      n_flipped          int   how many occurrences were on minus strand
      n_all_same_strand  int   occurrences where left/target/right all same norm_strand
    """
    target_cog_votes    = Counter()
    left_cog_votes      = Counter()
    right_cog_votes     = Counter()
    left_strand_votes   = Counter()
    right_strand_votes  = Counter()
    target_strand_votes = Counter()
    n_annotated         = 0
    n_unannotated       = 0
    n_total             = 0
    n_flipped           = 0
    n_all_same_strand   = 0

    for contig in smorf_dir.iterdir():
        if not contig.is_dir():
            continue
        for gff in contig.glob("*.gff"):
            seq, flipped = read_gff_ordered(gff)
            if flipped:
                n_flipped += 1

            for i, feat in enumerate(seq):
                if not feat["is_target"]:
                    continue
                n_total += 1

                tcog      = feat["cog"]
                left_feat  = seq[i - 1] if i - 1 >= 0       else None
                right_feat = seq[i + 1] if i + 1 < len(seq) else None

                left_cog   = left_feat["cog"]   if left_feat  else frozenset()
                right_cog  = right_feat["cog"]  if right_feat else frozenset()
                left_str   = left_feat["norm_strand"]  if left_feat  else None
                right_str  = right_feat["norm_strand"] if right_feat else None
                target_str = feat["norm_strand"]

                # COG votes
                if is_annotated(tcog):
                    n_annotated += 1
                    target_cog_votes[tcog] += 1
                else:
                    n_unannotated += 1

                if is_annotated(left_cog):
                    left_cog_votes[left_cog] += 1
                if is_annotated(right_cog):
                    right_cog_votes[right_cog] += 1

                # Strand votes
                target_strand_votes[target_str] += 1
                if left_str:
                    left_strand_votes[left_str] += 1
                if right_str:
                    right_strand_votes[right_str] += 1

                # All-same-strand check
                if (left_str and right_str and
                        left_str == target_str == right_str):
                    n_all_same_strand += 1

    return dict(
        target_cog_votes    = target_cog_votes,
        left_cog_votes      = left_cog_votes,
        right_cog_votes     = right_cog_votes,
        left_strand_votes   = left_strand_votes,
        right_strand_votes  = right_strand_votes,
        target_strand_votes = target_strand_votes,
        n_annotated         = n_annotated,
        n_unannotated       = n_unannotated,
        n_total             = n_total,
        n_flipped           = n_flipped,
        n_all_same_strand   = n_all_same_strand,
    )


# ─────────────────────────────────────────────────────────────────────────────
#  Annotation status
# ─────────────────────────────────────────────────────────────────────────────

def smorf_annotation_status(evidence: dict):
    if evidence["n_total"] == 0:
        return None
    if evidence["n_annotated"] > 0 and evidence["n_unannotated"] == 0:
        return "all_annotated"
    if evidence["n_unannotated"] > 0 and evidence["n_annotated"] == 0:
        return "all_unannotated"
    return "mixed"


# ─────────────────────────────────────────────────────────────────────────────
#  Per-smORF rule evaluation
# ─────────────────────────────────────────────────────────────────────────────

def evaluate_rules_per_smorf(evidence: dict):
    """
    ONE decision per unique smORF using majority-vote representative COGs.
    Returns ({rule: {eligible, predicted_correct}}, is_ann_target)
    """
    target_cog = majority(evidence["target_cog_votes"]) or frozenset()
    left_cog   = majority(evidence["left_cog_votes"])   or frozenset()
    right_cog  = majority(evidence["right_cog_votes"])  or frozenset()
    is_ann_target = evidence["n_annotated"] > evidence["n_unannotated"]

    results = {}

    # R1
    r1_eligible  = is_annotated(left_cog) and is_annotated(right_cog) and left_cog == right_cog
    r1_predicted = r1_eligible and is_ann_target and target_cog == left_cog
    results["R1_both_same"] = {"eligible": r1_eligible, "predicted_correct": r1_predicted}

    # R2
    r2_eligible  = is_annotated(left_cog)
    r2_predicted = r2_eligible and is_ann_target and target_cog == left_cog
    results["R2_left_match"] = {"eligible": r2_eligible, "predicted_correct": r2_predicted}

    # R3
    r3_eligible  = is_annotated(right_cog)
    r3_predicted = r3_eligible and is_ann_target and target_cog == right_cog
    results["R3_right_match"] = {"eligible": r3_eligible, "predicted_correct": r3_predicted}

    # R4
    r4_eligible  = is_annotated(left_cog) or is_annotated(right_cog)
    if r4_eligible and is_ann_target:
        r4_predicted = (is_annotated(left_cog)  and target_cog == left_cog) or \
                       (is_annotated(right_cog) and target_cog == right_cog)
    else:
        r4_predicted = False
    results["R4_either_match"] = {"eligible": r4_eligible, "predicted_correct": r4_predicted}

    return results, is_ann_target


# ─────────────────────────────────────────────────────────────────────────────
#  NEW: Right vs Left neighbour match analysis
# ─────────────────────────────────────────────────────────────────────────────

def analyse_right_vs_left(evidence: dict, is_ann: bool):
    """
    For an annotated smORF, check whether:
      - target COG matches LEFT neighbour only
      - target COG matches RIGHT neighbour only
      - target COG matches BOTH neighbours
      - target COG matches NEITHER
    Returns a category string or None if not annotated.
    """
    if not is_ann:
        return None

    target_cog = majority(evidence["target_cog_votes"]) or frozenset()
    left_cog   = majority(evidence["left_cog_votes"])   or frozenset()
    right_cog  = majority(evidence["right_cog_votes"])  or frozenset()

    left_match  = is_annotated(left_cog)  and target_cog == left_cog
    right_match = is_annotated(right_cog) and target_cog == right_cog

    if left_match and right_match:
        return "matches_both"
    elif right_match:
        return "matches_right_only"
    elif left_match:
        return "matches_left_only"
    else:
        return "matches_neither"


# ─────────────────────────────────────────────────────────────────────────────
#  NEW: Strand consistency & orientation pattern
# ─────────────────────────────────────────────────────────────────────────────

def get_orientation_pattern(evidence: dict):
    """
    Determine majority-vote orientation pattern for this smORF.
    After normalisation, target is always '+'.
    Pattern = arrow(left) + arrow(target) + arrow(right)

    Also returns:
      is_same_strand : bool  all three genes share the same norm_strand
    """
    left_strand   = majority(evidence["left_strand_votes"])
    right_strand  = majority(evidence["right_strand_votes"])
    target_strand = majority(evidence["target_strand_votes"])  # should be '+'

    pattern = orientation_label(left_strand, target_strand, right_strand)

    # All-same-strand: after normalisation, target is '+', so same means all '+'
    is_same_strand = (left_strand == target_strand == right_strand == "+")

    return pattern, is_same_strand


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

def plot_bar(data_dict, title, ylabel="Number of smORFs", color=None):
    labels = list(data_dict.keys())
    values = list(data_dict.values())
    colors = (color if color else plt.cm.Dark2(np.linspace(0, 1, len(labels))))
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.4), 4.5))
    bars = ax.bar(labels, values, color=colors, width=0.5,
                  edgecolor="white", linewidth=0.8)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(axis="x", labelrotation=30, labelsize=9)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    vmax = max(values) if values else 1
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + vmax * 0.012,
                str(val), ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.show()


def plot_right_vs_left(counts, title):
    """
    Horizontal bar chart showing how many annotated smORFs match
    right-only, left-only, both, or neither.
    """
    order  = ["matches_both", "matches_right_only", "matches_left_only", "matches_neither"]
    labels = ["Matches both\nneighbours", "Matches right\n(downstream) only",
              "Matches left\n(upstream) only", "Matches neither\nneighbour"]
    values = [counts.get(k, 0) for k in order]
    colors = ["#2196F3", "#FF9800", "#4CAF50", "#9E9E9E"]
    total  = sum(values)

    fig, ax = plt.subplots(figsize=(9, 4))
    bars = ax.barh(labels, values, color=colors, edgecolor="white", height=0.5)
    ax.set_xlabel("Number of annotated smORFs", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    for bar, val in zip(bars, values):
        pct = val / total * 100 if total else 0
        ax.text(bar.get_width() + total * 0.005, bar.get_y() + bar.get_height() / 2,
                f"{val}  ({pct:.1f}%)", va="center", fontsize=9)
    ax.set_xlim(0, max(values) * 1.25)
    plt.tight_layout()
    plt.show()


def plot_orientation_accuracy(orientation_stats, title):
    """
    For each orientation pattern, show:
      - how many smORFs have that pattern
      - R1 accuracy within that pattern
    Helps answer: do co-directional genes (→→→) have higher accuracy?
    """
    # Filter to patterns with at least 3 smORFs for reliability
    patterns = [p for p, s in orientation_stats.items()
                if s["total"] >= 3]
    patterns.sort(key=lambda p: orientation_stats[p]["total"], reverse=True)

    if not patterns:
        print("  [!] Not enough data for orientation accuracy plot.")
        return

    totals   = [orientation_stats[p]["total"]   for p in patterns]
    r1_acc   = [orientation_stats[p]["r1_correct"] / orientation_stats[p]["r1_eligible"] * 100
                if orientation_stats[p]["r1_eligible"] > 0 else 0
                for p in patterns]

    x, w = np.arange(len(patterns)), 0.35
    fig, ax1 = plt.subplots(figsize=(max(8, len(patterns) * 1.2), 5))
    ax2 = ax1.twinx()

    b1 = ax1.bar(x - w/2, totals, w, label="# smORFs", color="#90CAF9", edgecolor="white")
    b2 = ax2.bar(x + w/2, r1_acc,  w, label="R1 Accuracy (%)", color="#FF8F00", edgecolor="white")
    ax2.axhline(ACCURACY_TARGET, color="#E53935", linestyle="--", linewidth=1.3,
                label=f"{ACCURACY_TARGET}% target")

    ax1.set_xticks(x)
    ax1.set_xticklabels(patterns, fontsize=11)
    ax1.set_ylabel("Number of smORFs", fontsize=10)
    ax2.set_ylabel("R1 Accuracy (%)", fontsize=10)
    ax2.set_ylim(0, 115)
    ax1.set_title(title, fontsize=12, fontweight="bold")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc="upper right")
    ax1.grid(axis="y", linestyle="--", alpha=0.25)
    plt.tight_layout()
    plt.show()


def plot_strand_consistency(same_strand_count, diff_strand_count, title):
    """
    Pie chart: what fraction of smORFs have all three genes on the same strand?
    """
    sizes  = [same_strand_count, diff_strand_count]
    labels = [f"All same strand\n({same_strand_count})",
              f"Mixed strands\n({diff_strand_count})"]
    colors = ["#4CAF50", "#FF9800"]
    fig, ax = plt.subplots(figsize=(5, 5))
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, colors=colors,
        autopct="%1.1f%%", startangle=90,
        wedgeprops=dict(edgecolor="white", linewidth=1.5)
    )
    for t in autotexts:
        t.set_fontsize(11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.show()


def plot_coverage_accuracy_bars(rule_counters, total_smorf, title):
    coverages, accuracies = [], []
    print(f"\n{'─'*72}")
    print(f"  Rule Evaluation (per unique smORF)  --  {title}")
    print(f"{'─'*72}")
    print(f"  {'Rule':<38} {'Eligible':>9} {'Predicted':>10} {'Cov%':>8} {'Acc%':>8}")
    print(f"{'─'*72}")
    for r in RULE_KEYS:
        eligible          = rule_counters[r]["eligible"]
        predicted_correct = rule_counters[r]["predicted_correct"]
        cov = eligible          / total_smorf * 100 if total_smorf else 0
        acc = predicted_correct / eligible    * 100 if eligible    else 0
        coverages.append(cov)
        accuracies.append(acc)
        flag = "  <<< hits target" if acc >= ACCURACY_TARGET else ""
        print(f"  {RULE_LABELS[r].replace(chr(10),' '):<38} "
              f"{eligible:>9} {predicted_correct:>10} {cov:>8.1f} {acc:>8.1f}{flag}")
    print(f"{'─'*72}")
    print(f"  Total unique annotated smORFs: {total_smorf}")

    labels = [RULE_LABELS[r] for r in RULE_KEYS]
    x, w = np.arange(len(RULE_KEYS)), 0.35
    fig, ax = plt.subplots(figsize=(10, 5))
    b1 = ax.bar(x - w/2, coverages,  w, label="Coverage (%)",
                color="#2196F3", alpha=0.85, edgecolor="white")
    b2 = ax.bar(x + w/2, accuracies, w, label="Accuracy (%)",
                color="#4CAF50", alpha=0.85, edgecolor="white")
    ax.axhline(ACCURACY_TARGET, color="#E53935", linestyle="--", linewidth=1.5,
               label=f"{ACCURACY_TARGET}% accuracy target")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=9)
    ax.set_ylabel("Percentage (%)", fontsize=10)
    ax.set_ylim(0, 110)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    for bar in list(b1) + list(b2):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=7.5)
    plt.tight_layout()
    plt.show()


def plot_pareto_frontier(rule_counters, total_smorf, title):
    points = {}
    for r in RULE_KEYS:
        eligible          = rule_counters[r]["eligible"]
        predicted_correct = rule_counters[r]["predicted_correct"]
        cov = eligible          / total_smorf * 100 if total_smorf else 0
        acc = predicted_correct / eligible    * 100 if eligible    else 0
        points[r] = (cov, acc)

    sorted_pts = sorted(points.values(), key=lambda p: p[0])
    frontier, best_acc = [], -1
    for cov, acc in sorted_pts:
        if acc >= best_acc:
            frontier.append((cov, acc))
            best_acc = acc

    fig, ax = plt.subplots(figsize=(7, 5))
    if frontier:
        fx, fy = zip(*frontier)
        ax.fill_between([0] + list(fx), [0] + list(fy), 0,
                        alpha=0.07, color="#2196F3", label="Dominated region (avoid)")
        if len(frontier) >= 2:
            ax.plot(fx, fy, "--", color="#E53935", linewidth=1.8,
                    label="Pareto frontier", zorder=3)

    on_frontier = set(frontier)
    for r, (cov, acc) in points.items():
        ec = "black" if (cov, acc) in on_frontier else "none"
        ax.scatter(cov, acc, s=130, color=COLORS[r], edgecolors=ec,
                   linewidths=1.8, zorder=5)
        ax.annotate(RULE_LABELS[r].replace("\n", " "),
                    (cov, acc), textcoords="offset points",
                    xytext=(8, 4), fontsize=8.5)

    ax.axhline(ACCURACY_TARGET, color="#FF9800", linestyle=":",
               linewidth=1.5, label=f"{ACCURACY_TARGET}% accuracy target")
    ax.set_xlabel("Coverage (% of unique annotated smORFs)", fontsize=10)
    ax.set_ylabel("Accuracy (%)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(-2, 105)
    ax.set_ylim(-2, 105)
    ax.legend(fontsize=8.5, loc="lower left")
    ax.grid(linestyle="--", alpha=0.35)
    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    for d in RANGE_DIRS:
        ROOT = BASE / d

        total_smorfs      = 0
        annotation_counts = {"all_annotated": 0, "all_unannotated": 0, "mixed": 0}

        # Rule counters (per unique smORF)
        global_rule_counters = {r: {"eligible": 0, "predicted_correct": 0} for r in RULE_KEYS}
        global_total_ann     = 0

        # Unannotated neighbour stats
        unan_neighbour_stats = defaultdict(int)

        # ── NEW: Right vs Left match counts ──────────────────────────────────
        right_vs_left_counts = defaultdict(int)

        # ── NEW: Strand consistency ───────────────────────────────────────────
        strand_same_count = 0   # all 3 genes same strand (after normalisation)
        strand_diff_count = 0   # mixed strands
        total_flipped     = 0   # how many smORFs were on minus strand

        # ── NEW: Orientation pattern stats ───────────────────────────────────
        # pattern -> {total, r1_eligible, r1_correct}
        orientation_stats = defaultdict(lambda: {"total": 0,
                                                  "r1_eligible": 0,
                                                  "r1_correct": 0})

        # ── Main loop ─────────────────────────────────────────────────────────
        for smorf_dir in ROOT.iterdir():
            if not smorf_dir.is_dir() or not smorf_dir.name.startswith("SHD1_SM.100AA"):
                continue

            total_smorfs += 1

            evidence = aggregate_smorf_evidence(smorf_dir)

            # Annotation status
            status = smorf_annotation_status(evidence)
            if status:
                annotation_counts[status] += 1

            # Rule evaluation
            rule_results, is_ann = evaluate_rules_per_smorf(evidence)

            if is_ann:
                global_total_ann += 1
                for r in RULE_KEYS:
                    if rule_results[r]["eligible"]:
                        global_rule_counters[r]["eligible"]          += 1
                        global_rule_counters[r]["predicted_correct"] += int(
                            rule_results[r]["predicted_correct"])

                # ── Right vs Left match ───────────────────────────────────────
                rvl = analyse_right_vs_left(evidence, is_ann)
                if rvl:
                    right_vs_left_counts[rvl] += 1

            else:
                left_cog  = majority(evidence["left_cog_votes"])  or frozenset()
                right_cog = majority(evidence["right_cog_votes"]) or frozenset()
                if is_annotated(left_cog) or is_annotated(right_cog):
                    unan_neighbour_stats["any_annotated_neighbor"] += 1
                    if is_annotated(left_cog):
                        unan_neighbour_stats["left_neighbor_annotated"] += 1
                    if is_annotated(right_cog):
                        unan_neighbour_stats["right_neighbor_annotated"] += 1
                    if is_annotated(left_cog) and is_annotated(right_cog):
                        if left_cog == right_cog:
                            unan_neighbour_stats["both_neighbors_same_function"] += 1
                        else:
                            unan_neighbour_stats["both_neighbors_different"] += 1

            # ── Strand consistency (all smORFs) ──────────────────────────────
            if evidence["n_flipped"] > evidence["n_total"] / 2:
                total_flipped += 1

            pattern, is_same = get_orientation_pattern(evidence)
            if is_same:
                strand_same_count += 1
            else:
                strand_diff_count += 1

            # ── Orientation pattern accuracy (annotated smORFs only) ─────────
            if is_ann:
                orientation_stats[pattern]["total"] += 1
                if rule_results["R1_both_same"]["eligible"]:
                    orientation_stats[pattern]["r1_eligible"] += 1
                    orientation_stats[pattern]["r1_correct"]  += int(
                        rule_results["R1_both_same"]["predicted_correct"])

        # ─── Print summary ─────────────────────────────────────────────────────
        print(f"\n{'='*70}")
        print(f"  Folder : {d}")
        print(f"{'='*70}")
        print(f"  Total unique smORFs  : {total_smorfs}")
        print(f"  Unique annotated     : {global_total_ann}")
        print(f"  Unique unannotated   : {total_smorfs - global_total_ann}")

        print("\n  Annotation status (per unique smORF, majority-vote):")
        for k, v in annotation_counts.items():
            pct = v / total_smorfs * 100 if total_smorfs else 0
            print(f"    {k:<20}: {v:>6}  ({pct:.1f}%)")

        print(f"\n  Rule evaluation  [denominator = {global_total_ann} annotated smORFs]")
        print(f"  {'Rule':<38} {'Eligible':>9} {'Predicted':>10} {'Cov%':>7} {'Acc%':>7}")
        print(f"  {'-'*68}")
        for r in RULE_KEYS:
            eligible          = global_rule_counters[r]["eligible"]
            predicted_correct = global_rule_counters[r]["predicted_correct"]
            cov = eligible          / global_total_ann * 100 if global_total_ann else 0
            acc = predicted_correct / eligible          * 100 if eligible         else 0
            flag = "  <<< hits 80% target" if acc >= ACCURACY_TARGET else ""
            print(f"  {RULE_LABELS[r].replace(chr(10),' '):<38} "
                  f"{eligible:>9} {predicted_correct:>10} {cov:>7.1f} {acc:>7.1f}{flag}")

        # ── NEW: Right vs Left summary ──────────────────────────────────────
        print(f"\n  {'─'*60}")
        print(f"  Right vs Left Neighbour Match  (annotated smORFs = {global_total_ann})")
        print(f"  {'─'*60}")
        rvl_order = ["matches_both", "matches_right_only",
                     "matches_left_only", "matches_neither"]
        for k in rvl_order:
            v   = right_vs_left_counts.get(k, 0)
            pct = v / global_total_ann * 100 if global_total_ann else 0
            print(f"    {k:<30}: {v:>5}  ({pct:.1f}%)")

        right_only = right_vs_left_counts.get("matches_right_only", 0)
        left_only  = right_vs_left_counts.get("matches_left_only", 0)
        if right_only > left_only:
            print(f"\n  >> RIGHT (downstream) neighbour matches more often "
                  f"({right_only} vs {left_only})  -- confirms R3 > R2 observation")
        elif left_only > right_only:
            print(f"\n  >> LEFT (upstream) neighbour matches more often "
                  f"({left_only} vs {right_only})")
        else:
            print(f"\n  >> Left and Right match equally ({left_only} each)")

        # ── NEW: Strand consistency summary ────────────────────────────────
        print(f"\n  {'─'*60}")
        print(f"  Strand Consistency  (all {total_smorfs} unique smORFs, after normalisation)")
        print(f"  {'─'*60}")
        pct_same = strand_same_count / total_smorfs * 100 if total_smorfs else 0
        pct_flip = total_flipped     / total_smorfs * 100 if total_smorfs else 0
        print(f"    All 3 genes same strand : {strand_same_count:>5}  ({pct_same:.1f}%)")
        print(f"    Mixed strands           : {strand_diff_count:>5}  ({100-pct_same:.1f}%)")
        print(f"    smORFs on minus strand  : {total_flipped:>5}  ({pct_flip:.1f}%)  [were flipped]")

        # ── NEW: Orientation pattern summary ───────────────────────────────
        print(f"\n  {'─'*60}")
        print(f"  Orientation Patterns  (annotated smORFs, sorted by frequency)")
        print(f"  {'─'*60}")
        print(f"  {'Pattern':<10} {'# smORFs':>9} {'R1 Eligible':>12} "
              f"{'R1 Correct':>11} {'R1 Acc%':>9}")
        sorted_patterns = sorted(orientation_stats.items(),
                                 key=lambda x: x[1]["total"], reverse=True)
        for pat, s in sorted_patterns:
            r1_acc = (s["r1_correct"] / s["r1_eligible"] * 100
                      if s["r1_eligible"] > 0 else 0)
            print(f"  {pat:<10} {s['total']:>9} {s['r1_eligible']:>12} "
                  f"{s['r1_correct']:>11} {r1_acc:>9.1f}%")

        # ─── Plots ──────────────────────────────────────────────────────────
        plot_bar(annotation_counts,
                 f"smORF Annotation Status  -  {d}")

        plot_bar(dict(unan_neighbour_stats),
                 f"Unannotated smORFs: Neighbour Potential  -  {d}")

        plot_coverage_accuracy_bars(
            global_rule_counters, global_total_ann,
            f"Coverage & Accuracy per Rule  -  {d}")

        plot_pareto_frontier(
            global_rule_counters, global_total_ann,
            f"Pareto Frontier: Coverage vs Accuracy  -  {d}")

        # NEW plots
        plot_right_vs_left(
            right_vs_left_counts,
            f"Right vs Left Neighbour Match  -  {d}")

        plot_strand_consistency(
            strand_same_count, strand_diff_count,
            f"Strand Consistency After Normalisation  -  {d}")

        plot_orientation_accuracy(
            dict(orientation_stats),
            f"R1 Accuracy by Gene Orientation Pattern  -  {d}")


if __name__ == "__main__":
    main()
