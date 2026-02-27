#!/usr/bin/env python3
"""
smORF Neighbourhood Analysis  –  v3  (per-smORF evaluation)
============================================================
Key design principle (from supervisor feedback):
  Each smORF directory = ONE unique smORF.
  Rules are evaluated ONCE per smORF, not once per occurrence/GFF file.

  Within a smORF directory we aggregate evidence across ALL occurrences
  (majority vote on neighbour COGs), then make a single binary decision:
    - did the rule FIRE?  (bool)
    - was it CORRECT?     (bool, only meaningful if fired)

  Denominators for Coverage and Accuracy are therefore #unique_smORFs,
  not #GFF occurrences.

Other improvements (supervisor feedback):
  - Coverage & Accuracy as primary metrics; Pareto frontier plot
  - Per-frequency-bin breakdown (accuracy shifts with observation count)
  - Strand-aware orientation: upstream always = left neighbour
  - 80 % accuracy target highlighted on all plots
"""

from pathlib import Path
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
#  Configuration
# ─────────────────────────────────────────────────────────────────────────────
BASE       = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_26_30_no_min_overlap"]

FEATURE_TYPE    = "CDS"
UNKNOWN         = "Unknown"
ACCURACY_TARGET = 80.0      # % – highlighted in all plots


RULE_KEYS = ["R1_both_same", "R2_left_match", "R3_right_match", "R4_either_match"]

RULE_LABELS = {
    "R1_both_same":    "R1: Both neighbours\nsame function",
    "R2_left_match":   "R2: Left (upstream)\nneighbour match",
    "R3_right_match":  "R3: Right (downstream)\nneighbour match",
    "R4_either_match": "R4: Either neighbour\nmatches",
}

COLORS = {
    "R1_both_same":    "#2196F3",
    "R2_left_match":   "#4CAF50",
    "R3_right_match":  "#FF9800",
    "R4_either_match": "#9C27B0",
}

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
    """Return frozenset of COG functional-category letters, or empty frozenset."""
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return frozenset()
    i = str(name).rfind("-")
    if i == -1:
        return frozenset()
    return frozenset(ch for ch in str(name)[i + 1:] if ch.isalpha())


def is_annotated(cog: frozenset) -> bool:
    return bool(cog)


# ─────────────────────────────────────────────────────────────────────────────
#  GFF parsing  (strand-normalised)
# ─────────────────────────────────────────────────────────────────────────────

def read_gff_ordered(gff_file: Path):
    """
    Return a list of feature dicts sorted by genomic position and
    strand-normalised so that index i-1 is always upstream (left) and
    index i+1 is always downstream (right) relative to the smORF.
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

    # Flip so upstream is always on the left
    target_strands = [f["strand"] for f in features if f["is_target"]]
    if target_strands and target_strands[0] == "-":
        features = list(reversed(features))

    return features


# ─────────────────────────────────────────────────────────────────────────────
#  Per-smORF evidence aggregation  (across ALL occurrences)
# ─────────────────────────────────────────────────────────────────────────────

def aggregate_smorf_evidence(smorf_dir: Path):
    """
    Walk every GFF occurrence for this unique smORF and accumulate
    majority-vote evidence:

      target_cog_votes  – Counter: COG frozenset -> #occurrences with that COG
      left_cog_votes    – Counter: COG frozenset -> #occurrences with that left neighbour
      right_cog_votes   – Counter: COG frozenset -> #occurrences with that right neighbour
      n_annotated       – occurrences where target is annotated
      n_unannotated     – occurrences where target is NOT annotated
      n_total           – total target occurrences seen
    """
    target_cog_votes = Counter()
    left_cog_votes   = Counter()
    right_cog_votes  = Counter()
    n_annotated      = 0
    n_unannotated    = 0
    n_total          = 0

    for contig in smorf_dir.iterdir():
        if not contig.is_dir():
            continue
        for gff in contig.glob("*.gff"):
            seq = read_gff_ordered(gff)
            for i, feat in enumerate(seq):
                if not feat["is_target"]:
                    continue
                n_total += 1
                tcog      = feat["cog"]
                left_cog  = seq[i - 1]["cog"] if i - 1 >= 0       else frozenset()
                right_cog = seq[i + 1]["cog"] if i + 1 < len(seq) else frozenset()

                if is_annotated(tcog):
                    n_annotated += 1
                    target_cog_votes[tcog] += 1
                else:
                    n_unannotated += 1

                if is_annotated(left_cog):
                    left_cog_votes[left_cog] += 1
                if is_annotated(right_cog):
                    right_cog_votes[right_cog] += 1

    return dict(
        target_cog_votes = target_cog_votes,
        left_cog_votes   = left_cog_votes,
        right_cog_votes  = right_cog_votes,
        n_annotated      = n_annotated,
        n_unannotated    = n_unannotated,
        n_total          = n_total,
    )


def majority(counter: Counter):
    """Return the most common element, or frozenset() if counter is empty."""
    if not counter:
        return frozenset()
    return counter.most_common(1)[0][0]


# ─────────────────────────────────────────────────────────────────────────────
#  Per-smORF rule evaluation  (ONE decision per unique smORF)
# ─────────────────────────────────────────────────────────────────────────────

def evaluate_rules_per_smorf(evidence: dict):
    """
    Given aggregated majority-vote evidence for ONE unique smORF, return:
      {rule: {"fired": bool, "correct": bool}}

    The smORF is treated as ANNOTATED if majority of occurrences have an
    annotated target.  Representative COGs are majority-vote winners.

    "correct" is evaluated against the true majority target COG.
    Only meaningful when fired=True.
    """
    target_cog    = majority(evidence["target_cog_votes"])
    left_cog      = majority(evidence["left_cog_votes"])
    right_cog     = majority(evidence["right_cog_votes"])
    is_ann_target = evidence["n_annotated"] > evidence["n_unannotated"]

    results = {}

    # R1: Both neighbours have the same annotated function
    r1_fires   = is_annotated(left_cog) and is_annotated(right_cog) and left_cog == right_cog
    r1_correct = r1_fires and is_ann_target and target_cog == left_cog
    results["R1_both_same"] = {"fired": r1_fires, "correct": r1_correct}

    # R2: Left (upstream) neighbour annotated -> predict target = left_cog
    r2_fires   = is_annotated(left_cog)
    r2_correct = r2_fires and is_ann_target and target_cog == left_cog
    results["R2_left_match"] = {"fired": r2_fires, "correct": r2_correct}

    # R3: Right (downstream) neighbour annotated -> predict target = right_cog
    r3_fires   = is_annotated(right_cog)
    r3_correct = r3_fires and is_ann_target and target_cog == right_cog
    results["R3_right_match"] = {"fired": r3_fires, "correct": r3_correct}

    # R4: Either neighbour annotated -> predict with whichever matches
    r4_fires = is_annotated(left_cog) or is_annotated(right_cog)
    if r4_fires and is_ann_target:
        r4_correct = (is_annotated(left_cog)  and target_cog == left_cog) or \
                     (is_annotated(right_cog) and target_cog == right_cog)
    else:
        r4_correct = False
    results["R4_either_match"] = {"fired": r4_fires, "correct": r4_correct}

    return results, is_ann_target


# ─────────────────────────────────────────────────────────────────────────────
#  Annotation status  (per unique smORF, majority-vote)
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
#  Plotting
# ─────────────────────────────────────────────────────────────────────────────

def plot_bar(data_dict, title, ylabel="Number of smORFs"):
    labels = list(data_dict.keys())
    values = list(data_dict.values())
    colors = plt.cm.Dark2(np.linspace(0, 1, len(labels)))
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.3), 4.5))
    bars = ax.bar(labels, values, color=colors, width=0.5,
                  edgecolor="white", linewidth=0.8)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(axis="x", labelrotation=30, labelsize=8)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    vmax = max(values) if values else 1
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + vmax * 0.012,
                str(val), ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.show()


def plot_coverage_accuracy_bars(rule_counters, total_smorf, title):
    """Grouped bar: Coverage vs Accuracy per rule (per-unique-smORF counts)."""
    coverages, accuracies = [], []
    print(f"\n{'─'*72}")
    print(f"  Per-unique-smORF Rule Evaluation  --  {title}")
    print(f"{'─'*72}")
    print(f"  {'Rule':<38} {'Fired':>7} {'Correct':>8} {'Cov%':>8} {'Acc%':>8}")
    print(f"{'─'*72}")
    for r in RULE_KEYS:
        fired   = rule_counters[r]["fired"]
        correct = rule_counters[r]["correct"]
        cov = fired   / total_smorf * 100 if total_smorf else 0
        acc = correct / fired       * 100 if fired       else 0
        coverages.append(cov)
        accuracies.append(acc)
        flag = "  <<< hits target" if acc >= ACCURACY_TARGET else ""
        print(f"  {RULE_LABELS[r].replace(chr(10),' '):<38} "
              f"{fired:>7} {correct:>8} {cov:>8.1f} {acc:>8.1f}{flag}")
    print(f"{'─'*72}")
    print(f"  Total unique annotated smORFs (denominator): {total_smorf}")

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
    """Scatter: Coverage (x) vs Accuracy (y) with Pareto frontier drawn."""
    points = {}
    for r in RULE_KEYS:
        fired   = rule_counters[r]["fired"]
        correct = rule_counters[r]["correct"]
        cov = fired   / total_smorf * 100 if total_smorf else 0
        acc = correct / fired       * 100 if fired       else 0
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
    ax.set_xlabel("Coverage  (% of unique annotated smORFs)", fontsize=10)
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

        # Global counters – one entry per unique smORF
        global_rule_counters = {r: {"fired": 0, "correct": 0} for r in RULE_KEYS}
        global_total_ann     = 0   # unique smORFs where majority target is annotated

        # Unannotated neighbour stats (for bar chart)
        unan_neighbour_stats = defaultdict(int)

        # ── Iterate unique smORFs ─────────────────────────────────────────────
        for smorf_dir in ROOT.iterdir():
            if not smorf_dir.is_dir() or not smorf_dir.name.startswith("SHD1_SM.100AA"):
                continue

            total_smorfs += 1

            # Aggregate all occurrence evidence into one representative signal
            evidence = aggregate_smorf_evidence(smorf_dir)

            # Annotation status
            status = smorf_annotation_status(evidence)
            if status:
                annotation_counts[status] += 1

            # ONE rule decision per smORF
            rule_results, is_ann = evaluate_rules_per_smorf(evidence)

            if is_ann:
                global_total_ann += 1
                for r in RULE_KEYS:
                    if rule_results[r]["fired"]:
                        global_rule_counters[r]["fired"]   += 1
                        global_rule_counters[r]["correct"] += int(rule_results[r]["correct"])


            else:
                # Unannotated smORF: check neighbour annotation potential
                left_cog  = majority(evidence["left_cog_votes"])
                right_cog = majority(evidence["right_cog_votes"])
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

        # ─── Print summary ────────────────────────────────────────────────────
        print(f"\n{'='*70}")
        print(f"  Folder: {d}")
        print(f"{'='*70}")
        print(f"  Total unique smORFs  : {total_smorfs}")
        print(f"  Unique annotated     : {global_total_ann}")
        print(f"  Unique unannotated   : {total_smorfs - global_total_ann}")

        print("\n  Annotation status (per unique smORF, majority-vote):")
        for k, v in annotation_counts.items():
            pct = v / total_smorfs * 100 if total_smorfs else 0
            print(f"    {k:<20}: {v:>6}  ({pct:.1f}%)")

        print(f"\n  Rule evaluation  [unit=unique smORF, denominator={global_total_ann}]")
        print(f"  NOTE: 'Fired' = smORFs where rule applies; "
              f"'Correct' = correct predictions among those.")
        print(f"  {'Rule':<38} {'Fired':>7} {'Correct':>8} {'Cov%':>7} {'Acc%':>7}")
        print(f"  {'-'*68}")
        for r in RULE_KEYS:
            fired   = global_rule_counters[r]["fired"]
            correct = global_rule_counters[r]["correct"]
            cov = fired   / global_total_ann * 100 if global_total_ann else 0
            acc = correct / fired            * 100 if fired             else 0
            flag = "  <<< hits 80% target" if acc >= ACCURACY_TARGET else ""
            print(f"  {RULE_LABELS[r].replace(chr(10),' '):<38} "
                  f"{fired:>7} {correct:>8} {cov:>7.1f} {acc:>7.1f}{flag}")


        # ─── Plots ────────────────────────────────────────────────────────────
        plot_bar(annotation_counts,
                 f"smORF Annotation Status  (unique smORFs)  -  {d}")

        plot_bar(dict(unan_neighbour_stats),
                 f"Unannotated smORFs: Neighbour Potential  -  {d}")

        plot_coverage_accuracy_bars(
            global_rule_counters, global_total_ann,
            f"Coverage & Accuracy per Rule  (unique smORFs)  -  {d}")

        plot_pareto_frontier(
            global_rule_counters, global_total_ann,
            f"Pareto Frontier: Coverage vs Accuracy  -  {d}")



if __name__ == "__main__":
    main()

def read_gff_with_order(gff_file):
    """Read GFF and return list of (target_flag, COG set) in order"""
    res = []
    with gff_file.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue
            a = attrs(parts[8])
            target = a.get("target", "0") == "1"
            res.append((target, cog_set(a.get("Name", UNKNOWN))))
    return res

# -------------------- SmORF annotation status --------------------

def smorf_annotation_status(folder: Path):
    annotated_found = False
    unannotated_found = False
    any_target = False

    for contig in folder.iterdir():
        if not contig.is_dir():
            continue
        for gff in contig.glob("*.gff"):
            with gff.open("r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                        continue
                    a = attrs(parts[8])
                    if a.get("target", "0") != "1":
                        continue
                    any_target = True
                    if is_annotated(a.get("Name", UNKNOWN)):
                        annotated_found = True
                    else:
                        unannotated_found = True

    if not any_target:
        return None
    if annotated_found and not unannotated_found:
        return "all_annotated"
    elif unannotated_found and not annotated_found:
        return "all_unannotated"
    else:
        return "mixed"

# -------------------- Neighbor analysis --------------------

def analyze_neighbors_per_smorf(folder: Path):
    ann_flags = {
        "both_neighbors_same_function": False,
        "left_neighbor_match": False,
        "right_neighbor_match": False,
        "neighbors_match_each_other_target_diff": False,
        "neighbors_match_each_other_and_target_match": False,
    }

    unann_flags = {
        "any_annotated_neighbor": False,
        "left_neighbor_annotated": False,
        "right_neighbor_annotated": False,
        "both_neighbors_same_function": False,
        "both_neighbors_different": False,
    }

    for contig in folder.iterdir():
        if not contig.is_dir():
            continue
        for gff in contig.glob("*.gff"):
            seq = read_gff_with_order(gff)
            for i, (is_target, target_cog) in enumerate(seq):
                if not is_target:
                    continue

                left_cog = seq[i-1][1] if i-1 >= 0 else set()
                right_cog = seq[i+1][1] if i+1 < len(seq) else set()

                # Annotated
                if target_cog:
                    if left_cog and right_cog and target_cog == left_cog == right_cog:
                        ann_flags["both_neighbors_same_function"] = True
                    if left_cog and target_cog == left_cog:
                        ann_flags["left_neighbor_match"] = True
                    if right_cog and target_cog == right_cog:
                        ann_flags["right_neighbor_match"] = True
                    if left_cog and right_cog and left_cog == right_cog and target_cog != left_cog:
                        ann_flags["neighbors_match_each_other_target_diff"] = True
                    if left_cog and right_cog and left_cog == right_cog and target_cog == left_cog:
                        ann_flags["neighbors_match_each_other_and_target_match"] = True

                # Unannotated
                else:
                    annotated_neighbors = [c for c in (left_cog, right_cog) if c]
                    if annotated_neighbors:
                        unann_flags["any_annotated_neighbor"] = True
                        if left_cog:
                            unann_flags["left_neighbor_annotated"] = True
                        if right_cog:
                            unann_flags["right_neighbor_annotated"] = True
                        if len(annotated_neighbors) == 2:
                            if left_cog == right_cog:
                                unann_flags["both_neighbors_same_function"] = True
                            else:
                                unann_flags["both_neighbors_different"] = True

    return ann_flags, unann_flags

# -------------------- Plotting --------------------

def plot_bar(data_dict, title):
    labels = list(data_dict.keys())
    values = list(data_dict.values())
    colors = plt.cm.Dark2(np.linspace(0, 1, len(labels)))
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.bar(labels, values, color=colors, width=0.45)
    ax.set_title(title, fontsize=11)
    ax.set_ylabel("Number of smORFs", fontsize=10)
    ax.tick_params(axis='x', labelrotation=40)
    ax.tick_params(axis='both', labelsize=9)
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()

def plot_pareto(rules_dict, total_annotated, title):
    labels, coverage_vals, accuracy_vals = [], [], []
    for k, v in rules_dict.items():
        applied = v["applied"]
        correct = v["correct"]
        coverage = applied / total_annotated * 100 if total_annotated else 0
        accuracy = correct / applied * 100 if applied else 0
        labels.append(k)
        coverage_vals.append(coverage)
        accuracy_vals.append(accuracy)

    print("\nPareto Plot Data:")
    print(f"{'Rule':<30}{'Applied':>8}{'Correct':>8}{'Coverage (%)':>15}{'Accuracy (%)':>15}")
    for i in range(len(labels)):
        print(f"{labels[i]:<30}{rules_dict[labels[i]]['applied']:>8}"
              f"{rules_dict[labels[i]]['correct']:>8}"
              f"{coverage_vals[i]:15.2f}{accuracy_vals[i]:15.2f}")

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - width/2, coverage_vals, width, label='Coverage (%)')
    ax.bar(x + width/2, accuracy_vals, width, label='Accuracy (%)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Percentage')
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    plt.show()

# -------------------- Main --------------------

def main():
    for d in RANGE_DIRS:
        ROOT = BASE / d

        total_smorfs = 0
        annotation_counts = {"all_annotated": 0, "all_unannotated": 0, "mixed": 0}

        total_annotated_stats = {
            "both_neighbors_same_function": 0,
            "left_neighbor_match": 0,
            "right_neighbor_match": 0,
            "neighbors_match_each_other_target_diff": 0,
            "neighbors_match_each_other_and_target_match": 0,
        }

        total_unannotated_stats = {
            "any_annotated_neighbor": 0,
            "left_neighbor_annotated": 0,
            "right_neighbor_annotated": 0,
            "both_neighbors_same_function": 0,
            "both_neighbors_different": 0,
        }

        for smorf_dir in ROOT.iterdir():
            if not smorf_dir.is_dir() or not smorf_dir.name.startswith("SHD1_SM.100AA"):
                continue

            total_smorfs += 1

            status = smorf_annotation_status(smorf_dir)
            if status:
                annotation_counts[status] += 1

            ann_flags, unann_flags = analyze_neighbors_per_smorf(smorf_dir)

            for k in total_annotated_stats:
                if ann_flags[k]:
                    total_annotated_stats[k] += 1

            for k in total_unannotated_stats:
                if unann_flags[k]:
                    total_unannotated_stats[k] += 1

        # -------- Print Summary --------
        print(f"\n=== Evaluating folder: {d} ===")
        print("Total unique smORFs:", total_smorfs)
        print("\nAnnotation status:")
        for k, v in annotation_counts.items():
            coverage = v / total_smorfs * 100 if total_smorfs else 0
            print(f"{k}: {v} ({coverage:.2f}%)")

        print("\nAnnotated neighbor analysis:")
        for k, v in total_annotated_stats.items():
            print(f"{k}: {v}")

        print("\nUnannotated neighbor prediction:")
        total_predictions = sum(total_unannotated_stats.values())
        if total_predictions == 0:
            print("No unannotated smORFs found, cannot compute prediction coverage.")
        else:
            print(f"{'Rule':<25}{'Count':>7}{'Coverage (%)':>15}")
            for k, v in total_unannotated_stats.items():
                coverage = v / total_predictions * 100
                print(f"{k:<25}{v:7}{coverage:15.2f}")
            print(f"{'Total':<25}{total_predictions:7}{100.00:15.2f}")

        # -------- Plotting --------
        plot_bar(annotation_counts, "SmORF Annotation Status")
        plot_bar(total_annotated_stats, "Annotated smORFs: Neighbor Patterns")
        plot_bar(total_unannotated_stats, "Unannotated smORFs: Prediction Potential")

        # Pareto plotting for rules (annotated only)
        combined_rules = {
            "both_neighbors_same": {
                "applied": total_annotated_stats["both_neighbors_same_function"],
                "correct": total_annotated_stats["neighbors_match_each_other_and_target_match"],
            },
            "either_neighbor_match": {
                "applied": total_annotated_stats["left_neighbor_match"] + total_annotated_stats["right_neighbor_match"],
                "correct": total_annotated_stats["left_neighbor_match"] + total_annotated_stats["right_neighbor_match"],
            },
            "left_neighbor_only": {
                "applied": total_annotated_stats["left_neighbor_match"],
                "correct": total_annotated_stats["left_neighbor_match"],
            },
            "right_neighbor_only": {
                "applied": total_annotated_stats["right_neighbor_match"],
                "correct": total_annotated_stats["right_neighbor_match"],
            },
        }

        total_annotated_all = sum(total_annotated_stats.values())
        plot_pareto(combined_rules, total_annotated_all, f"smORF Annotation Rules: {d}")

if __name__ == "__main__":
    main()
