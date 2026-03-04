#!/usr/bin/env python3
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
      {rule: {"eligible": bool, "predicted_correct": bool}}

    The smORF is treated as ANNOTATED if majority of occurrences have an
    annotated target.  Representative COGs are majority-vote winners.

    "predicted_correct" is evaluated against the true majority target COG.
    Only meaningful when eligible=True.
    """
    target_cog    = majority(evidence["target_cog_votes"])
    left_cog      = majority(evidence["left_cog_votes"])
    right_cog     = majority(evidence["right_cog_votes"])
    is_ann_target = evidence["n_annotated"] > evidence["n_unannotated"]

    results = {}

    # R1: Both neighbours have the same annotated function
    r1_eligible   = is_annotated(left_cog) and is_annotated(right_cog) and left_cog == right_cog
    r1_predicted  = r1_eligible and is_ann_target and target_cog == left_cog
    results["R1_both_same"] = {"eligible": r1_eligible, "predicted_correct": r1_predicted}

    # R2: Left (upstream) neighbour annotated -> predict target = left_cog
    r2_eligible   = is_annotated(left_cog)
    r2_predicted  = r2_eligible and is_ann_target and target_cog == left_cog
    results["R2_left_match"] = {"eligible": r2_eligible, "predicted_correct": r2_predicted}

    # R3: Right (downstream) neighbour annotated -> predict target = right_cog
    r3_eligible   = is_annotated(right_cog)
    r3_predicted  = r3_eligible and is_ann_target and target_cog == right_cog
    results["R3_right_match"] = {"eligible": r3_eligible, "predicted_correct": r3_predicted}

    # R4: Either neighbour annotated -> predict with whichever matches
    r4_eligible  = is_annotated(left_cog) or is_annotated(right_cog)
    if r4_eligible and is_ann_target:
        r4_predicted = (is_annotated(left_cog)  and target_cog == left_cog) or \
                       (is_annotated(right_cog) and target_cog == right_cog)
    else:
        r4_predicted = False
    results["R4_either_match"] = {"eligible": r4_eligible, "predicted_correct": r4_predicted}

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
        global_rule_counters = {r: {"eligible": 0, "predicted_correct": 0} for r in RULE_KEYS}
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
                    if rule_results[r]["eligible"]:
                        global_rule_counters[r]["eligible"]           += 1
                        global_rule_counters[r]["predicted_correct"] += int(rule_results[r]["predicted_correct"])


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
        print(f"  NOTE: 'Eligible' = smORFs where the rule had enough info to make a prediction;")
        print(f"          'Predicted correct' = predictions that matched the true COG.")
        print(f"  {'Rule':<38} {'Eligible':>9} {'Predicted':>10} {'Cov%':>7} {'Acc%':>7}")
        print(f"  {'-'*68}")
        for r in RULE_KEYS:
            eligible          = global_rule_counters[r]["eligible"]
            predicted_correct = global_rule_counters[r]["predicted_correct"]
            cov = eligible          / global_total_ann * 100 if global_total_ann  else 0
            acc = predicted_correct / eligible          * 100 if eligible           else 0
            flag = "  <<< hits 80% target" if acc >= ACCURACY_TARGET else ""
            print(f"  {RULE_LABELS[r].replace(chr(10),' '):<38} "
                  f"{eligible:>9} {predicted_correct:>10} {cov:>7.1f} {acc:>7.1f}{flag}")


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
