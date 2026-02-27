#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

BASE = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_26_30_no_min_overlap"]

FEATURE_TYPE = "CDS"
UNKNOWN = "Unknown"

# -------------------- Utilities --------------------

def attrs(s: str) -> dict:
    """Parse GFF attributes into a dictionary"""
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def is_annotated(name: str):
    if not name or name == UNKNOWN:
        return False
    return str(name).startswith("COG")

def cog_set(name: str):
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return set()
    i = str(name).rfind("-")
    if i == -1:
        return set()
    return set(ch for ch in str(name)[i + 1:] if ch.isalpha())

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
