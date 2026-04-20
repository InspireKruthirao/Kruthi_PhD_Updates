#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv
import matplotlib.pyplot as plt

# ---------------- CONFIG ----------------
BASE = Path("/work/microbiome/users/kruthi")
INPUT_DIR = BASE / "SmORF_neighbourhoods_26_30_no_min_overlap"

# ---------------- HELPERS ----------------
def parse_attrs(s):
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def cog_label(name):
    if not name or not name.startswith("COG") or "-" not in name:
        return ""
    return name.split("-")[-1]

def read_gff(gff):
    feats = []
    with open(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue

            attrs = parse_attrs(parts[8])

            feats.append({
                "is_target": attrs.get("target") == "1",
                "cog": cog_label(attrs.get("Name", ""))
            })

    return feats

# ---------------- CORE ----------------
def process():

    results = []
    decision_counts = Counter()

    for cluster in INPUT_DIR.iterdir():
        if not cluster.name.startswith("SHD1_SM"):
            continue

        votes = []

        for gff in cluster.rglob("*.gff"):
            feats = read_gff(gff)

            for i, f in enumerate(feats):
                if not f["is_target"]:
                    continue

                up = feats[i-1]["cog"] if i > 0 else ""
                down = feats[i+1]["cog"] if i < len(feats)-1 else ""

                # 🔥 NEW LOGIC: use BOTH neighbours
                if up:
                    votes.append(up)
                if down:
                    votes.append(down)

        total_votes = len(votes)

        if total_votes == 0:
            result = {
                "cluster": cluster.name,
                "predicted": "NA",
                "distribution": "-",
                "confidence_%": 0,
                "decision": "no_signal"
            }
            decision_counts["no_signal"] += 1

        else:
            counts = Counter(votes)

            # percentages
            dist = {k: round(v/total_votes*100, 2) for k,v in counts.items()}

            # sort by frequency
            sorted_labels = counts.most_common()

            top_label, top_count = sorted_labels[0]
            top_pct = dist[top_label]

            # check tie
            if len(sorted_labels) > 1 and sorted_labels[0][1] == sorted_labels[1][1]:
                predicted = f"{sorted_labels[0][0]}/{sorted_labels[1][0]}"
                decision = "tie"
            else:
                predicted = top_label
                decision = "majority"

            # build distribution string
            dist_str = ", ".join([f"{k}:{v}%" for k,v in dist.items()])

            result = {
                "cluster": cluster.name,
                "predicted": predicted,
                "distribution": dist_str,
                "confidence_%": top_pct,
                "decision": decision
            }

            decision_counts[decision] += 1

        results.append(result)

    return results, decision_counts


# ---------------- SAVE ----------------
def save_results(results):

    with open("smorf_predictions.tsv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print("Saved → smorf_predictions.tsv")


def save_summary(decision_counts):

    total = sum(decision_counts.values())

    with open("summary.txt", "w") as f:
        f.write("Summary\n========\n")
        for k, v in decision_counts.items():
            f.write(f"{k}: {v} ({v/total*100:.1f}%)\n")

    print("Saved → summary.txt")


# ---------------- PLOT ----------------
def plot_stats(decision_counts):

    plt.figure(figsize=(6,4))
    plt.bar(decision_counts.keys(), decision_counts.values())
    plt.title("Prediction Types")
    plt.ylabel("Number of smORFs")
    plt.tight_layout()
    plt.show()


# ---------------- MAIN ----------------
def main():

    results, decision_counts = process()

    print("\nTotal smORFs:", len(results))

    print("\nSample predictions:")
    for r in results[:10]:
        print(r)

    print("\nSummary:")
    total = sum(decision_counts.values())
    for k, v in decision_counts.items():
        print(f"{k}: {v} ({v/total*100:.1f}%)")

    save_results(results)
    save_summary(decision_counts)

    plot_stats(decision_counts)


if __name__ == "__main__":
    main()
