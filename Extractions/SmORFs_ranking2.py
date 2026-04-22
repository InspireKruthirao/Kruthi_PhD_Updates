#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv

# ---------------- CONFIG ----------------
BASE = Path("/work/microbiome/users/kruthi")

RANGE_DIRS = [
    "SmORF_neighbourhoods_1_5_no_min_overlap",
    "SmORF_neighbourhoods_6_10_no_min_overlap",
    "SmORF_neighbourhoods_11_15_no_min_overlap",
    "SmORF_neighbourhoods_16_20_no_min_overlap",
    "SmORF_neighbourhoods_21_25_no_min_overlap",
    "SmORF_neighbourhoods_26_30_no_min_overlap",
]

OUTPUT_DIR = BASE / "22_April"
OUTPUT_DIR.mkdir(exist_ok=True)

# ---------------- HELPERS ----------------
def parse_attrs(s):
    d = {}
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def cog_label(name):
    if not name or not name.startswith("COG") or "-" not in name:
        return ""
    return name.split("-")[-1][0]

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
def process(input_dir):

    results = []
    match_counts = Counter()
    annotated_total = 0

    for cluster in input_dir.iterdir():
        if not cluster.name.startswith("SHD1_SM"):
            continue

        votes = []
        target_labels = []

        for gff in cluster.rglob("*.gff"):
            feats = read_gff(gff)

            for i, f in enumerate(feats):
                if not f["is_target"]:
                    continue

                tc = f["cog"]
                if tc:
                    target_labels.append(tc)

                up = feats[i-1]["cog"] if i > 0 else ""
                down = feats[i+1]["cog"] if i < len(feats)-1 else ""

                if up:
                    votes.append(up)
                if down:
                    votes.append(down)

        # -------- TRUE LABEL --------
        if target_labels:
            true_label = Counter(target_labels).most_common(1)[0][0]
        else:
            true_label = "NA"

        total_votes = len(votes)

        if total_votes == 0:
            predicted = "NA"
            distribution = "-"
            confidence = 0
            decision = "no_signal"

        else:
            counts = Counter(votes)
            sorted_labels = counts.most_common()

            top_label, top_count = sorted_labels[0]
            confidence = round(top_count / total_votes * 100, 2)

            if len(sorted_labels) > 1 and sorted_labels[0][1] == sorted_labels[1][1]:
                predicted = f"{sorted_labels[0][0]}/{sorted_labels[1][0]}"
                decision = "tie"
            else:
                predicted = top_label
                decision = "majority"

            top3 = sorted_labels[:3]
            distribution = ", ".join(
                [f"{k}:{round(v/total_votes*100,2)}%" for k, v in top3]
            )

        # -------- MATCH --------
        if true_label == "NA":
            match = "NA"
        else:
            annotated_total += 1

            if predicted == "NA":
                match = "NA"
            elif "/" in predicted:
                match = "PARTIAL" if true_label in predicted.split("/") else "NOT_MATCH"
            elif predicted == true_label:
                match = "MATCH"
            else:
                match = "NOT_MATCH"

            match_counts[match] += 1

        results.append({
            "cluster": cluster.name,
            "true_label": true_label,
            "predicted": predicted,
            "distribution": distribution,
            "confidence_%": confidence,
            "decision": decision,
            "match": match
        })

    return results, match_counts, annotated_total

# ---------------- SAVE ----------------
def save_results(results, folder):

    out_file = OUTPUT_DIR / f"{folder}_predictions.tsv"

    with open(out_file, "w") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Saved → {out_file}")

def save_summary(match_counts, annotated_total, folder):

    out_file = OUTPUT_DIR / f"{folder}_accuracy.txt"

    with open(out_file, "w") as f:
        f.write(f"Accuracy Summary ({folder})\n")
        f.write("=============================\n\n")

        for k, v in match_counts.items():
            pct = round(v/annotated_total*100, 2)
            f.write(f"{k}: {v} ({pct}%)\n")

        strict_acc = match_counts["MATCH"] / annotated_total * 100
        f.write(f"\nStrict accuracy: {round(strict_acc,2)}%\n")

        if "PARTIAL" in match_counts:
            partial_acc = (match_counts["MATCH"] + match_counts["PARTIAL"]) / annotated_total * 100
            f.write(f"Accuracy (with partial): {round(partial_acc,2)}%\n")

    print(f"Saved → {out_file}")

# ---------------- MAIN ----------------
def main():

    for folder in RANGE_DIRS:

        print(f"\nProcessing {folder}")

        input_dir = BASE / folder

        results, match_counts, annotated_total = process(input_dir)

        print(f"Annotated smORFs: {annotated_total}")

        for k, v in match_counts.items():
            print(f"{k}: {v} ({v/annotated_total*100:.1f}%)")

        save_results(results, folder)
        save_summary(match_counts, annotated_total, folder)

    print(f"\nAll outputs saved in → {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
