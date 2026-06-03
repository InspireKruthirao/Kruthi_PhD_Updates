#!/usr/bin/env python3
"""
One-Hot Weighted Vote — Batch Prediction
Run as: python3 <this_script>.py
No arguments needed. Just edit the two paths at the top.
"""

import re
import sys
from pathlib import Path
from collections import defaultdict, Counter

# ── EDIT THESE TWO LINES ───────────────────────────────────────────────────────
BIN_DIR  = Path("/work/microbiome/users/kruthi/SmORF_neighbourhoods_21_25_no_min_overlap")
OUT_FILE = BIN_DIR / "smorf_onehot_predictions.tsv"
# ──────────────────────────────────────────────────────────────────────────────

COG_DESC = {
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control and division',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation and ribosome',
    'K': 'Transcription',
    'L': 'DNA replication and repair',
    'M': 'Cell wall and membrane',
    'N': 'Cell motility',
    'O': 'Protein turnover and chaperones',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolite biosynthesis',
    'R': 'General function prediction',
    'S': 'Function unknown',
    'T': 'Signal transduction',
    'U': 'Intracellular trafficking',
    'V': 'Defense mechanisms',
}

POSITION_WEIGHTS = {
    -5: 0.2, -4: 0.4, -3: 0.6, -2: 0.8, -1: 1.0,
     1: 1.0,  2: 0.8,  3: 0.6,  4: 0.4,  5: 0.2,
}


def parse_gff(gff_file):
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs    = parts[8]
            name_m   = re.search(r'Name=([^;]+)', attrs)
            is_tgt   = 'target=1' in attrs
            name     = name_m.group(1) if name_m else "Unknown"
            strand   = parts[6]
            start    = int(parts[3])
            cog_cats = name.split("-")[1] if (
                name != "Unknown" and "-" in name
            ) else None
            genes.append({
                "name": name, "cog_cats": cog_cats,
                "strand": strand, "start": start,
                "is_target": is_tgt,
            })
    genes.sort(key=lambda g: g["start"])
    return genes


def normalize(genes):
    target_strand = next(
        (g["strand"] for g in genes if g["is_target"]), "+"
    )
    if target_strand == "-":
        return list(reversed(genes)), True
    return genes, False


def vote_one_occurrence(norm_genes, norm_target_idx):
    scores = defaultdict(float)
    for i, gene in enumerate(norm_genes):
        if gene["is_target"]:
            continue
        rel_pos    = i - norm_target_idx
        pos_weight = POSITION_WEIGHTS.get(rel_pos, 0.0)
        if pos_weight == 0:
            continue
        cog_cats = gene["cog_cats"]
        if not cog_cats or set(cog_cats) <= {'S', 'R'}:
            continue
        valid = [c for c in cog_cats
                 if c in COG_DESC and c not in ('S', 'R')]
        if not valid:
            continue
        weight_each = pos_weight / len(valid)
        for cat in valid:
            scores[cat] += weight_each
    return scores


def predict_one_smorf(smorf_dir):
    smorf_id  = smorf_dir.name
    gff_files = sorted(smorf_dir.glob("*/*.gff"))
    if not gff_files:
        return None

    total_scores = defaultdict(float)
    left_counts  = Counter()
    right_counts = Counter()
    target_cog   = "Unknown"
    n_parsed     = 0

    for gff_file in gff_files:
        genes = parse_gff(gff_file)
        target_idx = next(
            (i for i, g in enumerate(genes) if g["is_target"]), None
        )
        if target_idx is None:
            continue
        norm, _ = normalize(genes)
        norm_ti  = next(i for i, g in enumerate(norm) if g["is_target"])

        if target_cog == "Unknown":
            target_cog = norm[norm_ti]["name"]

        if norm_ti > 0:
            left_counts[norm[norm_ti - 1]["name"]] += 1
        if norm_ti < len(norm) - 1:
            right_counts[norm[norm_ti + 1]["name"]] += 1

        scores = vote_one_occurrence(norm, norm_ti)
        for cat, score in scores.items():
            total_scores[cat] += score

        n_parsed += 1

    if n_parsed == 0 or not total_scores:
        return None

    n          = n_parsed
    total_vote = sum(total_scores.values())
    score_pct  = {cat: round(v / total_vote * 100, 1)
                  for cat, v in total_scores.items()}
    ranked     = sorted(score_pct.items(), key=lambda x: -x[1])
    top_cat    = ranked[0][0]
    top_score  = ranked[0][1]
    top3_str   = " | ".join(f"{c}:{p}%" for c, p in ranked[:3])

    cumsum   = 0.0
    pred_set = []
    for cat, pct in ranked:
        pred_set.append(cat)
        cumsum += pct
        if cumsum >= 90:
            break

    left_pct     = (left_counts.most_common(1)[0][1]  / n * 100
                    if left_counts  else 0)
    right_pct    = (right_counts.most_common(1)[0][1] / n * 100
                    if right_counts else 0)
    conservation = (left_pct + right_pct) / 2
    top_left     = left_counts.most_common(1)[0][0]  if left_counts  else "None"
    top_right    = right_counts.most_common(1)[0][0] if right_counts else "None"

    if top_score >= 50 and conservation >= 90 and n >= 15:
        tier = "VERY HIGH"
    elif top_score >= 40 and conservation >= 75 and n >= 8:
        tier = "HIGH"
    elif top_score >= 30 and conservation >= 50 and n >= 4:
        tier = "MEDIUM"
    else:
        tier = "LOW"

    return {
        "smorf_id":           smorf_id,
        "current_annotation": target_cog,
        "n_occurrences":      n,
        "conservation_pct":   round(conservation, 1),
        "confidence_tier":    tier,
        "predicted_cog_cat":  top_cat,
        "vote_score_pct":     top_score,
        "predicted_function": COG_DESC.get(top_cat, "?"),
        "top3":               top3_str,
        "pred_set_90pct":     ",".join(pred_set),
        "top_left_neighbor":  top_left,
        "top_left_pct":       round(left_pct,  1),
        "top_right_neighbor": top_right,
        "top_right_pct":      round(right_pct, 1),
        "all_scores":         " | ".join(f"{c}:{p}%" for c, p in ranked),
    }


def main():
    if not BIN_DIR.exists():
        print(f"ERROR: {BIN_DIR} does not exist")
        sys.exit(1)

    smorf_dirs = sorted([
        d for d in BIN_DIR.iterdir()
        if d.is_dir()
        and d.name.startswith("SHD1_SM")
        and any(d.glob("*/*.gff"))
    ])

    print(f"Directory : {BIN_DIR}")
    print(f"smORFs    : {len(smorf_dirs)}")
    print(f"Output    : {OUT_FILE}")
    print()

    results = []
    failed  = []
    done    = 0

    for smorf_dir in smorf_dirs:
        result = predict_one_smorf(smorf_dir)
        if result:
            results.append(result)
        else:
            failed.append(smorf_dir.name)
        done += 1
        if done % 100 == 0 or done == len(smorf_dirs):
            print(f"  {done:>5} / {len(smorf_dirs)}  "
                  f"predicted: {len(results)}  failed: {len(failed)}")

    tier_order = {"VERY HIGH": 0, "HIGH": 1, "MEDIUM": 2, "LOW": 3}
    results.sort(key=lambda r: (
        tier_order.get(r["confidence_tier"], 4),
        -r["vote_score_pct"],
        -r["n_occurrences"]
    ))

    cols = [
        "rank", "smorf_id", "current_annotation",
        "n_occurrences", "conservation_pct", "confidence_tier",
        "predicted_cog_cat", "vote_score_pct", "predicted_function",
        "top3_candidates", "prediction_set_90pct",
        "top_left_neighbor", "top_left_pct",
        "top_right_neighbor", "top_right_pct",
        "all_scores",
    ]
    with OUT_FILE.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for rank, r in enumerate(results, 1):
            row = [
                rank, r["smorf_id"], r["current_annotation"],
                r["n_occurrences"], r["conservation_pct"],
                r["confidence_tier"], r["predicted_cog_cat"],
                r["vote_score_pct"], r["predicted_function"],
                r["top3"], r["pred_set_90pct"],
                r["top_left_neighbor"], r["top_left_pct"],
                r["top_right_neighbor"], r["top_right_pct"],
                r["all_scores"],
            ]
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"\nTotal predicted : {len(results)}")
    print(f"Failed/no data  : {len(failed)}")
    print()

    tier_counts = Counter(r["confidence_tier"] for r in results)
    total = len(results)
    print("CONFIDENCE TIER DISTRIBUTION:")
    for tier in ["VERY HIGH", "HIGH", "MEDIUM", "LOW"]:
        n   = tier_counts.get(tier, 0)
        bar = "█" * int(n / total * 40)
        print(f"  {tier:<12}  {n:5d}  ({n/total*100:5.1f}%)  {bar}")

    print()
    print("TOP 10 PREDICTIONS:")
    print(f"  {'smORF ID':<35}  {'Tier':<10}  {'Cat':>4}  {'Score':>7}  Function")
    print(f"  {'-'*35}  {'-'*10}  {'---':>4}  {'-----':>7}  {'--------'}")
    for r in results[:10]:
        print(
            f"  {r['smorf_id']:<35}  "
            f"{r['confidence_tier']:<10}  "
            f"{r['predicted_cog_cat']:>4}  "
            f"{r['vote_score_pct']:>6.1f}%  "
            f"{r['predicted_function'][:35]}"
        )

    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    main()
