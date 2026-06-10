#!/usr/bin/env python3
"""
Batch runner — runs the one-hot weighted vote prediction over ALL
*_no_min_overlap directories found under a given root.

Run as:
    python3 run_all_no_min_overlap.py

Edit ROOT_DIR below if needed.
"""

import re
import sys
from pathlib import Path
from collections import defaultdict, Counter

# ── EDIT THIS LINE ─────────────────────────────────────────────────────────────
ROOT_DIR = Path("/work/microbiome/users/kruthi")
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
            attrs = parts[8]
            name_m = re.search(r'Name=([^;]+)', attrs)
            is_tgt = 'target=1' in attrs
            name = name_m.group(1) if name_m else "Unknown"
            strand = parts[6]
            start = int(parts[3])
            cog_cats = name.split("-")[1] if (name != "Unknown" and "-" in name) else None
            genes.append({
                "name": name, "cog_cats": cog_cats,
                "strand": strand, "start": start, "is_target": is_tgt,
            })
    genes.sort(key=lambda g: g["start"])
    return genes


def normalize(genes):
    target_strand = next((g["strand"] for g in genes if g["is_target"]), "+")
    if target_strand == "-":
        return list(reversed(genes)), True
    return genes, False


def extract_true_cogs(name):
    if not name or "-" not in name:
        return None
    cog_part = name.split("-", 1)[1]
    valid = {c for c in cog_part if c in COG_DESC and c not in ("S", "R")}
    return valid if valid else None


def vote_one_occurrence(norm_genes, norm_target_idx):
    scores = defaultdict(float)
    for i, gene in enumerate(norm_genes):
        if gene["is_target"]:
            continue
        rel_pos = i - norm_target_idx
        pos_weight = POSITION_WEIGHTS.get(rel_pos, 0.0)
        if pos_weight == 0:
            continue
        cog_cats = gene["cog_cats"]
        if not cog_cats or set(cog_cats) <= {"S", "R"}:
            continue
        valid = [c for c in cog_cats if c in COG_DESC and c not in ("S", "R")]
        if not valid:
            continue
        weight_each = pos_weight / len(valid)
        for cat in valid:
            scores[cat] += weight_each
    return scores


def predict_one_smorf(smorf_dir):
    smorf_id = smorf_dir.name
    gff_files = sorted(smorf_dir.glob("*/*.gff"))
    if not gff_files:
        return None

    total_scores = defaultdict(float)
    left_counts = Counter()
    right_counts = Counter()
    target_name = "Unknown"
    true_cogs = None
    n_parsed = 0

    for gff_file in gff_files:
        genes = parse_gff(gff_file)
        target_idx = next((i for i, g in enumerate(genes) if g["is_target"]), None)
        if target_idx is None:
            continue
        norm, _ = normalize(genes)
        norm_ti = next(i for i, g in enumerate(norm) if g["is_target"])
        if target_name == "Unknown":
            target_name = norm[norm_ti]["name"]
            true_cogs = extract_true_cogs(target_name)
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
    total_vote = sum(total_scores.values())
    if total_vote <= 0:
        return None

    score_pct = {cat: round(v / total_vote * 100, 1) for cat, v in total_scores.items()}
    ranked = sorted(score_pct.items(), key=lambda x: -x[1])
    top_cat, top_score = ranked[0]
    top3_str = " | ".join(f"{c}:{p}%" for c, p in ranked[:3])

    cumsum = 0.0
    pred_set = []
    for cat, pct in ranked:
        pred_set.append(cat)
        cumsum += pct
        if cumsum >= 90:
            break

    left_pct = (left_counts.most_common(1)[0][1] / n_parsed * 100) if left_counts else 0.0
    right_pct = (right_counts.most_common(1)[0][1] / n_parsed * 100) if right_counts else 0.0
    conservation = (left_pct + right_pct) / 2
    top_left = left_counts.most_common(1)[0][0] if left_counts else "None"
    top_right = right_counts.most_common(1)[0][0] if right_counts else "None"

    if top_score >= 50 and conservation >= 90 and n_parsed >= 15:
        tier = "VERY HIGH"
    elif top_score >= 40 and conservation >= 75 and n_parsed >= 8:
        tier = "HIGH"
    elif top_score >= 30 and conservation >= 50 and n_parsed >= 4:
        tier = "MEDIUM"
    else:
        tier = "LOW"

    correct = (true_cogs is not None and top_cat in true_cogs)

    return {
        "smorf_id": smorf_id,
        "current_annotation": target_name,
        "true_cogs": ",".join(sorted(true_cogs)) if true_cogs else None,
        "correct": correct,
        "n_occurrences": n_parsed,
        "conservation_pct": round(conservation, 1),
        "confidence_tier": tier,
        "predicted_cog_cat": top_cat,
        "vote_score_pct": top_score,
        "predicted_function": COG_DESC.get(top_cat, "?"),
        "top3": top3_str,
        "pred_set_90pct": ",".join(pred_set),
        "top_left_neighbor": top_left,
        "top_left_pct": round(left_pct, 1),
        "top_right_neighbor": top_right,
        "top_right_pct": round(right_pct, 1),
        "all_scores": " | ".join(f"{c}:{p}%" for c, p in ranked),
    }


def print_validation_summary(results, label=""):
    validated = [r for r in results if r["true_cogs"] is not None]
    header = f"\nVALIDATION{' — ' + label if label else ''}"
    print(header)
    print("-" * len(header.strip()))
    if not validated:
        print("No annotated smORFs with a usable true COG were found.")
        return
    correct = sum(1 for r in validated if r["correct"])
    acc = correct / len(validated)
    print(f"Annotated smORFs : {len(validated)}")
    print(f"Correct          : {correct}")
    print(f"Accuracy         : {acc:.3%}")

    print("\nACCURACY BY TIER")
    print("----------------")
    for tier in ["VERY HIGH", "HIGH", "MEDIUM", "LOW"]:
        subset = [r for r in validated if r["confidence_tier"] == tier]
        if not subset:
            continue
        tier_acc = sum(1 for r in subset if r["correct"]) / len(subset)
        print(f"{tier:<10}  N={len(subset):5d}  Acc={tier_acc:.3%}")


COLS = [
    "rank", "smorf_id", "true_cogs", "correct", "current_annotation",
    "n_occurrences", "conservation_pct", "confidence_tier",
    "predicted_cog_cat", "vote_score_pct", "predicted_function",
    "top3_candidates", "prediction_set_90pct",
    "top_left_neighbor", "top_left_pct",
    "top_right_neighbor", "top_right_pct", "all_scores",
]

TIER_ORDER = {"VERY HIGH": 0, "HIGH": 1, "MEDIUM": 2, "LOW": 3}


def write_tsv(results, out_file):
    results_sorted = sorted(results, key=lambda r: (
        TIER_ORDER.get(r["confidence_tier"], 4),
        -r["vote_score_pct"],
        -r["n_occurrences"],
    ))
    with out_file.open("w") as f:
        f.write("\t".join(COLS) + "\n")
        for rank, r in enumerate(results_sorted, 1):
            row = [
                rank, r["smorf_id"], r["true_cogs"], r["correct"],
                r["current_annotation"], r["n_occurrences"],
                r["conservation_pct"], r["confidence_tier"],
                r["predicted_cog_cat"], r["vote_score_pct"],
                r["predicted_function"], r["top3"], r["pred_set_90pct"],
                r["top_left_neighbor"], r["top_left_pct"],
                r["top_right_neighbor"], r["top_right_pct"], r["all_scores"],
            ]
            f.write("\t".join(str(x) for x in row) + "\n")
    return results_sorted


def run_one_bin_dir(bin_dir):
    smorf_dirs = sorted([
        d for d in bin_dir.iterdir()
        if d.is_dir() and d.name.startswith("SHD1_SM") and any(d.glob("*/*.gff"))
    ])
    if not smorf_dirs:
        print(f"  [skip] no SHD1_SM* subdirs with GFFs found")
        return []

    print(f"  smORFs found : {len(smorf_dirs)}")
    results, failed, done = [], [], 0
    for smorf_dir in smorf_dirs:
        result = predict_one_smorf(smorf_dir)
        if result:
            results.append(result)
        else:
            failed.append(smorf_dir.name)
        done += 1
        if done % 100 == 0 or done == len(smorf_dirs):
            print(f"    {done:>5} / {len(smorf_dirs)}  predicted: {len(results)}  failed: {len(failed)}")

    out_file = bin_dir / "smorf_onehot_predictions_validation.tsv"
    write_tsv(results, out_file)
    print(f"  Saved → {out_file}")
    print(f"  Total predicted : {len(results)}  |  Failed/no data : {len(failed)}")
    print_validation_summary(results, label=bin_dir.name)
    return results


def main():
    if not ROOT_DIR.exists():
        print(f"ERROR: {ROOT_DIR} does not exist")
        sys.exit(1)

    overlap_dirs = sorted([
        d for d in ROOT_DIR.iterdir()
        if d.is_dir() and "no_min_overlap" in d.name
    ])

    if not overlap_dirs:
        print(f"No *no_min_overlap* directories found under {ROOT_DIR}")
        sys.exit(1)

    print(f"Root      : {ROOT_DIR}")
    print(f"Matched directories ({len(overlap_dirs)}):")
    for d in overlap_dirs:
        print(f"  {d.name}")
    print()

    all_results = []
    for bin_dir in overlap_dirs:
        print(f"\n{'='*60}")
        print(f"Processing: {bin_dir.name}")
        print(f"{'='*60}")
        results = run_one_bin_dir(bin_dir)
        all_results.extend(results)

    # ── Combined summary across all overlap dirs ──────────────────────────────
    if all_results:
        print(f"\n{'='*60}")
        print(f"COMBINED SUMMARY — all no_min_overlap dirs")
        print(f"{'='*60}")
        print(f"Total smORFs predicted : {len(all_results)}")
        print_validation_summary(all_results, label="ALL no_min_overlap dirs")

        combined_out = ROOT_DIR / "smorf_onehot_predictions_ALL_no_min_overlap.tsv"
        write_tsv(all_results, combined_out)
        print(f"\nCombined TSV saved → {combined_out}")


if __name__ == "__main__":
    main()
