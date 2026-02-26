#!/usr/bin/env python3
from pathlib import Path

BASE = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_26_30_no_min_overlap"]

FEATURE_TYPE = "CDS"
UNKNOWN = "Unknown"

# -------------------- Utilities --------------------
def attrs(s: str) -> dict:
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
    """Return set of letters from COG annotation"""
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return set()
    i = str(name).rfind("-")
    if i == -1:
        return set()
    return set(ch for ch in str(name)[i + 1:] if ch.isalpha())

def read_gff_with_order(gff_file):
    """Return list of (is_target, cog_set) in order"""
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

                # Annotated target
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
                # Unannotated target
                else:
                    annotated_neighbors = [(c) for c in [left_cog, right_cog] if c]
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

# -------------------- Main --------------------
def main():
    for d in RANGE_DIRS:
        ROOT = BASE / d
        print("\n" + "="*80)
        print(d)
        print("="*80)

        total_smorfs = 0
        annotation_counts = {"all_annotated":0, "all_unannotated":0, "mixed":0}
        total_annotated_stats = {k:0 for k in ["both_neighbors_same_function","left_neighbor_match","right_neighbor_match","neighbors_match_each_other_target_diff","neighbors_match_each_other_and_target_match"]}
        total_unannotated_stats = {k:0 for k in ["any_annotated_neighbor","left_neighbor_annotated","right_neighbor_annotated","both_neighbors_same_function","both_neighbors_different"]}

        for smorf_dir in ROOT.iterdir():
            if not smorf_dir.is_dir() or not smorf_dir.name.startswith("SHD1_SM.100AA"):
                continue

            total_smorfs += 1

            # --- annotation status ---
            status = smorf_annotation_status(smorf_dir)
            if status:
                annotation_counts[status] += 1

            # --- neighbor patterns ---
            ann_flags, unann_flags = analyze_neighbors_per_smorf(smorf_dir)
            for k in total_annotated_stats:
                if ann_flags[k]:
                    total_annotated_stats[k] += 1
            for k in total_unannotated_stats:
                if unann_flags[k]:
                    total_unannotated_stats[k] += 1

        # --- Print results ---
        print(f"Total unique smORFs: {total_smorfs}\n")

        print("--- SmORF annotation status ---")
        for k,v in annotation_counts.items():
            print(f"{k}: {v} ({v/total_smorfs:.4f})")

        print("\n--- Annotated neighbor analysis (per smORF) ---")
        for k,v in total_annotated_stats.items():
            print(f"{k}: {v}")

        print("\n--- Unannotated neighbor-based prediction (per smORF) ---")
        for k,v in total_unannotated_stats.items():
            print(f"{k}: {v}")

if __name__ == "__main__":
    main()

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


def main():

    for d in RANGE_DIRS:
        ROOT = BASE / d

        print("\n" + "="*80)
        print(d)
        print("="*80)

        total_smorf = 0
        all_annotated = 0
        all_unannotated = 0
        mixed = 0

        for smorf_dir in ROOT.iterdir():
            if not smorf_dir.is_dir() or not smorf_dir.name.startswith("SHD1_SM.100AA"):
                continue

            total_smorf += 1

            status = smorf_annotation_status(smorf_dir)

            if status == "all_annotated":
                all_annotated += 1
            elif status == "all_unannotated":
                all_unannotated += 1
            elif status == "mixed":
                mixed += 1

        print(f"Total unique smORFs: {total_smorf}")
        print(f"All annotated: {all_annotated}")
        print(f"All unannotated: {all_unannotated}")
        print(f"Mixed annotation: {mixed}")

        if total_smorf > 0:
            print(f"\nAll annotated %: {all_annotated / total_smorf:.4f}")
            print(f"All unannotated %: {all_unannotated / total_smorf:.4f}")
            print(f"Mixed %: {mixed / total_smorf:.4f}")

        print("="*80)


if __name__ == "__main__":
    main()
