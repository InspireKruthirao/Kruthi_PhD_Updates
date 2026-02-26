#!/usr/bin/env python3

from pathlib import Path

BASE = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = [
    "SmORF_neighbourhoods_26_30_no_min_overlap",
]

FEATURE_TYPE = "CDS"
UNKNOWN = "Unknown"


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


def smorf_annotation_status(folder: Path):
    """
    For one smORF folder, determine:
    - all_annotated
    - all_unannotated
    - mixed
    """

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
