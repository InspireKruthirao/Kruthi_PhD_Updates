#!/usr/bin/env python3
from pathlib import Path
import matplotlib.pyplot as plt

BASE = Path("/work/microbiome/users/kruthi")

FEATURE_TYPE = "CDS"
UNKNOWN = "Unknown"

RANGE_DIRS = [
    "SmORF_neighbourhoods_6_10_no_min_overlap",
    "SmORF_neighbourhoods_11_15_no_min_overlap",
    "SmORF_neighbourhoods_16_20_no_min_overlap",
    "SmORF_neighbourhoods_21_25_no_min_overlap",
    "SmORF_neighbourhoods_26_30_no_min_overlap",
]

# simple, readable palette (no fancy config)
COL_ANNOT = "#1b9e77"
COL_UNANNOT = "#d95f02"
COL_ANNOT_BREAK = "#7570b3"
COL_UNANNOT_BREAK = "#66a61e"


def parse_attrs(attr_field: str) -> dict:
    d = {}
    for item in str(attr_field).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def cog_letters(name: str) -> set[str]:
    """
    From Name like 'COG1538-MU' -> {'M','U'}
    From 'COG0254-J' -> {'J'}
    Unknown / empty -> set()
    """
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return set()
    i = str(name).rfind("-")
    if i == -1:
        return set()
    return {ch for ch in str(name)[i + 1 :] if ch.isalpha()}


def iter_gffs(root: Path):
    # root/<SmORF_ID>/<contig_dir>/*.gff
    for smorf_dir in root.iterdir():
        if not smorf_dir.is_dir():
            continue
        sid = smorf_dir.name
        if sid.endswith((".tsv", ".csv")):
            continue
        for contig_dir in smorf_dir.iterdir():
            if contig_dir.is_dir():
                yield from ((sid, g) for g in contig_dir.glob("*.gff"))


def target_and_neighbors(gff: Path):
    """
    Returns (t_set, up_set, dn_set) for immediate neighbours (±1),
    or None if target not found.
    """
    feats = []
    with gff.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue

            a = parse_attrs(parts[8])
            name = a.get("Name", UNKNOWN)
            target = a.get("target", "0")
            feats.append((int(parts[3]), int(parts[4]), name, target))

    feats.sort(key=lambda x: (x[0], x[1]))
    t_idx = next((i for i, x in enumerate(feats) if x[3] == "1"), None)
    if t_idx is None:
        return None

    t_set = cog_letters(feats[t_idx][2])
    up_set = cog_letters(feats[t_idx - 1][2]) if t_idx > 0 else set()
    dn_set = cog_letters(feats[t_idx + 1][2]) if t_idx + 1 < len(feats) else set()
    return t_set, up_set, dn_set


def compute_stats(root: Path) -> dict:
    # total smORFs = top-level dirs under root
    smorf_ids = {p.name for p in root.iterdir() if p.is_dir() and not p.name.endswith((".tsv", ".csv"))}

    annotated = set()
    # annotated patterns
    a_left_match = set()
    a_right_match = set()
    a_neighbors_same = set()
    a_neighbors_same_target_match = set()
    a_neighbors_same_target_diff = set()

    # unannotated patterns
    u_any_annot_neighbor = set()
    u_left_annot = set()
    u_right_annot = set()
    u_both_same = set()
    u_both_diff = set()

    for sid, gff in iter_gffs(root):
        smorf_ids.add(sid)

        rec = target_and_neighbors(gff)
        if rec is None:
            continue
        t_set, up_set, dn_set = rec

        if t_set:
            annotated.add(sid)

            if up_set and (t_set & up_set):
                a_left_match.add(sid)
            if dn_set and (t_set & dn_set):
                a_right_match.add(sid)

            if up_set and dn_set and (up_set == dn_set):
                a_neighbors_same.add(sid)
                if t_set & up_set:
                    a_neighbors_same_target_match.add(sid)
                else:
                    a_neighbors_same_target_diff.add(sid)

        else:
            # unannotated
            if up_set or dn_set:
                u_any_annot_neighbor.add(sid)
            if up_set:
                u_left_annot.add(sid)
            if dn_set:
                u_right_annot.add(sid)

            if up_set and dn_set:
                if up_set == dn_set:
                    u_both_same.add(sid)
                else:
                    u_both_diff.add(sid)

    total = len(smorf_ids)
    annot = len(annotated)
    unannot = total - annot

    return {
        "Total_SmORFs": total,
        "Annotated_SmORFs": annot,
        "Unannotated_SmORFs": unannot,
        # annotated breakdown
        "A_left_match": len(a_left_match),
        "A_right_match": len(a_right_match),
        "A_neighbors_same": len(a_neighbors_same),
        "A_neighbors_same_target_diff": len(a_neighbors_same_target_diff),
        "A_neighbors_same_target_match": len(a_neighbors_same_target_match),
        # unannotated breakdown
        "U_any_annot_neighbor": len(u_any_annot_neighbor),
        "U_left_annot": len(u_left_annot),
        "U_right_annot": len(u_right_annot),
        "U_both_same": len(u_both_same),
        "U_both_diff": len(u_both_diff),
    }


def _beautify(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.set_axisbelow(True)


def plot_for_bin(bin_name: str, stats: dict):
    annot = stats["Annotated_SmORFs"]
    unannot = stats["Unannotated_SmORFs"]

    # 1) annotated vs unannotated
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    bars = ax.bar(["Annotated", "Unannotated"], [annot, unannot], color=[COL_ANNOT, COL_UNANNOT])
    ax.set_title(f"{bin_name}: Annotated vs Unannotated")
    ax.set_ylabel("Count")
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{int(b.get_height())}",
                ha="center", va="bottom", fontsize=10)
    _beautify(ax)
    fig.tight_layout()
    plt.show()

    # 2) annotated breakdown
    cats_a = [
        "Left neighbor match",
        "Right neighbor match",
        "Both neighbors same function",
        "Neighbors same, target different",
        "Neighbors same, target matches",
    ]
    vals_a = [
        stats["A_left_match"],
        stats["A_right_match"],
        stats["A_neighbors_same"],
        stats["A_neighbors_same_target_diff"],
        stats["A_neighbors_same_target_match"],
    ]

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.barh(cats_a, vals_a, color=COL_ANNOT_BREAK)
    ax.set_title(f"{bin_name}: Annotated smORFs — neighbor patterns")
    ax.set_xlabel("Count")
    for y, v in enumerate(vals_a):
        pct = (100 * v / annot) if annot else 0
        ax.text(v, y, f"{v} ({pct:.1f}%)", va="center", fontsize=9)
    _beautify(ax)
    fig.tight_layout()
    plt.show()

    # 3) unannotated breakdown
    cats_u = [
        "Has any annotated neighbor",
        "Left neighbor annotated",
        "Right neighbor annotated",
        "Both neighbors same function",
        "Both neighbors annotated but different",
    ]
    vals_u = [
        stats["U_any_annot_neighbor"],
        stats["U_left_annot"],
        stats["U_right_annot"],
        stats["U_both_same"],
        stats["U_both_diff"],
    ]

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.barh(cats_u, vals_u, color=COL_UNANNOT_BREAK)
    ax.set_title(f"{bin_name}: Unannotated smORFs — neighbor annotation")
    ax.set_xlabel("Count")
    for y, v in enumerate(vals_u):
        pct = (100 * v / unannot) if unannot else 0
        ax.text(v, y, f"{v} ({pct:.1f}%)", va="center", fontsize=9)
    _beautify(ax)
    fig.tight_layout()
    plt.show()


def main():
    for d in RANGE_DIRS:
        root = BASE / d
        print("\n" + "=" * 80)
        print(d)
        print("=" * 80)

        if not root.exists():
            print(f"[SKIP] Not found: {root}")
            continue

        stats = compute_stats(root)

        print(f"Total_SmORFs: {stats['Total_SmORFs']}")
        print(f"Annotated_SmORFs: {stats['Annotated_SmORFs']}")
        print(f"Unannotated_SmORFs: {stats['Unannotated_SmORFs']}\n")

        print("--- Annotated patterns ---")
        print(f"Left neighbor match: {stats['A_left_match']}")
        print(f"Right neighbor match: {stats['A_right_match']}")
        print(f"Both neighbors same function: {stats['A_neighbors_same']}")
        print(f"Neighbors same, target different: {stats['A_neighbors_same_target_diff']}")
        print(f"Neighbors same, target matches: {stats['A_neighbors_same_target_match']}\n")

        print("--- Unannotated patterns ---")
        print(f"Has any annotated neighbor: {stats['U_any_annot_neighbor']}")
        print(f"Left neighbor annotated: {stats['U_left_annot']}")
        print(f"Right neighbor annotated: {stats['U_right_annot']}")
        print(f"Both neighbors same function: {stats['U_both_same']}")
        print(f"Both neighbors annotated but different: {stats['U_both_diff']}")

        plot_for_bin(d, stats)


if __name__ == "__main__":
    main()
