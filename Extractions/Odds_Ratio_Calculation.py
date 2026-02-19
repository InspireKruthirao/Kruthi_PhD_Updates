#!/usr/bin/env python3
from pathlib import Path
import math

BASE = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_26_30_no_min_overlap"]

FEATURE_TYPE = "CDS"
UNKNOWN = "Unknown"


def attrs(s: str) -> dict:
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def cog_set(name: str):
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return set()
    i = str(name).rfind("-")
    if i == -1:
        return set()
    return set(ch for ch in str(name)[i + 1:] if ch.isalpha())


def iter_smorf_gffs(root: Path):
    for smorf in root.iterdir():
        if not smorf.is_dir():
            continue
        sid = smorf.name
        if sid.endswith(".tsv") or sid.endswith(".csv"):
            continue
        for contig in smorf.iterdir():
            if not contig.is_dir():
                continue
            for gff in contig.glob("*.gff"):
                yield sid, gff


def read_target_neighbors(gff: Path):
    feats = []
    with gff.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue
            a = attrs(parts[8])
            name = a.get("Name", UNKNOWN)
            target = a.get("target", "0")
            feats.append((int(parts[3]), int(parts[4]), name, target))

    feats.sort(key=lambda x: (x[0], x[1]))
    t = next((i for i, x in enumerate(feats) if x[3] == "1"), None)
    if t is None:
        return None

    t_set = cog_set(feats[t][2])
    up_set = cog_set(feats[t - 1][2]) if t > 0 else set()
    dn_set = cog_set(feats[t + 1][2]) if t + 1 < len(feats) else set()
    return t_set, up_set, dn_set


def odds_ratio_and_ci(a, b, c, d):
    # Haldane–Anscombe correction if any cell is zero
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    OR = (a * d) / (b * c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    lo = math.exp(math.log(OR) - 1.96 * se)
    hi = math.exp(math.log(OR) + 1.96 * se)
    return OR, lo, hi


def smorf_level_2x2(ROOT: Path):
    """
    One row per smORF ID (annotated only).
    Exposure:
      - neighbors ever same: ever had both neighbors annotated and up_set == dn_set
      - neighbors never same: otherwise, but must have both neighbors annotated at least once
    Outcome:
      - ever matches (same context): t_set overlaps neighbor set in SAME-neighbor occurrences
      - ever matches (diff context): t_set overlaps either neighbor in DIFF-neighbor occurrences
    """
    st = {}  # sid -> flags

    for sid, gff in iter_smorf_gffs(ROOT):
        rec = read_target_neighbors(gff)
        if rec is None:
            continue
        t_set, up_set, dn_set = rec

        # only evaluate where smORF itself is annotated
        if not t_set:
            continue

        if sid not in st:
            st[sid] = {
                "has_both_neighbors": False,
                "ever_same": False,
                "match_in_same": False,
                "match_in_diff": False,
            }

        if not up_set or not dn_set:
            continue

        st[sid]["has_both_neighbors"] = True

        if up_set == dn_set:
            st[sid]["ever_same"] = True
            if t_set & up_set:
                st[sid]["match_in_same"] = True
        else:
            if (t_set & up_set) or (t_set & dn_set):
                st[sid]["match_in_diff"] = True

    a = b = c = d = 0
    for sid, f in st.items():
        if not f["has_both_neighbors"]:
            continue

        if f["ever_same"]:
            if f["match_in_same"]:
                a += 1
            else:
                b += 1
        else:
            if f["match_in_diff"]:
                c += 1
            else:
                d += 1

    return a, b, c, d


def main():
    for dname in RANGE_DIRS:
        ROOT = BASE / dname
        print("\n" + "=" * 80)
        print(dname)
        print("=" * 80)

        if not ROOT.exists():
            print(f"[SKIP] Not found: {ROOT}")
            continue

        a, b, c, d = smorf_level_2x2(ROOT)

        print("\n2x2 (Odds Ratio)")
        print(f"a (ever same + match)    : {a}")
        print(f"b (ever same + no match) : {b}")
        print(f"c (never same + match)   : {c}")
        print(f"d (never same + no match): {d}")

        OR, lo, hi = odds_ratio_and_ci(a, b, c, d)
        print(f"\nOdds Ratio = {OR:.3f}")
        print(f"95% CI     = [{lo:.3f}, {hi:.3f}]")
        print(f"Total used = {a+b+c+d}")


if __name__ == "__main__":
    main()
