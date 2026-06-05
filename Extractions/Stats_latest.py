# %%
from pathlib import Path
from collections import Counter, defaultdict
from IPython.display import display
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# SET ROOT HERE
# ----------------------------
ROOT = Path("/work/microbiome/users/kruthi/SmORF_neighbourhoods_26_30_no_min_overlap")
FEATURE_TYPE = "CDS"
TOP_N = 20

def attrs(s: str) -> dict:
    d = {}
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def cog_set(name: str):
    """
    From Name like 'COG1538-MU' -> {'M','U'}
    From 'COG0254-J' -> {'J'}
    Unknown / empty -> set()
    """
    if not name or name == "Unknown" or not str(name).startswith("COG"):
        return set()
    i = str(name).rfind("-")
    if i == -1:
        return set()
    return set(ch for ch in str(name)[i + 1:] if ch.isalpha())

def iter_smorf_gffs(root: Path):
    """
    Expected structure:
      ROOT/
        SHD1_SM.../         <-- each is a smORF folder in your bin
          D###_contig_###/
            contig_###.gff
    """
    for smorf in root.iterdir():
        if not smorf.is_dir():
            continue
        sid = smorf.name
        if sid.endswith(".tsv") or sid.endswith(".csv"):
            continue

        for contig in smorf.iterdir():
            if contig.is_dir() and contig.name.startswith("D") and "_contig_" in contig.name:
                for gff in contig.glob("*.gff"):
                    yield sid, gff

def read_target_neighbors(gff: Path):
    """
    Returns:
      t_set, up_set, dn_set  (sets of COG letters)
    or None if no target=1 CDS found.
    """
    feats = []
    with gff.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue
            a = attrs(parts[8])
            feats.append((int(parts[3]), int(parts[4]), a.get("Name"), a.get("target")))

    feats.sort(key=lambda x: (x[0], x[1]))
    t = next((i for i, x in enumerate(feats) if x[3] == "1"), None)
    if t is None:
        return None

    t_set = cog_set(feats[t][2])
    up_set = cog_set(feats[t - 1][2]) if t > 0 else set()
    dn_set = cog_set(feats[t + 1][2]) if t + 1 < len(feats) else set()
    return t_set, up_set, dn_set

def pct(a, b):
    return (a / b * 100) if b else 0.0

def dominant_letter_from_counter(ctr: Counter) -> str:
    """Tie-break alphabetically."""
    if not ctr:
        return ""
    mx = max(ctr.values())
    return min([L for L, c in ctr.items() if c == mx])

def pick_single_letter_from_set(s: set) -> str:
    """For neighbour shared set in L==R, pick a single representative (sorted) for pair tables."""
    if not s:
        return ""
    return "".join(sorted(s))  # keep full multi-letter like 'MU' if present


# %%
# -------- MAIN COUNTS --------
smorf_dirs = [p for p in ROOT.iterdir() if p.is_dir()]
total_smorf = len(smorf_dirs)

annot_smorf = set()
target_votes = defaultdict(Counter)

# Occurrence-level tracking
occ_total = 0
occ_annot = 0
occ_unknown = 0
letter_counts_per_occ = Counter()

# Any-neighbour match (occ + smORF)
occ_annot_total = 0
occ_match_any = 0
smorf_occ_annot = Counter()
smorf_match_any = Counter()

# L==R (occurrence-level)
occ_LeqR_annot = 0
occ_LeqR_match = 0
occ_LeqR_diff = 0  # L==R but target different (no overlap)

# L==R (smORF-level)
smorf_LeqR_occ = Counter()
smorf_LeqR_match = Counter()
smorf_LeqR_diff = Counter()   # occurrences where L==R and target differs

# L!=R (smORF-level)
smorf_LneqR_occ = Counter()
smorf_LneqR_left_only = Counter()
smorf_LneqR_right_only = Counter()
smorf_LneqR_both = Counter()
smorf_LneqR_neither = Counter()

# For plotting / table
smorf_stats_rows = []

# ---- NEW: store mismatch function info for L==R and target differs ----
# Counts of neighbour shared letters and target letters in these mismatch cases
LeqR_diff_neighbor_letters = Counter()  # counts of shared neighbour letter-sets (e.g., 'J', 'MU')
LeqR_diff_target_letters = Counter()    # counts of target letters (per occurrence, multi-letter contributes multiple)
LeqR_diff_pair_counts = Counter()       # (neighbor_shared_str, target_dom_letter) -> count

# Also smORF-level mismatch summary (which smORFs show this)
smorf_LeqR_diff_smorfs = set()
smorf_LeqR_diff_pairs = defaultdict(Counter)  # sid -> Counter[(neighbor_shared, target_dom)] += 1

for sid, gff in iter_smorf_gffs(ROOT):
    rec = read_target_neighbors(gff)
    if rec is None:
        continue

    t_set, up_set, dn_set = rec
    occ_total += 1

    # Unknown/unannotated target
    if not t_set:
        occ_unknown += 1
        continue

    # Annotated occurrence
    occ_annot += 1
    annot_smorf.add(sid)

    # target vote letters (dominant letter per smORF)
    for L in t_set:
        letter_counts_per_occ[L] += 1
        target_votes[sid][L] += 1

    # any-neighbour match at occurrence
    occ_annot_total += 1
    smorf_occ_annot[sid] += 1
    has_any_match = bool((t_set & up_set) or (t_set & dn_set))
    if has_any_match:
        occ_match_any += 1
        smorf_match_any[sid] += 1

    # ---- L==R analysis (require both neighbours have letters) ----
    if up_set and dn_set and (up_set == dn_set):
        occ_LeqR_annot += 1
        smorf_LeqR_occ[sid] += 1

        shared = up_set  # == dn_set
        if bool(t_set & shared):
            occ_LeqR_match += 1
            smorf_LeqR_match[sid] += 1
        else:
            # L==R but different
            occ_LeqR_diff += 1
            smorf_LeqR_diff[sid] += 1
            smorf_LeqR_diff_smorfs.add(sid)

            neigh_str = pick_single_letter_from_set(shared)
            # target dominant letter for this smORF occurrence based on its current t_set:
            # (for per-occurrence, pick a single letter representative)
            target_dom = min(sorted(t_set))  # deterministic; could also keep full t_set string
            LeqR_diff_neighbor_letters[neigh_str] += 1
            for L in t_set:
                LeqR_diff_target_letters[L] += 1
            LeqR_diff_pair_counts[(neigh_str, target_dom)] += 1
            smorf_LeqR_diff_pairs[sid][(neigh_str, target_dom)] += 1

    # ---- L!=R analysis (require both neighbours have letters) ----
    if up_set and dn_set and (up_set != dn_set):
        smorf_LneqR_occ[sid] += 1
        match_left = bool(t_set & up_set)
        match_right = bool(t_set & dn_set)

        if match_left and (not match_right):
            smorf_LneqR_left_only[sid] += 1
        elif (not match_left) and match_right:
            smorf_LneqR_right_only[sid] += 1
        elif match_left and match_right:
            smorf_LneqR_both[sid] += 1
        else:
            smorf_LneqR_neither[sid] += 1


# %%
# -------- PER-smORF DOMINANT LETTER --------
dominant_target = {}
letter_counts_per_smorf = Counter()

for sid, ctr in target_votes.items():
    dom = dominant_letter_from_counter(ctr)
    if not dom:
        continue
    dominant_target[sid] = dom
    letter_counts_per_smorf[dom] += 1

annotated_smorf = len(annot_smorf)
unannotated_smorf = total_smorf - annotated_smorf
annot_ids = list(smorf_occ_annot.keys())

# %%
# -------- PER-smORF rate table --------
for sid in smorf_occ_annot:
    n = smorf_occ_annot[sid]
    m = smorf_match_any[sid]
    smorf_stats_rows.append({
        "smorf_id": sid,
        "n_annot_occ": n,
        "n_match_any_occ": m,
        "match_rate": (m / n) if n else 0.0,
        "dominant_letter": dominant_target.get(sid, "")
    })
df_smorf_rates = pd.DataFrame(smorf_stats_rows)

# %%
# -------- PRINT CORE SUMMARY --------
print(f"\nRoot: {ROOT}")
print(f"BIN smORFs (folders): {total_smorf}")
print(f"Annotated smORFs (>=1 annotated target occurrence): {annotated_smorf} ({pct(annotated_smorf, total_smorf):.2f}%)")
print(f"Unannotated smORFs (never annotated): {unannotated_smorf} ({pct(unannotated_smorf, total_smorf):.2f}%)")

print(f"\nTarget occurrences (all gffs with target=1 found): {occ_total}")
print(f"Annotated occurrences: {occ_annot} ({pct(occ_annot, occ_total):.2f}%)")
print(f"Unknown/unannotated occurrences: {occ_unknown} ({pct(occ_unknown, occ_total):.2f}%)")

print("\nCOG dominant letter per smORF (top)")
for i, (L, n) in enumerate(letter_counts_per_smorf.most_common(TOP_N), 1):
    print(f"{i:>2}. {L}\t{n}")

print("\nCOG letters per occurrence (top)")
for i, (L, n) in enumerate(letter_counts_per_occ.most_common(TOP_N), 1):
    print(f"{i:>2}. {L}\t{n}")

# %%
# -------- MATCH STRENGTH BINS (PER smORF) --------
bins = Counter()
for sid in annot_ids:
    rate = (smorf_match_any[sid] / smorf_occ_annot[sid]) if smorf_occ_annot[sid] else 0.0
    if rate == 0:
        bins["0%"] += 1
    elif rate <= 0.25:
        bins["1–25%"] += 1
    elif rate <= 0.50:
        bins["26–50%"] += 1
    elif rate <= 0.75:
        bins["51–75%"] += 1
    else:
        bins["76–100%"] += 1

smorf_any = sum(1 for sid in annot_ids if smorf_match_any[sid] > 0)
smorf_never = len(annot_ids) - smorf_any

print(f"\nAnnotated smORFs (neighbor-match analysis): {len(annot_ids)}")
print(f"SmORFs with >=1 neighbor match: {smorf_any} ({pct(smorf_any, len(annot_ids)):.2f}%)")
print(f"SmORFs with no match: {smorf_never} ({pct(smorf_never, len(annot_ids)):.2f}%)")

print("\nMatch strength (per smORF) [fraction of annotated occurrences that match]")
for k in ["0%", "1–25%", "26–50%", "51–75%", "76–100%"]:
    print(f"{k:>7}: {bins[k]} ({pct(bins[k], len(annot_ids)):.2f}%)")

print(f"\nAnnotated target occurrences (neighbor-match): {occ_annot_total}")
print(f"Occurrences with any neighbor match: {occ_match_any} ({pct(occ_match_any, occ_annot_total):.2f}%)")

# %%
# -------- YOUR L!=R RESULTS (vice versa) --------
smorf_has_LneqR = sum(1 for sid in annot_ids if smorf_LneqR_occ[sid] > 0)
smorf_any_left_only = sum(1 for sid in annot_ids if smorf_LneqR_left_only[sid] > 0)
smorf_any_right_only = sum(1 for sid in annot_ids if smorf_LneqR_right_only[sid] > 0)
smorf_any_both = sum(1 for sid in annot_ids if smorf_LneqR_both[sid] > 0)
smorf_any_neither = sum(1 for sid in annot_ids if smorf_LneqR_neither[sid] > 0)

print("\n============================")
print("SMORF COUNTS FOR 'NEIGHBOURS DIFFERENT' (L!=R)  [VICE VERSA]")
print("============================")
print(f"Annotated smORFs with >=1 L!=R occurrence: {smorf_has_LneqR} ({pct(smorf_has_LneqR, len(annot_ids)):.2f}%)")
print(f"smORFs that ever match LEFT only (in L!=R):   {smorf_any_left_only} ({pct(smorf_any_left_only, smorf_has_LneqR):.2f}%)")
print(f"smORFs that ever match RIGHT only (in L!=R):  {smorf_any_right_only} ({pct(smorf_any_right_only, smorf_has_LneqR):.2f}%)")
print(f"smORFs that ever match BOTH neighbours (in L!=R): {smorf_any_both} ({pct(smorf_any_both, smorf_has_LneqR):.2f}%)")
print(f"smORFs that ever match NEITHER (in L!=R):     {smorf_any_neither} ({pct(smorf_any_neither, smorf_has_LneqR):.2f}%)")

# %%
# -------- NEW: L==R but DIFFERENT -> which smORF functions are those? --------
smorf_has_LeqR = sum(1 for sid in annot_ids if smorf_LeqR_occ[sid] > 0)
smorf_any_LeqR_diff = sum(1 for sid in annot_ids if smorf_LeqR_diff[sid] > 0)

print("\n============================")
print("L==R CASES WHERE TARGET FUNCTION IS DIFFERENT")
print("============================")
print(f"Annotated smORFs with >=1 L==R occurrence: {smorf_has_LeqR} ({pct(smorf_has_LeqR, len(annot_ids)):.2f}%)")
print(f"smORFs that show at least one L==R DIFFERENT case: {smorf_any_LeqR_diff} ({pct(smorf_any_LeqR_diff, smorf_has_LeqR):.2f}%)")
print(f"Occurrence-level L==R annotated: {occ_LeqR_annot}")
print(f"Occurrence-level L==R match:     {occ_LeqR_match} ({pct(occ_LeqR_match, occ_LeqR_annot):.2f}%)")
print(f"Occurrence-level L==R DIFFERENT: {occ_LeqR_diff} ({pct(occ_LeqR_diff, occ_LeqR_annot):.2f}%)")

print("\nTop shared neighbour COG letter-sets in L==R DIFFERENT occurrences:")
for i, (x, n) in enumerate(LeqR_diff_neighbor_letters.most_common(TOP_N), 1):
    print(f"{i:>2}. {x}\t{n}")

print("\nTop target COG letters in L==R DIFFERENT occurrences (target letters counted individually):")
for i, (x, n) in enumerate(LeqR_diff_target_letters.most_common(TOP_N), 1):
    print(f"{i:>2}. {x}\t{n}")

# Pair table: (shared neighbour set) -> (target representative letter)
pair_rows = []
for (neigh, tdom), n in LeqR_diff_pair_counts.most_common(200):
    pair_rows.append({
        "shared_neighbour_letters": neigh,
        "target_rep_letter": tdom,
        "n_occurrences": n
    })
df_pairs = pd.DataFrame(pair_rows)
print("\nTop (neighbour_shared -> target_letter) pairs for L==R DIFFERENT:")
display(df_pairs.head(30))

# Also per-smORF summary: which smORFs contribute most to L==R DIFFERENT
smorf_diff_rows = []
for sid in sorted(smorf_LeqR_diff_smorfs):
    smorf_diff_rows.append({
        "smorf_id": sid,
        "dominant_target_letter_overall": dominant_target.get(sid, ""),
        "n_LeqR_occurrences": smorf_LeqR_occ[sid],
        "n_LeqR_different_occurrences": smorf_LeqR_diff[sid],
        "pct_LeqR_different": pct(smorf_LeqR_diff[sid], smorf_LeqR_occ[sid])
    })
df_smorf_diff = pd.DataFrame(smorf_diff_rows).sort_values(
    ["n_LeqR_different_occurrences", "pct_LeqR_different"],
    ascending=False
)
print("\nTop smORFs showing L==R DIFFERENT (per smORF):")
display(df_smorf_diff.head(30))

# %%
# -------- PLOTS --------

# 1) Annotated vs unannotated smORFs
plt.figure()
plt.bar(["Annotated", "Unannotated"], [annotated_smorf, unannotated_smorf])
plt.title("26–30 bin: smORFs annotated vs unannotated (per smORF)")
plt.ylabel("Number of smORFs")
plt.show()

# 2) Dominant letter per smORF (top N)
top_dom = letter_counts_per_smorf.most_common(TOP_N)
plt.figure()
plt.bar([x[0] for x in top_dom], [x[1] for x in top_dom])
plt.title(f"Dominant COG letter per smORF (Top {TOP_N})")
plt.xlabel("COG letter")
plt.ylabel("Number of smORFs")
plt.show()

# 3) Match strength bins (per smORF)
bin_order = ["0%", "1–25%", "26–50%", "51–75%", "76–100%"]
plt.figure()
plt.bar(bin_order, [bins[k] for k in bin_order])
plt.title("Match strength bins (per smORF): fraction of annotated occurrences with any neighbour match")
plt.xlabel("Bin")
plt.ylabel("Number of smORFs")
plt.show()

# 4) Scatter: n_annot_occ vs match_rate
if not df_smorf_rates.empty:
    plt.figure()
    plt.scatter(df_smorf_rates["n_annot_occ"], df_smorf_rates["match_rate"])
    plt.title("Per smORF: annotated occurrence count vs neighbour-match rate")
    plt.xlabel("Number of annotated occurrences (per smORF)")
    plt.ylabel("Match rate (any-neighbour match / annotated occurrences)")
    plt.show()

# 5) L!=R match-type presence (smORF level)
plt.figure()
plt.bar(["Left-only", "Right-only", "Both", "Neither"],
        [smorf_any_left_only, smorf_any_right_only, smorf_any_both, smorf_any_neither])
plt.title("smORF-level (L!=R): smORFs that ever show each match type")
plt.ylabel("Number of smORFs")
plt.show()

# 6) L==R: match vs different (occurrence-level)
plt.figure()
plt.bar(["L==R match", "L==R different"], [occ_LeqR_match, occ_LeqR_diff])
plt.title("Occurrence-level: when neighbours are same (L==R), does target match or differ?")
plt.ylabel("Number of occurrences")
plt.show()

# 7) L==R DIFFERENT: neighbour letters (top N)
top_neigh_diff = LeqR_diff_neighbor_letters.most_common(TOP_N)
if top_neigh_diff:
    plt.figure()
    plt.bar([x[0] for x in top_neigh_diff], [x[1] for x in top_neigh_diff])
    plt.title(f"L==R DIFFERENT: shared neighbour letter-sets (Top {TOP_N})")
    plt.xlabel("Shared neighbour COG letters")
    plt.ylabel("Number of occurrences")
    plt.show()

# 8) L==R DIFFERENT: target letters (top N)
top_target_diff = LeqR_diff_target_letters.most_common(TOP_N)
if top_target_diff:
    plt.figure()
    plt.bar([x[0] for x in top_target_diff], [x[1] for x in top_target_diff])
    plt.title(f"L==R DIFFERENT: target COG letters (Top {TOP_N})")
    plt.xlabel("Target COG letter")
    plt.ylabel("Number of occurrences")
    plt.show()

Root: /work/microbiome/users/kruthi/SmORF_neighbourhoods_26_30_no_min_overlap
BIN smORFs (folders): 1309
Annotated smORFs (>=1 annotated target occurrence): 644 (49.20%)
Unannotated smORFs (never annotated): 665 (50.80%)

Target occurrences (all gffs with target=1 found): 36298
Annotated occurrences: 16790 (46.26%)
Unknown/unannotated occurrences: 19508 (53.74%)

COG dominant letter per smORF (top)
 1. S	123
 2. J	116
 3. L	93
 4. K	83
 5. C	53
 6. G	44
 7. P	22
 8. U	20
 9. O	16
10. I	15
11. D	11
12. M	9
13. H	8
14. N	8
15. F	6
16. E	6
17. T	6
18. V	5

COG letters per occurrence (top)
 1. S	3285
 2. J	3243
 3. K	2210
 4. L	1877
 5. C	1491
 6. G	965
 7. P	672
 8. U	650
 9. O	527
10. Q	504
11. I	421
12. D	343
13. M	270
14. T	260
15. N	224
16. H	194
17. V	168
18. F	147
19. E	100
20. B	1

Annotated smORFs (neighbor-match analysis): 644
SmORFs with >=1 neighbor match: 320 (49.69%)
SmORFs with no match: 324 (50.31%)

Match strength (per smORF) [fraction of annotated occurrences that match]
     0%: 324 (50.31%)
  1–25%: 40 (6.21%)
 26–50%: 36 (5.59%)
 51–75%: 4 (0.62%)
76–100%: 240 (37.27%)

Annotated target occurrences (neighbor-match): 16790
Occurrences with any neighbor match: 6914 (41.18%)

============================
SMORF COUNTS FOR 'NEIGHBOURS DIFFERENT' (L!=R)  [VICE VERSA]
============================
Annotated smORFs with >=1 L!=R occurrence: 360 (55.90%)
smORFs that ever match LEFT only (in L!=R):   146 (40.56%)
smORFs that ever match RIGHT only (in L!=R):  145 (40.28%)
smORFs that ever match BOTH neighbours (in L!=R): 6 (1.67%)
smORFs that ever match NEITHER (in L!=R):     222 (61.67%)

============================
L==R CASES WHERE TARGET FUNCTION IS DIFFERENT
============================
Annotated smORFs with >=1 L==R occurrence: 132 (20.50%)
smORFs that show at least one L==R DIFFERENT case: 43 (32.58%)
Occurrence-level L==R annotated: 3014
Occurrence-level L==R match:     2156 (71.53%)
Occurrence-level L==R DIFFERENT: 858 (28.47%)

Top shared neighbour COG letter-sets in L==R DIFFERENT occurrences:
 1. L	239
 2. K	152
 3. S	135
 4. J	115
 5. P	55
 6. I	29
 7. F	28
 8. M	27
 9. O	26
10. E	23
11. G	19
12. C	10

Top target COG letters in L==R DIFFERENT occurrences (target letters counted individually):
 1. S	458
 2. J	86
 3. G	70
 4. K	61
 5. C	48
 6. L	32
 7. T	31
 8. U	30
 9. V	22
10. H	18
11. N	2

Top (neighbour_shared -> target_letter) pairs for L==R DIFFERENT:
