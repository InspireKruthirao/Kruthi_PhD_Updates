from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path("/work/microbiome/users/kruthi/SmORF_neighbourhoods_26_30")
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
    if not name or name == "Unknown" or not name.startswith("COG"):
        return set()
    i = name.rfind("-")
    if i == -1:
        return set()
    return set(ch for ch in name[i + 1:] if ch.isalpha())

def iter_smorf_gffs(root: Path):
    for smorf in root.iterdir():
        if not smorf.is_dir():
            continue
        sid = smorf.name
        for contig in smorf.iterdir():
            if contig.is_dir() and contig.name.startswith("D") and "_contig_" in contig.name:
                for gff in contig.glob("*.gff"):
                    yield sid, gff

def read_target_neighbors(gff: Path):
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

total_smorf = sum(1 for p in ROOT.iterdir() if p.is_dir())
annot_smorf = set()
letters_per_smorf = defaultdict(Counter)

occ_total = 0
occ_annot = 0
occ_unknown = 0
letter_counts_per_occ = Counter()

occ_annot_total = 0
occ_match_any = 0
smorf_occ_annot = Counter()
smorf_match_any = Counter()
target_votes = defaultdict(Counter)

for sid, gff in iter_smorf_gffs(ROOT):
    rec = read_target_neighbors(gff)
    if rec is None:
        continue

    t_set, up_set, dn_set = rec
    occ_total += 1

    if not t_set:
        occ_unknown += 1
        continue

    occ_annot += 1
    annot_smorf.add(sid)

    for L in t_set:
        letter_counts_per_occ[L] += 1
        letters_per_smorf[sid][L] += 1

    occ_annot_total += 1
    smorf_occ_annot[sid] += 1
    for L in t_set:
        target_votes[sid][L] += 1

    if (t_set & up_set) or (t_set & dn_set):
        occ_match_any += 1
        smorf_match_any[sid] += 1

letter_counts_per_smorf = Counter()
for sid, ctr in letters_per_smorf.items():
    if not ctr:
        continue
    m = max(ctr.values())
    dom = min([L for L, c in ctr.items() if c == m])
    letter_counts_per_smorf[dom] += 1

annotated_smorf = len(annot_smorf)
unannotated_smorf = total_smorf - annotated_smorf

def pct(a, b):
    return (a / b * 100) if b else 0.0

print(f"\nRoot: {ROOT}")
print(f"Total smORFs: {total_smorf}")
print(f"Annotated smORFs: {annotated_smorf} ({pct(annotated_smorf, total_smorf):.2f}%)")
print(f"Unannotated smORFs: {unannotated_smorf} ({pct(unannotated_smorf, total_smorf):.2f}%)")

print(f"\nTarget occurrences: {occ_total}")
print(f"Annotated occurrences: {occ_annot} ({pct(occ_annot, occ_total):.2f}%)")
print(f"Unknown/unannotated occurrences: {occ_unknown} ({pct(occ_unknown, occ_total):.2f}%)")

print("\nCOG dominant letter per smORF (top)")
for i, (L, n) in enumerate(letter_counts_per_smorf.most_common(TOP_N), 1):
    print(f"{i:>2}. {L}\t{n}")

print("\nCOG letters per occurrence (top)")
for i, (L, n) in enumerate(letter_counts_per_occ.most_common(TOP_N), 1):
    print(f"{i:>2}. {L}\t{n}")

annot_ids = list(smorf_occ_annot.keys())
bins = Counter()
for sid in annot_ids:
    rate = smorf_match_any[sid] / smorf_occ_annot[sid]
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

dominant_target = {}
for sid in annot_ids:
    ctr = target_votes[sid]
    mx = max(ctr.values())
    dominant_target[sid] = min([L for L, c in ctr.items() if c == mx])

letter_total = Counter()
letter_match = Counter()
for sid in annot_ids:
    L = dominant_target[sid]
    letter_total[L] += 1
    if smorf_match_any[sid] > 0:
        letter_match[L] += 1

smorf_any = sum(1 for sid in annot_ids if smorf_match_any[sid] > 0)
smorf_never = len(annot_ids) - smorf_any

print(f"\nAnnotated smORFs (neighbor-match analysis): {len(annot_ids)}")
print(f"SmORFs with >=1 neighbor match: {smorf_any} ({pct(smorf_any, len(annot_ids)):.2f}%)")
print(f"SmORFs with no match: {smorf_never} ({pct(smorf_never, len(annot_ids)):.2f}%)")

print("\nMatch strength (per smORF)")
for k in ["0%", "1–25%", "26–50%", "51–75%", "76–100%"]:
    print(f"{k:>7}: {bins[k]} ({pct(bins[k], len(annot_ids)):.2f}%)")

print(f"\nAnnotated target occurrences (neighbor-match): {occ_annot_total}")
print(f"Occurrences with any neighbor match: {occ_match_any} ({pct(occ_match_any, occ_annot_total):.2f}%)")

print("\nBy target dominant letter (per smORF)")
for i, (L, n) in enumerate(letter_total.most_common(TOP_N), 1):
    m = letter_match[L]
    print(f"{i:>2}. {L}\tTotal: {n}\tMatch>=1: {m} ({pct(m, n):.2f}%)")
