from pathlib import Path
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = Path("/work/microbiome/users/kruthi/SmORF_neighbourhoods_26_30")
OUT = ROOT / "plots"
OUT.mkdir(parents=True, exist_ok=True)
OUT_PNG = OUT / "cog_dominant_letter_per_smorf.png"

def parse_gff_target_lines(gff_path: Path):
    with gff_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            if "target=1" not in attrs:
                continue
            name_val = None
            for item in attrs.split(";"):
                if item.startswith("Name="):
                    name_val = item.split("=", 1)[1]
                    break
            if not name_val or name_val == "Unknown" or not name_val.startswith("COG") or "-" not in name_val:
                continue
            cat = name_val.rsplit("-", 1)[1].strip()
            yield [ch for ch in cat if ch.isalpha()]

def find_smorf_dirs(root: Path):
    for d in sorted([p for p in root.iterdir() if p.is_dir()]):
        if any(child.is_dir() and "_contig_" in child.name and child.name.startswith("D") for child in d.iterdir()):
            yield d

def find_gffs_in_smorf_dir(smorf_dir: Path):
    for contig_dir in smorf_dir.iterdir():
        if contig_dir.is_dir() and "_contig_" in contig_dir.name:
            yield from contig_dir.glob("*.gff")

letters_per_smorf = defaultdict(Counter)
letter_counts_per_smorf = Counter()

for smorf_dir in find_smorf_dirs(ROOT):
    smorf_id = smorf_dir.name
    for gff_path in find_gffs_in_smorf_dir(smorf_dir):
        for letters in parse_gff_target_lines(gff_path):
            letters_per_smorf[smorf_id].update(letters)

for smorf_id, ctr in letters_per_smorf.items():
    if not ctr:
        continue
    m = max(ctr.values())
    dom = sorted([L for L, c in ctr.items() if c == m])[0]
    letter_counts_per_smorf[dom] += 1

items = letter_counts_per_smorf.most_common()
letters = [k for k, _ in items][::-1]
counts = [v for _, v in items][::-1]

cmap = mpl.cm.get_cmap("Dark2")
colors = [cmap(i % cmap.N) for i in range(len(letters))]

plt.figure(figsize=(8, max(4, 0.35 * len(letters))))
plt.barh(letters, counts, color=colors)
plt.xlabel("smORFs")
plt.ylabel("COG")
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
plt.show()

print(OUT_PNG)
