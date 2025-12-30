import gzip
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns

CATALOG_DIR = Path("/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog")
fasta_file = CATALOG_DIR / "100AA_SmORFs_sequences.faa.gz"
tsv_file   = CATALOG_DIR / "100AA_SmORFs_origins.tsv.gz"

smorf_to_length = {}
with gzip.open(fasta_file, "rt") as f:
    seq = ""
    cid = None
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if cid:
                smorf_to_length[cid] = len(seq.replace("*", ""))
            cid = line[1:].split()[0]
            seq = ""
        else:
            seq += line.replace("*", "")
    if cid:
        smorf_to_length[cid] = len(seq.replace("*", ""))
smorf_to_samples = defaultdict(set)
with gzip.open(tsv_file, "rt") as f:
    next(f)  # header
    for line in f:
        smorf, sample = line.strip().split("\t")[:2]
        smorf_to_samples[smorf].add(sample)

sample_lengths = defaultdict(list)
for smorf_id, samples in smorf_to_samples.items():
    length = smorf_to_length.get(smorf_id)
    if length:
        for s in samples:
            sample_lengths[s].append(length)

samples = sorted(sample_lengths.keys())
all_lengths = [l for lst in sample_lengths.values() for l in lst]

if not all_lengths:
    raise ValueError("No smORF lengths found. Check FASTA/TSV paths and IDs.")

max_len = max(all_lengths)   # e.g., 99 (no misleading '...-100' label)

start = 8
step = 3

last_start = start + ((max_len - start) // step) * step
bins = list(range(start, last_start + step + 1, step))

bin_labels = [f"{i}-{min(i+2, max_len)}" for i in range(start, last_start + 1, step)]

data = []
for sample in samples:
    hist, _ = np.histogram(sample_lengths[sample], bins=bins)
    for label, count in zip(bin_labels, hist):
        data.append({"Sample": sample, "Length bin": label, "Count": int(count)})

df = pd.DataFrame(data)
df["Length bin"] = pd.Categorical(df["Length bin"], categories=bin_labels, ordered=True)


plt.figure(figsize=(15, 8))
sns.set_style("whitegrid")

main_color  = sns.color_palette("Dark2")[1]
short_color = sns.color_palette("Dark2")[2]

ax = sns.boxplot(
    data=df, x="Length bin", y="Count",
    color=main_color, showfliers=False,
    boxprops=dict(edgecolor="black"),
    medianprops=dict(color="white", linewidth=2),
    whiskerprops=dict(color="black"),
    capprops=dict(color="black")
)


means = df.groupby("Length bin")["Count"].mean()
for i, mean_val in enumerate(means):
    ax.plot(
        i, mean_val,
        marker="D", markersize=9, color="black",
        markeredgecolor="white", markeredgewidth=1.2
    )

plt.ylabel("Number of predicted smORFs per sample")
plt.xlabel("smORF length (amino acids)")
plt.xticks(rotation=90, ha="center")


n_long  = sum(1 for l in all_lengths if l > 50)
n_short = sum(1 for l in all_lengths if l <= 50)

pie_ax = plt.axes([0.30, 0.70, 0.19, 0.19])
wedges, texts, autotexts = pie_ax.pie(
    [n_long, n_short],
    labels=["", ""],
    autopct="%1.1f%%",
    startangle=90,
    colors=[main_color, short_color],
    wedgeprops=dict(edgecolor="white"),
    textprops={"fontsize": 11, "fontweight": "bold", "color": "white"}
)

pie_ax.legend(
    wedges,
    [f">50 aa  ({n_long:,})", f"≤50 aa  ({n_short:,})"],
    loc="lower center",
    bbox_to_anchor=(0.5, -0.22),
    fontsize=11,
    frameon=False
)

plt.tight_layout()
plt.subplots_adjust(left=0.22, bottom=0.15)
plt.savefig("smORF_length_distribution_FULL_PIE_CHART.png", dpi=400, bbox_inches="tight")
plt.show()
