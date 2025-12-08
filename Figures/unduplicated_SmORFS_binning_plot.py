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

# Parse smORF → length
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

# smORF → set of samples (unique per smORF)
smorf_to_samples = defaultdict(set)
with gzip.open(tsv_file, "rt") as f:
    next(f)  # skip header
    for line in f:
        smorf, sample = line.strip().split("\t")[:2]
        smorf_to_samples[smorf].add(sample)

# sample → set of unique smORFs
sample_to_smorfs = defaultdict(set)
for smorf, samples in smorf_to_samples.items():
    for sample in samples:
        sample_to_smorfs[sample].add(smorf)

# sample → list of smORF lengths (unique)
sample_lengths = {
    sample: [smorf_to_length[smorf] for smorf in smorfs]
    for sample, smorfs in sample_to_smorfs.items()
}

samples = sorted(sample_lengths.keys())

# Binning
bins = list(range(8, 104, 3))
bin_labels = [f"{i}-{i+2}" for i in range(8, 101, 3)]

data = []
for sample in samples:
    hist, _ = np.histogram(sample_lengths[sample], bins=bins)
    for label, count in zip(bin_labels, hist):
        data.append({"Sample": sample, "Length bin": label, "Count": int(count)})

df = pd.DataFrame(data)
df["Length bin"] = pd.Categorical(df["Length bin"], categories=bin_labels, ordered=True)

# Plot
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

# Mean diamonds
means = df.groupby("Length bin")["Count"].mean()
for i, mean_val in enumerate(means):
    ax.plot(i, mean_val, marker='D', markersize=9, color='black',
            markeredgecolor='white', markeredgewidth=1.2)

plt.ylabel("Number of unique smORFs per sample")
plt.xlabel("smORF length (amino acids)")
plt.xticks(rotation=90, ha='center')

# Pie chart inset
total_unique = len(smorf_to_length)
n_long  = sum(1 for l in smorf_to_length.values() if l > 50)
n_short = total_unique - n_long

pie_ax = plt.axes([0.30, 0.70, 0.19, 0.19])
wedges, _, autotexts = pie_ax.pie(
    [n_long, n_short],
    autopct='%1.1f%%',
    startangle=90,
    colors=[main_color, short_color],
    wedgeprops=dict(edgecolor='white'),
    textprops={'fontsize': 11, 'fontweight': 'bold', 'color': 'white'}
)

pie_ax.legend(
    wedges,
    [f">50 aa  ({n_long:,})", f"≤50 aa  ({n_short:,})"],
    loc="lower center",
    bbox_to_anchor=(0.5, -0.22),
    fontsize=11,
    frameon=False
)

plt.title(f"Total unique smORFs: {total_unique:,}", fontsize=10, loc='right', pad=10)
plt.tight_layout()
plt.subplots_adjust(left=0.22, bottom=0.15)
plt.savefig("smORF_length_distribution_unique.png", dpi=400, bbox_inches='tight')
plt.show()

# Summary printout
print(f"\nTotal unique smORFs: {total_unique:,}")
print(f">50 aa: {n_long:,} ({100*n_long/total_unique:.1f}%)")
print(f"≤50 aa: {n_short:,} ({100*n_short/total_unique:.1f}%)")
print("\nPer-sample unique smORF counts:")
for sample in samples:
    print(f"  {sample}: {len(sample_lengths[sample]):,}")
