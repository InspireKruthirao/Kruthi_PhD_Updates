import gzip
import pandas as pd
from collections import Counter
import numpy as np

base_path = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/"
sequences_file = base_path + "100AA_SmORFs_sequences.faa.gz"
origins_file = base_path + "100AA_SmORFs_origins.tsv.gz"

print("="*80)
print("SHANGHAI DOGS smORF CATALOG ANALYSIS")
print("="*80)

print("\n### SEQUENCE STATISTICS ###\n")

smorf_lengths = {}
with gzip.open(sequences_file, 'rt') as f:
    current_id = None
    current_seq = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                seq = ''.join(current_seq).replace('*', '')
                smorf_lengths[current_id] = len(seq)
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        seq = ''.join(current_seq).replace('*', '')
        smorf_lengths[current_id] = len(seq)

total_unique_smorfs = len(smorf_lengths)
lengths = list(smorf_lengths.values())

print(f"Total unique smORFs: {total_unique_smorfs:,}")
print(f"Length range: {min(lengths)} - {max(lengths)} aa")
print(f"Mean length: {np.mean(lengths):.1f} ± {np.std(lengths):.1f} aa")
print(f"Median length: {np.median(lengths):.0f} aa")

length_50_or_less = sum(1 for l in lengths if l <= 50)
length_greater_50 = sum(1 for l in lengths if l > 50)

print(f"\nsmORFs ≤50 aa: {length_50_or_less:,} ({100*length_50_or_less/total_unique_smorfs:.1f}%)")
print(f"smORFs >50 aa: {length_greater_50:,} ({100*length_greater_50/total_unique_smorfs:.1f}%)")

length_counter = Counter(lengths)
most_common = length_counter.most_common(5)
least_common = length_counter.most_common()[:-6:-1]

print(f"\nFive most common lengths: {', '.join(str(l) for l, c in most_common)} aa")
print(f"Five least common lengths: {', '.join(str(l) for l, c in least_common)} aa")

print("\n### GENOMIC ORIGINS ANALYSIS ###\n")

df_origins = pd.read_csv(origins_file, sep='\t', compression='gzip')

print(f"Total smORF occurrences: {len(df_origins):,}")

total_samples = df_origins['Sample ID'].nunique()
print(f"Total samples: {total_samples}")

smorfs_per_sample = df_origins.groupby('Sample ID')['SmORF ID'].nunique()
print(f"\nsmORFs per sample:")
print(f"  Mean: {smorfs_per_sample.mean():.1f} ± {smorfs_per_sample.std():.1f}")
print(f"  Range: {smorfs_per_sample.min():,} - {smorfs_per_sample.max():,}")

smorf_samples = df_origins.groupby('SmORF ID')['Sample ID'].nunique()
sample_prevalence = smorf_samples.value_counts().sort_index()

unique_to_one = sample_prevalence.get(1, 0)
present_in_2plus = smorf_samples[smorf_samples >= 2].count()
core_smorfs = smorf_samples[smorf_samples == total_samples].count()

print(f"\nSample prevalence:")
print(f"  Unique to one sample: {unique_to_one:,} ({100*unique_to_one/total_unique_smorfs:.1f}%)")
print(f"  Present in ≥2 samples: {present_in_2plus:,} ({100*present_in_2plus/total_unique_smorfs:.1f}%)")
print(f"  Core (all samples): {core_smorfs:,}")

print("\n### CORE smORF CHARACTERISTICS ###\n")

core_smorf_ids = smorf_samples[smorf_samples == total_samples].index.tolist()
core_lengths = [smorf_lengths[sid] for sid in core_smorf_ids if sid in smorf_lengths]

print(f"Number of core smORFs: {len(core_lengths)}")
print(f"Length range: {min(core_lengths)} - {max(core_lengths)} aa")
print(f"Mean length: {np.mean(core_lengths):.1f} ± {np.std(core_lengths):.1f} aa")
print(f"Median length: {np.median(core_lengths):.0f} aa")

print("\n" + "="*80)
