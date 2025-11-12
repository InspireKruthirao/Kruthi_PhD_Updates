
#!/usr/bin/env python3
import pandas as pd

# === File paths ===
input_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
target_smorf = "SHD.ORF.008_361_757"
output_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/smorf_neighbors.tsv"

# Number of neighbors on each side
n_neighbors = 5

# === Step 1: Load ORF annotation file ===
df = pd.read_csv(input_file, sep="\t", compression="xz", header=None,
                 names=["ORF_ID", "Sample", "Contig", "Start", "End", "Strand", "Partial"])

# Sort by contig and start coordinate for reliable ordering
df = df.sort_values(by=["Contig", "Start"]).reset_index(drop=True)

# === Step 2: Locate target SmORF ===
if target_smorf not in df["ORF_ID"].values:
    raise ValueError(f"Target SmORF {target_smorf} not found in {input_file}")

target_row = df[df["ORF_ID"] == target_smorf].iloc[0]
contig = target_row["Contig"]

# Restrict to same contig
same_contig = df[df["Contig"] == contig].reset_index(drop=True)

# Index of the target within the contig
target_idx = same_contig.index[same_contig["ORF_ID"] == target_smorf][0]

# === Step 3: Extract neighbors ===
start_idx = max(0, target_idx - n_neighbors)
end_idx = min(len(same_contig), target_idx + n_neighbors + 1)
region = same_contig.iloc[start_idx:end_idx].copy()

# === Step 4: Assign position labels ===
region["Position"] = None
for i, row in region.iterrows():
    offset = i - target_idx
    if offset == 0:
        region.loc[i, "Position"] = "Target"
    elif offset < 0:
        region.loc[i, "Position"] = f"Left_{abs(offset)}"
    else:
        region.loc[i, "Position"] = f"Right_{offset}"

# Reorder columns
region = region[["ORF_ID", "Sample", "Contig", "Start", "End", "Strand", "Partial", "Position"]]

# === Step 5: Save output ===
region.to_csv(output_file, sep="\t", index=False)

print(f"✅ SmORF neighborhood written to: {output_file}")
