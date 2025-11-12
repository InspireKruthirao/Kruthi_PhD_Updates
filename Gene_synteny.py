import pandas as pd

input_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
target_smorf = "SHD.ORF.008_361_757"
output_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/smorf_neighbors.tsv"

#Number of neighbors on each side
n_neighbors = 5

# Load ORF file ===
df = pd.read_csv(
    input_file,
    sep="\t",
    compression="xz",
    dtype={"Start": int, "End": int, "Strand": int, "Partial": str},
    low_memory=False
)

# contig_133 
df_contig = df[df["Original_ID"].str.startswith("contig_133")].copy()

# start coordinate for ordering
df_contig = df_contig.sort_values(by="Start").reset_index(drop=True)

# Locate target SmORF 
if target_smorf not in df_contig["ORF"].values:
    raise ValueError(f"Target SmORF {target_smorf} not found in contig_133")

target_idx = df_contig.index[df_contig["ORF"] == target_smorf][0]

# Extract neighbors 
start_idx = max(0, target_idx - n_neighbors)
end_idx = min(len(df_contig), target_idx + n_neighbors + 1)
region = df_contig.iloc[start_idx:end_idx].copy()

# Assign position labels
region["Position"] = None
for i, row in region.iterrows():
    offset = i - target_idx
    if offset == 0:
        region.loc[i, "Position"] = "Target"
    elif offset < 0:
        region.loc[i, "Position"] = f"Left_{abs(offset)}"
    else:
        region.loc[i, "Position"] = f"Right_{offset}"


region = region[["ORF", "Sample", "Original_ID", "Start", "End", "Strand", "Partial", "Position"]]

#  display output
region.to_csv(output_file, sep="\t", index=False)

print(f"✅ SmORF neighborhood written to: {output_file}\n")
print("--- SmORF Neighborhood ---")
print(region.to_string(index=False))

