import pandas as pd
import importlib

clusters_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.clusters.tsv.xz"
target_orf = "SHD.ORF.008_361_776"

print(f"Searching for 95NT cluster of {target_orf} ...")

df = pd.read_csv(
    clusters_file,
    sep="\t",          
    compression="xz",
    header=None,
    low_memory=False
)

row = df[df.iloc[:, 0] == target_orf]

if not row.empty:
    parts = row.iloc[0].dropna().tolist()
    if "R" in parts:
        cluster_95nt = parts[parts.index("R") + 1]
    else:
        cluster_95nt = parts[-1]
    print(f"{target_orf} → {cluster_95nt}")
else:
    raise ValueError(f"{target_orf} not found in clusters file!")

# Find all ORFs in the same 95NT cluster
mask = df.iloc[:, -1] == cluster_95nt
orfs_in_cluster = df[mask].iloc[:, 0].tolist()

print(f"\nFound {len(orfs_in_cluster)} ORFs in cluster {cluster_95nt}:\n")
for orf in orfs_in_cluster:
    print(orf)
