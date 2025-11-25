import pandas as pd
import os
smorf_id = "SHD1_SM.100AA.003_862"   # CHANGE 
output_root = "/work/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods"

smorf_origins_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
orf_catalog_file   = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
clusters_file      = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.clusters.tsv.xz"
eggNOG_base        = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"

print("Loading SmORF origins...")
smorf_df = pd.read_csv(
    smorf_origins_file,
    sep="\t",
    compression="gzip",
    names=["SmORF_ID", "Sample", "Contig", "Coordinates", "Strand"],
    dtype=str
)

print("Loading ORF catalog...")
orf_df = pd.read_csv(
    orf_catalog_file,
    sep="\t",
    compression="xz",
    names=["ORF", "Sample", "Contig", "Start", "End", "Strand", "Partial"],
    dtype=str
)

print("Loading cluster table...")
clusters = pd.read_csv(
    clusters_file,
    sep=r"\s+",
    header=None,
    compression="xz",
    engine="python"
)

### remove _polypolish suffix
def clean_contig(c):
    return c.split("_polypolish")[0]

## find cluster of ORF
def get_cluster_for_orf(orf_id):
    row = clusters[clusters[0] == orf_id].iloc[0]
    parts = row.dropna().tolist()
    # cluster is = last column unless "R" present
    cluster = parts[parts.index("R") + 1] if "R" in parts else parts[-1]
    return cluster

# get all ORFs in that cluster
def get_all_orfs_in_cluster(cluster_name):
    all_orfs = clusters[clusters.iloc[:, -1] == cluster_name][0].tolist()
    return all_orfs
occ = smorf_df[smorf_df["SmORF_ID"] == smorf_id]

if occ.empty:
    raise ValueError(f"No occurrences of {smorf_id} found.")

print(f"Found {len(occ)} occurrences of {smorf_id}.")

smorf_out = os.path.join(output_root, smorf_id)
os.makedirs(smorf_out, exist_ok=True)

for i, row in occ.iterrows():
    sample = row["Sample"]
    contig_raw = row["Contig"]
    contig_clean = clean_contig(contig_raw)
    start_s, end_s = map(int, row["Coordinates"].split("-"))

    print(f"\n=== Processing: {sample} / {contig_clean} ===")

    # Create subfolder
    occ_dir = os.path.join(smorf_out, f"{sample}_{contig_clean}")
    os.makedirs(occ_dir, exist_ok=True)

    # Filter ORFs matching sample & contig
    subset = orf_df[
        (orf_df["Sample"] == sample) &
        (orf_df["Contig"].str.startswith(contig_clean))
    ].copy()

    subset["Start"] = subset["Start"].astype(int)
    subset["End"] = subset["End"].astype(int)
    subset = subset.sort_values("Start").reset_index(drop=True)

    # Find ORF matching coordinates
    match = subset[(subset["Start"] == start_s) & (subset["End"] == end_s)]

    if match.empty:
        print("WARNING: No ORF matched coordinates. Skipping this occurrence.")
        continue

    target_orf = match["ORF"].iloc[0]
    idx = match.index[0]

    print(f"Target ORF found: {target_orf}")

    ### Extract ±5 neighbours
    L = max(0, idx - 5)
    R = min(len(subset), idx + 6)

    neighbors = subset.iloc[L:R].copy()
    neighbors["Note"] = ""
    neighbors.loc[neighbors["ORF"] == target_orf, "Note"] = f"TARGET ({smorf_id})"

    neighbors_outfile = os.path.join(occ_dir, "neighbors.tsv")
    neighbors.to_csv(neighbors_outfile, sep="\t", index=False)
    print(f"Saved neighbours → {neighbors_outfile}")
    
    ### CLUSTER + EggNOG annotation
    egg_file = os.path.join(occ_dir, "EggNOG_annotations.tsv")
    with open(egg_file, "w") as OUT:
        OUT.write(
            "SmORF_ID\tOccurrence_ORF\tCluster\tCopy_ORF\tSample\tContig\t" +
            "seed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t" +
            "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t" +
            "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t" +
            "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
        )

        for _, nrow in neighbors.iterrows():
            orf_id = nrow["ORF"]

            # Find cluster
            cluster_name = get_cluster_for_orf(orf_id)
            copies = get_all_orfs_in_cluster(cluster_name)

            copies_info = orf_df[orf_df["ORF"].isin(copies)]

            # For each copy, extract EggNOG annotation
            for _, crow in copies_info.iterrows():
                smp = crow["Sample"]
                cont = crow["Contig"]

                annot_file = os.path.join(
                    eggNOG_base, smp, f"{smp}.emapper.annotations"
                )

                if not os.path.exists(annot_file):
                    continue

                with open(annot_file) as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        if cont in line:
                            OUT.write(
                                f"{smorf_id}\t{orf_id}\t{cluster_name}\t" +
                                f"{crow['ORF']}\t{smp}\t{cont}\t{line}"
                            )
                            break

    print(f"Saved EggNOG annotations → {egg_file}")

print("\nDONE.")
