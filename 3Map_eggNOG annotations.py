import pandas as pd
import os

# Load clusters
clusters = pd.read_csv(
    "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.clusters.tsv.xz",
    sep=r'\s+', header=None, compression='xz', engine='python'
)
orf_info = pd.read_csv(
    "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz",
    sep='\t', compression='xz',
    names=['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial'],
    dtype=str
)

# Function to get all ORFs in the same cluster
def get_cluster_orfs(orf_id):
    row = clusters[clusters[0] == orf_id].iloc[0]
    parts = row.dropna().tolist()
    cluster = parts[parts.index('R') + 1] if 'R' in parts else parts[-1]
    copies = clusters[clusters.iloc[:, -1] == cluster][0].tolist()
    result = orf_info[orf_info['ORF'].isin(copies)][['ORF', 'Sample', 'Contig']]
    result = result.sort_values(['Sample', 'Contig'])
    return result, cluster

# Target ORF
target_orf = "SHD.ORF.008_361_776"
orf_table, cluster_name = get_cluster_orfs(target_orf)

# Output path
BASE_DIR = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"
output_file = f"eggNOG_annotations_{cluster_name}.tsv"

# Write annotations
with open(output_file, "w") as out:
    out.write("#query\tORF\tSample\tContig\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
    for _, row in orf_table.iterrows():
        orf_id, sample_id, contig_id = row['ORF'], row['Sample'], row['Contig']
        annot_file = os.path.join(BASE_DIR, sample_id, f"{sample_id}.emapper.annotations")
        with open(annot_file, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                if contig_id in line:
                    out.write(f"{orf_id}\t{sample_id}\t{contig_id}\t{line}")
                    break

print(f"Annotations saved → {output_file}")
