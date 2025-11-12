import pandas as pd
import os
import re

print("Loading data...")
clusters = pd.read_csv(
    "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.clusters.tsv.xz",
    sep=r'\s+', header=None, compression='xz', engine='python'
)
orf_info = pd.read_csv(
    "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz",
    sep='\t', compression='xz',
    names=['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial'],
    dtype=str, low_memory=False
)

def get_cluster_orfs(orf_id):
    row = clusters[clusters[0] == orf_id]
    if row.empty:
        return None, None
    parts = row.iloc[0].tolist()
    cluster = parts[parts.index('R') + 1] if 'R' in parts else parts[-1]
    copies = clusters[clusters.iloc[:, -1] == cluster][0].tolist()
    result = orf_info[orf_info['ORF'].isin(copies)][['ORF', 'Sample', 'Contig']].copy()
    result = result.sort_values(['Sample', 'Contig'])
    print(f"{orf_id} → {cluster} | {len(copies)} copies found")
    return result, cluster

orf_table, cluster_name = get_cluster_orfs("SHD.ORF.008_361_776")

BASE_DIR = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"
output_file = f"eggNOG_annotations_{cluster_name}.tsv"

hits = []
missing = []

print(f"\nSearching eggNOG annotations for {len(orf_table)} ORFs...")

with open(output_file, "w") as out:
    out.write("#query\tORF\tSample\tContig\tseed_ortholog\tevalue\tscore\t eggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
    
    for _, row in orf_table.iterrows():
        orf_id = row['ORF']
        sample_id = row['Sample']
        contig_id = row['Contig']
        
        annot_file = os.path.join(BASE_DIR, sample_id, f"{sample_id}.emapper.annotations")
        
        if not os.path.exists(annot_file):
            missing.append(sample_id)
            continue
        
        found = False
        with open(annot_file, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                if contig_id in line:
                    fields = line.strip().split("\t")
                    if len(fields) >= 11:
                        out.write(f"{orf_id}\t{sample_id}\t{contig_id}\t" + "\t".join(fields) + "\n")
                        hits.append(orf_id)
                        found = True
                        break
        if not found:
            out.write(f"{orf_id}\t{sample_id}\t{contig_id}\tNO_ANNOTATION_IN_CONTIG\n")

print(f"\nDone!")
print(f"eggNOG hits found: {len(hits)} / {len(orf_table)} ORFs")
print(f"Saved full annotations → {output_file}")

if missing:
    print(f"Warning: Missing eggNOG files for samples: {sorted(set(missing))}")
