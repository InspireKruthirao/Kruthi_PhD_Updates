import pandas as pd

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

print("Ready.")

def get_cluster_copies(orf_id):
    row = clusters[clusters[0] == orf_id]
    if row.empty: 
        return print("Not found")
    
    parts = row.iloc[0].tolist()
    cluster = parts[parts.index('R') + 1] if 'R' in parts else parts[-1]
    copies = clusters[clusters.iloc[:, -1] == cluster][0].tolist()
    
    result = orf_info[orf_info['ORF'].isin(copies)][['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand']]
    result = result.sort_values(['Sample', 'Contig', 'Start'])
    
    result.to_csv(f"cluster_{cluster}.tsv", sep='\t', index=False)
    print(f"{orf_id} → {cluster} | {len(copies)} copies → saved")
    
    return result

get_cluster_copies("SHD.ORF.008_361_776")
