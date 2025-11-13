import pandas as pd
#paths
clusters_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.clusters.tsv.xz"
orf_file      = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"

#data
clusters = pd.read_csv(clusters_file, sep=r'\s+', header=None, compression='xz', low_memory=False)
orf_info = pd.read_csv(orf_file, sep='\t', compression='xz',
                       names=['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial'],
                       dtype=str, low_memory=False)

# Target ORF
target_orf = "SHD.ORF.008_361_776"

# Find cluster of target ORF
row = clusters[clusters[0] == target_orf].iloc[0]
parts = row.dropna().tolist()
cluster = parts[parts.index('R') + 1] if 'R' in parts else parts[-1]

# Get all ORFs in the same cluster
copies = clusters[clusters.iloc[:, -1] == cluster][0].tolist()

# Get ORF info for these copies
cluster_df = orf_info[orf_info['ORF'].isin(copies)][['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand']]
cluster_df = cluster_df.sort_values(['Sample', 'Contig', 'Start'])

# Display
cluster_df
