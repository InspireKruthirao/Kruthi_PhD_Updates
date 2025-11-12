import pandas as pd

smorf_id      = "SHD1_SM.100AA.003_862"   
target_sample = "D044"                    
target_contig = "contig_133_polypolish"  

smorf_file    = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
gene_catalog  = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"

# Load smORF catalog
smorf_df = pd.read_csv(smorf_file, sep='\t', compression='gzip',
                       names=['SmORF_ID', 'Sample', 'Contig', 'Coordinates', 'Strand'])

# Get target smORF coordinates
row = smorf_df[
    (smorf_df['SmORF_ID'] == smorf_id) &
    (smorf_df['Sample'] == target_sample) &
    (smorf_df['Contig'].str.startswith(target_contig))
].iloc[0]

start, end = map(int, row['Coordinates'].split('-'))

# Load gene catalog
df = pd.read_csv(gene_catalog, sep='\t', compression='xz', header=None,
                 names=['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial'],
                 usecols=[0,1,2,3,4,5,6])

# Filter for target contig and sort
contig_df = df[
    (df['Sample'] == target_sample) &
    (df['Contig'].str.startswith(target_contig))
].sort_values('Start').reset_index(drop=True)

# Find matching ORF
target_orf_row = contig_df[
    (contig_df['Start'] == start) &
    (contig_df['End'] == end)
]

idx = target_orf_row.index[0]
orf_id = target_orf_row['ORF'].values[0]

# Extract ±5 ORFs around target
left = max(0, idx - 5)
right = min(len(contig_df), idx + 6)
neighbors = contig_df.iloc[left:right].copy()

# Standardize contig name and mark target
neighbors['Contig'] = "contig_133"
neighbors['Note'] = ""
neighbors.loc[neighbors['ORF'] == orf_id, 'Note'] = f"TARGET_smORF ({smorf_id})"

# Display final table
neighbors[['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial', 'Note']]
