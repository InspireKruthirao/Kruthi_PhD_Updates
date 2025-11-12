import pandas as pd
import os

smorf_id         = "SHD1_SM.100AA.003_862"   
target_sample    = "D044"                    
target_contig    = "contig_133_polypolish"  

smorf_file       = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
gene_catalog     = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
output_file      = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/smorf_neighbors.tsv"


print(f"Looking for smORF: {smorf_id}")
smorf_df = pd.read_csv(smorf_file, sep='\t', compression='gzip',
                       names=['SmORF_ID', 'Sample', 'Contig', 'Coordinates', 'Strand'])

# Filter target
match = smorf_df[
    (smorf_df['SmORF_ID'] == smorf_id) &
    (smorf_df['Sample'] == target_sample) &
    (smorf_df['Contig'].str.startswith(target_contig))
]

if len(match) == 0:
    print("Not found! Available occurrences:")
    print(smorf_df[smorf_df['SmORF_ID'] == smorf_id][['Sample', 'Contig']].to_string(index=False))
    exit(1)
elif len(match) > 1:
    print("Multiple hits — picking first one")

row = match.iloc[0]
coords = row['Coordinates']
start, end = map(int, coords.split('-'))
strand = row['Strand']

print(f"Found in {target_sample} {row['Contig']} {coords} {strand}")

print("Loading full gene catalog...")
df = pd.read_csv(gene_catalog, sep='\t', compression='xz', header=None,
                 names=['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial'],
                 usecols=[0,1,2,3,4,5,6])

contig_mask = (df['Sample'] == target_sample) & (df['Contig'].str.startswith(target_contig))
contig_df = df[contig_mask].sort_values('Start').reset_index(drop=True)

target_orf_row = contig_df[
    (contig_df['Start'] == start) &
    (contig_df['End'] == end)
]

if len(target_orf_row) == 0:
    print("smORF not found in gene catalog! Possible mismatch.")
    exit(1)

idx = target_orf_row.index[0]
orf_id = target_orf_row['ORF'].values[0]
print(f"Matches ORF: {orf_id} at index {idx}")


left = max(0, idx - 5)
right = min(len(contig_df), idx + 6)
neighbors = contig_df.iloc[left:right].copy()

#contig name
neighbors['Contig'] = "contig_133"

# target
neighbors['Note'] = ""
neighbors.loc[neighbors['ORF'] == orf_id, 'Note'] = f"TARGET_smORF ({smorf_id})"

#Reorder
neighbors = neighbors[['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial', 'Note']]

#Save
neighbors.to_csv(output_file, sep='\t', index=False)

print(f"\nsaved! {len(neighbors)} ORFs saved to:")
print(output_file)
print("\nNeighborhood:")
print(neighbors.to_string(index=False))
