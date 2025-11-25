#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os

SMORF_FILE      = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
GENE_CATALOG    = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
EGGNOG_DIR      = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"
OUTPUT_DIR      = "/work/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods"

smorf_df = pd.read_csv(SMORF_FILE, sep='\t', compression='gzip',
                       names=['SmORF_ID', 'Sample', 'Contig', 'Coordinates', 'Strand'])
gene_df = pd.read_csv(GENE_CATALOG, sep='\t', compression='xz',
                      names=['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial'], dtype=str)

counts = smorf_df['SmORF_ID'].value_counts()
selected_smorfs = counts[(counts >= 25) & (counts <= 30)].sample(5, random_state=42).index.tolist()

def get_neighbours(sample, base_contig, smorf_start, smorf_end, window=5, smorf_id=""):
    df = gene_df[
        (gene_df['Sample'] == sample) &
        (gene_df['Contig'].str.split('_polypolish').str[0] == base_contig)
    ].copy()
    if df.empty:
        return pd.DataFrame()

    df['Start'] = df['Start'].astype(int)
    df['End']   = df['End'].astype(int)
    df = df.sort_values('Start').reset_index(drop=True)

    df['overlap'] = (np.minimum(df['End'], smorf_end) - np.maximum(df['Start'], smorf_start)).clip(lower=0)
    if df['overlap'].max() < 80:
        return pd.DataFrame()

    best_idx = df['overlap'].idxmax()
    idx = best_idx
    left = max(0, idx - window)
    right = min(len(df), idx + window + 1)

    neigh = df.iloc[left:right].copy()
    neigh['Note'] = ""
    neigh.loc[best_idx, 'Note'] = f"TARGET ({smorf_id})"
    return neigh.drop(columns='overlap')

def extract_eggnog(df, smorf_id, path):
    with open(path, "w") as out:
        out.write("#query\tORF\tSample\tContig\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
        for _, row in df.iterrows():
            orf, sample, contig = row['ORF'], row['Sample'], row['Contig']
            annot_file = os.path.join(EGGNOG_DIR, sample, f"{sample}.emapper.annotations")
            if not os.path.exists(annot_file):
                continue
            with open(annot_file) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    if contig in line:
                        out.write(f"{smorf_id}\t{orf}\t{sample}\t{contig}\t{line}")
                        break

os.makedirs(OUTPUT_DIR, exist_ok=True)

for smorf_id in selected_smorfs:
    for _, row in smorf_df[smorf_df['SmORF_ID'] == smorf_id].iterrows():
        sample = row['Sample']
        base_contig = row['Contig'].split("_polypolish")[0]
        start, end = map(int, row['Coordinates'].split('-'))
        out_dir = os.path.join(OUTPUT_DIR, smorf_id, f"{sample}_{base_contig}")
        os.makedirs(out_dir, exist_ok=True)

        neigh = get_neighbours(sample, base_contig, start, end, window=5, smorf_id=smorf_id)
        if neigh.empty:
            continue

        neigh[['ORF','Sample','Contig','Start','End','Strand','Partial','Note']] \
            .to_csv(f"{out_dir}/neighbors.tsv", sep='\t', index=False)
        extract_eggnog(neigh, smorf_id, f"{out_dir}/EggNOG_annotations.tsv")

print("Done")
