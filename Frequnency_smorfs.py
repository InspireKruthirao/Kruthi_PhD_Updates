#!/usr/bin/env python3
import pandas as pd
import os

SMORF_FILE = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
GENE_CATALOG_FILE = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
CLUSTERS_FILE = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.clusters.tsv.xz"
EGGNOG_BASE_DIR = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"
OUTPUT_DIR = "/work/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods"
RANK_FILE = "/work/microbiome/shanghai_dogs/intermediate-outputs/Ranked_SmORFs_optionC.tsv"

smorf_df = pd.read_csv(SMORF_FILE, sep='\t', compression='gzip',
                       names=['SmORF_ID', 'Sample', 'Contig', 'Coordinates', 'Strand'])

occurrences = smorf_df.groupby('SmORF_ID').size().rename('total_occurrences')
samples = smorf_df.groupby('SmORF_ID')['Sample'].nunique().rename('n_samples')

def is_good_contig(contig):
    return len(contig.split("_polypolish")) == 2

from collections import defaultdict
good_bad_counts = defaultdict(lambda: {'good':0, 'bad':0})
for _, row in smorf_df.iterrows():
    if is_good_contig(row['Contig']):
        good_bad_counts[row['SmORF_ID']]['good'] += 1
    else:
        good_bad_counts[row['SmORF_ID']]['bad'] += 1

good_counts = pd.Series({k:v['good'] for k,v in good_bad_counts.items()}, name='good_contigs')
bad_counts = pd.Series({k:v['bad'] for k,v in good_bad_counts.items()}, name='bad_contigs')

rank_df = pd.concat([occurrences, good_counts, bad_counts, samples], axis=1).fillna(0)
rank_df['good_fraction'] = rank_df['good_contigs'] / rank_df['total_occurrences']
os.makedirs(os.path.dirname(RANK_FILE), exist_ok=True)
rank_df.reset_index().rename(columns={'index':'SmORF_ID'}).to_csv(RANK_FILE, sep='\t', index=False)

subset = rank_df[(rank_df['total_occurrences'] >= 20) & (rank_df['total_occurrences'] <= 30)]
selected_smorfs = subset.sample(5, random_state=42).index.tolist()

gene_df = pd.read_csv(GENE_CATALOG_FILE, sep='\t', compression='xz',
                      names=['ORF', 'Sample', 'Contig', 'Start', 'End', 'Strand', 'Partial'], dtype=str)
clusters_df = pd.read_csv(CLUSTERS_FILE, sep=r'\s+', header=None, compression='xz', low_memory=False)

def get_neighbours(sample, contig, start, end, window=5):
    contig_df = gene_df[(gene_df['Sample'] == sample) & (gene_df['Contig'].str.startswith(contig))].copy()
    contig_df['Start'] = contig_df['Start'].astype(int)
    contig_df['End'] = contig_df['End'].astype(int)
    contig_df = contig_df.sort_values('Start').reset_index(drop=True)
    target_idx = contig_df[(contig_df['Start'] == start) & (contig_df['End'] == end)].index
    if len(target_idx) == 0:
        return pd.DataFrame()
    idx = target_idx[0]
    left = max(0, idx - window)
    right = min(len(contig_df), idx + window + 1)
    neighbours = contig_df.iloc[left:right].copy()
    neighbours['Note'] = ""
    neighbours.loc[(neighbours['Start'] == start) & (neighbours['End'] == end), 'Note'] = f"TARGET ({smorf_id})"
    return neighbours

def extract_eggnog(neighbour_df, smorf_id, output_path):
    with open(output_path, "w") as out:
        out.write("#query\tORF\tSample\tContig\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
        for _, row in neighbour_df.iterrows():
            orf_id, sample_id, contig_id = row['ORF'], row['Sample'], row['Contig']
            annot_file = os.path.join(EGGNOG_BASE_DIR, sample_id, f"{sample_id}.emapper.annotations")
            if not os.path.exists(annot_file):
                continue
            with open(annot_file, "r") as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    if contig_id in line:
                        out.write(f"{smorf_id}\t{orf_id}\t{sample_id}\t{contig_id}\t{line}")
                        break

for smorf_id in selected_smorfs:
    smorf_occurrences = smorf_df[smorf_df['SmORF_ID'] == smorf_id]
    for _, occ in smorf_occurrences.iterrows():
        sample = occ['Sample']
        contig_full = occ['Contig']
        contig = contig_full.split("_polypolish")[0]
        start, end = map(int, occ['Coordinates'].split('-'))

        out_dir = os.path.join(OUTPUT_DIR, smorf_id, f"{sample}_{contig}")
        os.makedirs(out_dir, exist_ok=True)

        neighbours_df = get_neighbours(sample, contig, start, end)
        if neighbours_df.empty:
            continue

        neighbours_df.to_csv(os.path.join(out_dir, "neighbors.tsv"), sep="\t", index=False)
        extract_eggnog(neighbours_df, smorf_id, os.path.join(out_dir, "EggNOG_annotations.tsv"))

print("success")
