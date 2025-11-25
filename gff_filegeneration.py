#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

BASE_DIR = "/work/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods"

def create_gff_for_contig(contig_dir):
    contig_path = Path(contig_dir)
    neighbors_file = contig_path / "neighbors.tsv"
    eggnog_file = contig_path / "EggNOG_annotations.tsv"
    
    if not neighbors_file.exists() or not eggnog_file.exists():
        return
    
    neighbors = pd.read_csv(neighbors_file, sep='\t')
    
    with open(eggnog_file, 'r') as f:
        header = f.readline().lstrip('#').strip().split('\t')
    
    eggnog = pd.read_csv(eggnog_file, sep='\t', skiprows=1, names=header)
    
    annotation_map = {}
    for _, row in eggnog.iterrows():
        cog_name = "Unknown"
        if pd.notna(row['eggNOG_OGs']) and 'COG' in str(row['eggNOG_OGs']):
            cog_parts = str(row['eggNOG_OGs']).split('@')[0].split(',')[0]
            if cog_parts.startswith('COG'):
                cog_name = cog_parts
        annotation_map[row['Contig']] = cog_name
    
    contig_name = contig_path.name.split('_', 1)[1]
    output_file = contig_path / f"{contig_name}.gff"
    
    with open(output_file, 'w') as gff:
        gff.write("##gff-version 3\n")
        for _, row in neighbors.iterrows():
            name = annotation_map.get(row['Contig'], "Unknown")
            score = "1" if name != "Unknown" else "0.0"
            strand = '+' if row['Strand'] == 1 else '-'
            gff.write(f"{row['Contig']}\tProdigal_v2.6.3\tCDS\t{row['Start']}\t{row['End']}\t{score}\t{strand}\t0\tID={row['Contig']};Name={name}\n")

def process_all_smorfs():
    base_path = Path(BASE_DIR)
    smorf_dirs = sorted(base_path.glob("SHD1_SM.100AA.*"))
    
    for smorf_dir in smorf_dirs:
        if smorf_dir.is_dir():
            for contig_dir in sorted(smorf_dir.glob("D*_contig_*")):
                if contig_dir.is_dir():
                    create_gff_for_contig(contig_dir)

if __name__ == "__main__":
    process_all_smorfs()
