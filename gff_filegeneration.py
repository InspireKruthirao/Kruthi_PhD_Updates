import pandas as pd
from pathlib import Path

BASE_DIR = "/work/microbiome/users/kruthi/SmORF_neighbourhoods_26_30"

def create_gff_for_contig(contig_dir):
    contig_path = Path(contig_dir)
    neighbors = pd.read_csv(contig_path / "neighbors.tsv", sep='\t')
    
    target_contig = None
    if 'Note' in neighbors.columns:
        target_row = neighbors[neighbors['Note'].str.contains('TARGET', na=False)]
        if not target_row.empty:
            target_contig = target_row.iloc[0]['Contig']
    
    with open(contig_path / "EggNOG_annotations.tsv", 'r') as f:
        header = f.readline().lstrip('#').strip().split('\t')
    eggnog = pd.read_csv(contig_path / "EggNOG_annotations.tsv", sep='\t', skiprows=1, names=header)
   
    annotation_map = {}
    for _, row in eggnog.iterrows():
        cog_name = "Unknown"
        if pd.notna(row['eggNOG_OGs']) and 'COG' in str(row['eggNOG_OGs']):
            cog_parts = str(row['eggNOG_OGs']).split('@')[0].split(',')[0]
            if cog_parts.startswith('COG'):
                cog_category = str(row['COG_category']).strip() if pd.notna(row['COG_category']) else ""
                cog_name = f"{cog_parts}-{cog_category}" if cog_category and cog_category != '-' else cog_parts
        annotation_map[row['Contig']] = cog_name
   
    contig_name = contig_path.name.split('_', 1)[1]
    
    with open(contig_path / f"{contig_name}.gff", 'w') as gff:
        gff.write("##gff-version 3\n")
        for _, row in neighbors.iterrows():
            name = annotation_map.get(row['Contig'], "Unknown")
            score = "1" if name != "Unknown" else "0.0"
            strand = '+' if row['Strand'] == 1 else '-'
            
            attributes = [f"ID={row['Contig']}", f"Name={name}"]
            
            if row['Contig'] == target_contig:
                attributes.extend(["target=1", "Note=TARGET_sMORF", "colour=1"])
            else:
                attributes.append("target=0")
            
            gff.write(f"{row['Contig']}\tProdigal_v2.6.3\tCDS\t{row['Start']}\t{row['End']}\t{score}\t{strand}\t0\t{'；'.join(attributes)}\n")

def process_all_smorfs():
    base_path = Path(BASE_DIR)
    for smorf_dir in sorted(base_path.glob("SHD1_SM.100AA.*")):
        print(f"Processing {smorf_dir.name}...")
        for contig_dir in sorted(smorf_dir.glob("D*_contig_*")):
            create_gff_for_contig(contig_dir)
    print("Don")

if __name__ == "__main__":
    process_all_smorfs()
