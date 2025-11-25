#!/usr/bin/env python3
"""
Generate GFF files for each contig from neighbors.tsv and EggNOG_annotations.tsv
"""

import pandas as pd
import os
from pathlib import Path

def create_gff_for_contig(contig_dir):
    """
    Create a GFF file for a single contig directory
    
    Args:
        contig_dir: Path to the contig directory (e.g., 'D000_contig_614')
    """
    contig_path = Path(contig_dir)
    
    # Read the files
    neighbors_file = contig_path / "neighbors.tsv"
    eggnog_file = contig_path / "EggNOG_annotations.tsv"
    
    if not neighbors_file.exists() or not eggnog_file.exists():
        print(f"  Skipping: Missing required files")
        return
    
    # Read neighbors data
    neighbors = pd.read_csv(neighbors_file, sep='\t')

    with open(eggnog_file, 'r') as f:
        first_line = f.readline()
        # Remove the leading # and any whitespace
        header = first_line.lstrip('#').strip().split('\t')
    
    eggnog = pd.read_csv(eggnog_file, sep='\t', skiprows=1, names=header)
    
    print(f"  EggNOG columns: {list(eggnog.columns)[:5]}...")
    print(f"  Found {len(eggnog)} EggNOG entries")
    
    annotation_map = {}
    
    for _, row in eggnog.iterrows():
        # The contig ID is in the 'Contig' column (4th column, index 3)
        contig_id = row['Contig']
        
        cog_name = "Unknown"
        
        if 'eggNOG_OGs' in row and pd.notna(row['eggNOG_OGs']) and 'COG' in str(row['eggNOG_OGs']):
            # Extract first COG identifier (e.g., COG2878)
            cog_parts = str(row['eggNOG_OGs']).split('@')[0].split(',')[0]
            if cog_parts.startswith('COG'):
                cog_name = cog_parts
        
        annotation_map[contig_id] = cog_name
    
    print(f"  Mapped {len(annotation_map)} annotations")
    
    contig_name = contig_path.name.split('_', 1)[1]  # Gets 'contig_614'
    
    output_file = contig_path / f"{contig_name}.gff"
    
    gff_lines = 0
    with open(output_file, 'w') as gff:
        # Write GFF header
        gff.write("##gff-version 3\n")
        
        for _, row in neighbors.iterrows():
            contig_id = row['Contig']
            start = row['Start']
            end = row['End']
            strand = '+' if row['Strand'] == 1 else '-'
            
            name = annotation_map.get(contig_id, "Unknown")
            
            score = "1" if name != "Unknown" else "0.0"
            
            gff_line = (
                f"{contig_id}\tProdigal_v2.6.3\tCDS\t{start}\t{end}\t"
                f"{score}\t{strand}\t0\tID={contig_id};Name={name}\n"
            )
            
            gff.write(gff_line)
            gff_lines += 1
    
    print(f"  Created GFF file: {output_file.name} ({gff_lines} genes)")


def process_all_contigs(parent_dir="."):
    """
    Process all contig directories in the parent directory
    
    Args:
        parent_dir: Directory containing all D*_contig_* subdirectories
    """
    parent_path = Path(parent_dir)
    
    contig_dirs = sorted(parent_path.glob("D*_contig_*"))
    
    if not contig_dirs:
        print(f"No contig directories found in {parent_dir}")
        return
    
    print(f"Found {len(contig_dirs)} contig directories")
    print("="*60)
    
    success_count = 0
    error_count = 0
    
    for contig_dir in contig_dirs:
        if contig_dir.is_dir():
            print(f"\nProcessing: {contig_dir.name}")
            try:
                create_gff_for_contig(contig_dir)
                success_count += 1
            except Exception as e:
                print(f"  ERROR: {e}")
                error_count += 1
                continue
    
    print("\n" + "="*60)
    print(f"Summary: {success_count} successful, {error_count} errors")
    print("All GFF files created successfully!")


if __name__ == "__main__":
    process_all_contigs()
