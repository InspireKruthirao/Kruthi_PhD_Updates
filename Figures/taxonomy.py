import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import re

GMSC_DIR = Path("/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER")
PLOT_OUTPUT = Path("/home/n12228516/taxonomic_level_coverage.png")

PLOT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

COLORS = sns.color_palette("Dark2")

def find_sample_directories(base_dir):
    pattern = re.compile(r'D\d{3}_PP1_PolcaCorr')
    return sorted([d for d in base_dir.iterdir() if d.is_dir() and pattern.match(d.name)])

def parse_taxonomy_string(tax_string):
    if pd.isna(tax_string) or not tax_string:
        return set()
    levels = set()
    for part in tax_string.split(';'):
        part = part.strip()
        if part.startswith('d__'): levels.add('domain')
        elif part.startswith('p__'): levels.add('phylum')
        elif part.startswith('c__'): levels.add('class')
        elif part.startswith('o__'): levels.add('order')
        elif part.startswith('f__'): levels.add('family')
        elif part.startswith('g__'): levels.add('genus')
        elif part.startswith('s__'): levels.add('species')
    return levels

def count_taxonomic_coverage(base_dir):
    sample_dirs = find_sample_directories(base_dir)
    smorfs_per_level = {lvl: set() for lvl in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']}
    total_smorfs = set()

    for sample_dir in sample_dirs:
        taxonomy_file = sample_dir / "taxonomy.out.smorfs.tsv"
        if not taxonomy_file.exists():
            continue
        df = pd.read_csv(taxonomy_file, sep='\t')
        for _, row in df.iterrows():
            smorf_id = row.get('qseqid') or row.get('q_seqid')
            taxonomy = row['taxonomy']
            total_smorfs.add(smorf_id)
            for level in parse_taxonomy_string(taxonomy):
                smorfs_per_level[level].add(smorf_id)

    counts = {lvl: len(smorfs) for lvl, smorfs in smorfs_per_level.items()}
    return counts, len(total_smorfs)

def plot_taxonomic_coverage(counts, output_path):
    level_order = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    levels = [lvl.capitalize() for lvl in level_order if counts[lvl] > 0]
    values = [counts[lvl] for lvl in level_order if counts[lvl] > 0]

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(len(levels))
    ax.barh(y_pos, values, color=[COLORS[i % len(COLORS)] for i in range(len(levels))],
            edgecolor='black', linewidth=1, alpha=0.85)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(levels, fontsize=11)
    ax.set_xlabel('Number of smORFs', fontsize=12, fontweight='bold')
    ax.set_ylabel('Taxonomic Level', fontsize=12, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(False)

    # Full border around the plot
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=400, bbox_inches='tight')
    plt.close()

def print_summary(counts, total_smorfs):
    level_order = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    print("\n" + "="*60)
    print("TAXONOMIC ANNOTATION COVERAGE SUMMARY")
    print("="*60)
    print(f"Total unique smORFs with taxonomy: {total_smorfs:,}\n")
    
    for lvl in level_order:
        if counts[lvl] > 0:
            pct = counts[lvl] / total_smorfs * 100
            print(f"{lvl.capitalize():10s}: {counts[lvl]:,} smORFs ({pct:5.1f}%)")
    
    # Paragraph for reports
    annotated_levels = [f"{lvl.capitalize()} ({counts[lvl]:,}; {counts[lvl]/total_smorfs*100:.1f}%)" 
                        for lvl in level_order if counts[lvl] > 0]
    paragraph = ("Out of {:,} unique smORFs detected across the Shanghai dog metagenomes, "
                 "taxonomic annotation coverage was as follows: " + 
                 ", ".join(annotated_levels) + ".")
    
    print("\nSuggested report paragraph:")
    print("-" * 60)
    print(paragraph.format(total_smorfs))
    print("-" * 60)

def main():
    counts, total_smorfs = count_taxonomic_coverage(GMSC_DIR)
    plot_taxonomic_coverage(counts, PLOT_OUTPUT)
    print(f"\nPlot with border saved to: {PLOT_OUTPUT}")
    print_summary(counts, total_smorfs)

if __name__ == "__main__":
    main()
