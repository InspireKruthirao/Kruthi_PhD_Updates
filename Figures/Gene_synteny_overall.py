#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from pathlib import Path

BASE_DIR = Path("/work/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods/SHD1_SM.100AA.006_627")
SMORF_ID = "SHD1_SM.100AA.006_627"
WINDOW = 15
ROW_HEIGHT = 0.55
OUT_PNG = "SHD1_SM.100AA.006_627_context_arrows.png"
OUT_PDF = "SHD1_SM.100AA.006_627_context_arrows.pdf"

COG_COLORS = {
    'J': '#6ec4e8', 'A': '#ff9999', 'K': '#66c2a5', 'L': '#fc8d62', 'B': '#e78ac3',
    'D': '#ffd92f', 'V': '#b3b3b3', 'T': '#8da0cb', 'M': '#e5c494', 'N': '#b3de69',
    'U': '#fccde5', 'O': '#d9d9d9', 'X': '#bc80bd', 'C': '#fb8072', 'G': '#bebada',
    'E': '#ffffb3', 'F': '#80b1d3', 'H': '#fdb462', 'I': '#b3de69', 'P': '#fccde5',
    'Q': '#d9d9d9', 'R': '#cccccc', 'S': '#999999', '-': '#ffffff'
}
UNKNOWN = "#dddddd"

def load_neigh(p): return pd.read_csv(p, sep='\t', comment='#')
def load_egg(p):
    df = pd.read_csv(p, sep='\t', comment=None)
    if df.columns[0].startswith('#'): df.columns = [c.lstrip('#') for c in df.columns]
    if '#query' in df.columns: df = df.rename(columns={'#query': 'query'})
    return df[['query', 'COG_category']].dropna(subset=['query'])

dirs = sorted([p for p in BASE_DIR.iterdir() if p.is_dir() and p.name.startswith("D")])
contexts = []

for d in dirs:
    neigh_f = d / "neighbors.tsv"
    egg_f   = d / "EggNOG_annotations.tsv"
    if not (neigh_f.exists() and egg_f.exists()): continue

    neigh = load_neigh(neigh_f)
    egg   = load_egg(egg_f)

    if not neigh['Note'].astype(str).str.contains(SMORF_ID).any(): continue

    tgt_idx = neigh['Note'].astype(str).str.contains(SMORF_ID).idxmax()
    start = max(0, tgt_idx - WINDOW)
    end   = min(len(neigh), tgt_idx + WINDOW + 1)
    win   = neigh.iloc[start:end].reset_index(drop=True)
    tgt_pos = win[win['Note'].astype(str).str.contains(SMORF_ID)].index[0]

    egg_map = dict(zip(egg['query'], egg['COG_category']))
    def get_cog(q):
        if pd.isna(q): return "-"
        c = egg_map.get(q, "-")
        return "-" if pd.isna(c) or c in ("", "-") else str(c).strip()[0]

    win['COG'] = win.get('query', win['ORF']).apply(get_cog)
    contexts.append({'name': d.name.split("_")[0], 'win': win, 'pos': tgt_pos})

contexts.sort(key=lambda x: x['name'])

n = len(contexts)
fig, ax = plt.subplots(figsize=(16, max(6, n * ROW_HEIGHT)))

for i, ctx in enumerate(reversed(contexts)):
    y = i
    for j, (cog, strand) in enumerate(zip(ctx['win']['COG'], ctx['win']['Strand'])):
        x = j - ctx['pos']
        dir_arrow = '->' if str(strand) == '1' else '<-'
        color = COG_COLORS.get(cog, UNKNOWN)

        if j == ctx['pos']:
            arrow = FancyArrowPatch((x-0.48, y), (x+0.48, y), arrowstyle='->,head_length=10,head_width=10,widthB=1.4',
                                    mutation_scale=22, lw=3.5, edgecolor='black', facecolor='#e31a1c')
            ax.text(x, y+0.18, "TARGET", ha='center', va='bottom', fontsize=8, fontweight='bold')
        else:
            arrow = FancyArrowPatch((x-0.45, y), (x+0.45, y), arrowstyle=f'{dir_arrow},head_length=7,head_width=7,widthA=0.7',
                                    mutation_scale=16, lw=0.9, edgecolor='#333333', facecolor=color)
            if cog != "-": ax.text(x, y, cog, ha='center', va='center', fontsize=9, fontweight='bold')

        ax.add_patch(arrow)

    ax.text(-WINDOW-1.5, y, ctx['name'], ha='right', va='center', fontsize=11, fontweight='bold')

ax.set_ylim(-0.8, n)
ax.set_xlim(-WINDOW-2, WINDOW+2)
ax.axis('off')
plt.suptitle(f"Conserved genomic context of {SMORF_ID} (n = {n})", fontsize=18, fontweight='bold', y=0.96)

handles = [FancyArrowPatch((0,0),(1,0), arrowstyle='->,head_length=10,head_width=10,widthB=1.4', mutation_scale=22,
                           lw=3.5, edgecolor='black', facecolor='#e31a1c', label='Target sMORF'),
           plt.Rectangle((0,0),1,1, facecolor=UNKNOWN, edgecolor='gray', label='No annotation')]
for c in 'LHGCE':
    handles.append(plt.Rectangle((0,0),1,1, facecolor=COG_COLORS[c],
                 label=f"[{c}] {'Replication/repair' if c=='L' else 'Coenzyme' if c=='H' else 'Carbohydrate' if c=='G' else 'Energy' if c=='C' else 'Amino acid'}"))

ax.legend(handles=handles, bbox_to_anchor=(1.01, 0.98), loc='upper left', fontsize=11, title="COG Category")

plt.tight_layout()
plt.subplots_adjust(right=0.82)
plt.savefig(OUT_PNG, dpi=400, bbox_inches='tight')
plt.savefig(OUT_PDF, bbox_inches='tight')
