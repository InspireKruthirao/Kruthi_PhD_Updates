import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np

# ── PATHS ──────────────────────────────────────────────────────────────────────
INPUT_TSV = (
    "/work/microbiome/users/kruthi"
    "/SmORF_neighbourhoods_21_25_no_min_overlap"
    "/smorf_onehot_predictions.tsv"
)
OUTPUT_PNG = (
    "/work/microbiome/users/kruthi"
    "/SmORF_neighbourhoods_21_25_no_min_overlap"
    "/smorf_predictions_plots.png"
)
# ──────────────────────────────────────────────────────────────────────────────

COG_DESC = {
    'C': 'Energy production',        'D': 'Cell cycle',
    'E': 'Amino acid metabolism',    'F': 'Nucleotide metabolism',
    'G': 'Carbohydrate metabolism',  'H': 'Coenzyme metabolism',
    'I': 'Lipid metabolism',         'J': 'Translation/ribosome',
    'K': 'Transcription',            'L': 'DNA replication/repair',
    'M': 'Cell wall/membrane',       'N': 'Cell motility',
    'O': 'Protein turnover',         'P': 'Ion transport',
    'Q': 'Secondary metabolites',    'T': 'Signal transduction',
    'U': 'Intracellular trafficking','V': 'Defense mechanisms',
}

# Dark2 palette
DARK2 = [
    '#1B9E77', '#D95F02', '#7570B3', '#E7298A',
    '#66A61E', '#E6AB02', '#A6761D', '#666666',
]

TIER_COLORS = {
    'VERY HIGH': DARK2[0],
    'HIGH':      DARK2[2],
    'MEDIUM':    DARK2[5],
    'LOW':       DARK2[7],
}

BG        = '#1C1C1C'
PANEL_BG  = '#2A2A2A'
GRID_COL  = '#444444'
TICK_COL  = '#AAAAAA'
LABEL_COL = '#CCCCCC'
TITLE_COL = '#EEEEEE'


def style_ax(ax):
    ax.set_facecolor(PANEL_BG)
    ax.tick_params(colors=TICK_COL, labelsize=9)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID_COL)


def main():
    df = pd.read_csv(INPUT_TSV, sep='\t')
    print(f"Loaded {len(df)} smORFs")

    fig = plt.figure(figsize=(16, 14), facecolor=BG)
    gs  = GridSpec(3, 2, figure=fig,
                   hspace=0.45, wspace=0.35,
                   left=0.07, right=0.97,
                   top=0.93, bottom=0.06)

    ax_donut   = fig.add_subplot(gs[0, 0])
    ax_cogbar  = fig.add_subplot(gs[0, 1])
    ax_score   = fig.add_subplot(gs[1, 0])
    ax_scatter = fig.add_subplot(gs[1, 1])
    ax_stack   = fig.add_subplot(gs[2, :])

    for ax in [ax_donut, ax_cogbar, ax_score, ax_scatter, ax_stack]:
        style_ax(ax)

    # ── Plot 1: Confidence tier donut ──────────────────────────────────────────
    tier_order  = ['VERY HIGH', 'HIGH', 'MEDIUM', 'LOW']
    tier_counts = [df['confidence_tier'].value_counts().get(t, 0)
                   for t in tier_order]
    tier_cols   = [TIER_COLORS[t] for t in tier_order]

    wedges, _, autotexts = ax_donut.pie(
        tier_counts,
        colors=tier_cols,
        autopct='%1.1f%%',
        startangle=90,
        pctdistance=0.78,
        wedgeprops=dict(width=0.5, edgecolor=BG, linewidth=1.5)
    )
    for at in autotexts:
        at.set_color('white')
        at.set_fontsize(8)
        at.set_fontweight('bold')

    ax_donut.set_title('Confidence tier distribution',
                       color=TITLE_COL, fontsize=11, pad=10)

    patches = [
        mpatches.Patch(color=TIER_COLORS[t],
                       label=f'{t}  (n={c})')
        for t, c in zip(tier_order, tier_counts)
    ]
    ax_donut.legend(handles=patches, loc='lower center',
                    bbox_to_anchor=(0.5, -0.18), ncol=2,
                    fontsize=8, framealpha=0, labelcolor=LABEL_COL)

    # ── Plot 2: COG horizontal bar ─────────────────────────────────────────────
    top_cogs   = df['predicted_cog_cat'].value_counts().head(12)
    cog_labels = [f"{cat} — {COG_DESC.get(cat, cat)}"
                  for cat in top_cogs.index]
    cog_colors = [DARK2[i % len(DARK2)] for i in range(len(top_cogs))]

    bars = ax_cogbar.barh(
        cog_labels[::-1],
        top_cogs.values[::-1],
        color=cog_colors[::-1],
        edgecolor='none',
        height=0.65
    )
    for bar, val in zip(bars, top_cogs.values[::-1]):
        ax_cogbar.text(
            val + 3, bar.get_y() + bar.get_height() / 2,
            str(val), va='center', ha='left',
            color=LABEL_COL, fontsize=8
        )

    ax_cogbar.set_title('Predicted COG categories (top 12)',
                        color=TITLE_COL, fontsize=11, pad=10)
    ax_cogbar.set_xlabel('Number of smORFs', color=TICK_COL, fontsize=9)
    ax_cogbar.tick_params(axis='y', labelsize=8.5, labelcolor=LABEL_COL)
    ax_cogbar.tick_params(axis='x', labelcolor=TICK_COL)
    ax_cogbar.xaxis.grid(True, color=GRID_COL, linewidth=0.5, linestyle='--')
    ax_cogbar.set_axisbelow(True)

    # ── Plot 3: Vote score distribution ────────────────────────────────────────
    scores     = df['vote_score_pct'].values
    bins       = np.arange(0, 105, 5)
    counts, edges = np.histogram(scores, bins=bins)

    bar_colors = []
    for e in edges[:-1]:
        if e >= 80:
            bar_colors.append(DARK2[0])
        elif e >= 50:
            bar_colors.append(DARK2[2])
        elif e >= 30:
            bar_colors.append(DARK2[5])
        else:
            bar_colors.append(DARK2[7])

    ax_score.bar(range(len(counts)), counts,
                 color=bar_colors, edgecolor=BG,
                 linewidth=0.5, width=0.85)

    bin_labels = [f'{int(e)}-{int(e+5)}' for e in edges[:-1]]
    ax_score.set_xticks(range(0, len(counts), 2))
    ax_score.set_xticklabels(bin_labels[::2],
                              rotation=45, ha='right',
                              color=TICK_COL, fontsize=8)
    ax_score.set_ylabel('Number of smORFs', color=TICK_COL, fontsize=9)
    ax_score.set_xlabel('Vote score (%)', color=TICK_COL, fontsize=9)
    ax_score.set_title('Vote score distribution',
                       color=TITLE_COL, fontsize=11, pad=10)
    ax_score.yaxis.grid(True, color=GRID_COL, linewidth=0.5, linestyle='--')
    ax_score.set_axisbelow(True)
    ax_score.tick_params(axis='y', labelcolor=TICK_COL)

    # ── Plot 4: Scatter conservation vs vote score ─────────────────────────────
    for tier in ['LOW', 'MEDIUM', 'HIGH', 'VERY HIGH']:
        sub = df[df['confidence_tier'] == tier]
        ax_scatter.scatter(
            sub['conservation_pct'],
            sub['vote_score_pct'],
            c=TIER_COLORS[tier],
            alpha=0.55,
            s=18,
            edgecolors='none',
            label=tier,
            zorder=3
        )

    ax_scatter.set_xlabel('Neighborhood conservation (%)',
                           color=TICK_COL, fontsize=9)
    ax_scatter.set_ylabel('Vote score (%)', color=TICK_COL, fontsize=9)
    ax_scatter.set_title('Conservation vs vote score',
                          color=TITLE_COL, fontsize=11, pad=10)
    ax_scatter.legend(fontsize=8, framealpha=0.2,
                      facecolor='#333333', labelcolor=LABEL_COL,
                      loc='lower right', edgecolor='#555555')
    ax_scatter.xaxis.grid(True, color=GRID_COL, linewidth=0.5, linestyle='--')
    ax_scatter.yaxis.grid(True, color=GRID_COL, linewidth=0.5, linestyle='--')
    ax_scatter.set_axisbelow(True)
    ax_scatter.tick_params(axis='both', labelcolor=TICK_COL)

    # ── Plot 5: Stacked bar — unknown smORFs ───────────────────────────────────
    unk   = df[df['current_annotation'] == 'Unknown']
    n_unk = len(unk)
    ucats = unk['predicted_cog_cat'].value_counts()

    top8    = ucats.head(8).index.tolist()
    other_n = ucats.iloc[8:].sum() if len(ucats) > 8 else 0
    vals    = [ucats.get(c, 0) for c in top8]

    if other_n > 0:
        top8.append('Other')
        vals.append(other_n)

    cat_labels = [f'{c} — {COG_DESC.get(c, "Other")}' for c in top8]
    stk_colors = DARK2[:len(top8)]

    left = 0
    for cat, label, val, col in zip(top8, cat_labels, vals, stk_colors):
        ax_stack.barh(
            [f'Unknown smORFs (n={n_unk})'],
            val, left=left,
            color=col, edgecolor=BG,
            linewidth=0.8, height=0.45,
            label=label
        )
        if val > 20:
            ax_stack.text(
                left + val / 2, 0,
                f'{cat}\n{val}',
                ha='center', va='center',
                color='white', fontsize=8, fontweight='bold'
            )
        left += val

    ax_stack.set_xlim(0, n_unk + 20)
    ax_stack.set_xlabel('Number of smORFs', color=TICK_COL, fontsize=9)
    ax_stack.set_title(
        f'Predicted functions — previously unknown smORFs (n={n_unk})',
        color=TITLE_COL, fontsize=11, pad=10
    )
    ax_stack.tick_params(axis='x', labelcolor=TICK_COL)
    ax_stack.tick_params(axis='y', labelcolor=LABEL_COL, labelsize=9)
    ax_stack.xaxis.grid(True, color=GRID_COL, linewidth=0.5, linestyle='--')
    ax_stack.set_axisbelow(True)

    ax_stack.legend(
        loc='lower center', bbox_to_anchor=(0.5, -0.52),
        ncol=5, fontsize=8, framealpha=0.2,
        facecolor='#333333', edgecolor='#555555',
        labelcolor=LABEL_COL
    )

    # ── Title and save ─────────────────────────────────────────────────────────
    fig.suptitle(
        'smORF function prediction — one-hot weighted vote (bin 21-25)',
        color=TITLE_COL, fontsize=13, fontweight='bold', y=0.97
    )

    plt.savefig(OUTPUT_PNG, dpi=150,
                bbox_inches='tight', facecolor=BG)
    print(f"Saved → {OUTPUT_PNG}")


if __name__ == '__main__':
    main()
