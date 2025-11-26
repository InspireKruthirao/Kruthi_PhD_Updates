import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")
sns.set_palette("Dark2")
plt.rcParams['figure.figsize'] = (16, 8)

df = pd.read_csv("/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/SHD_SMORF_resource/100AA_SmORFs_origins.tsv.gz",
                 sep="\t", names=["smorf", "sample"], usecols=[0,1])

total_per_sample = df['sample'].value_counts()

smorf_prevalence = df.groupby('smorf')['sample'].nunique()
unique_smORFs = smorf_prevalence[smorf_prevalence == 1].index
unique_per_sample = df[df['smorf'].isin(unique_smORFs)]['sample'].value_counts()
unique_per_sample = unique_per_sample.reindex(total_per_sample.index, fill_value=0)
shared_per_sample = total_per_sample - unique_per_sample

order = total_per_sample.sort_values(ascending=False).index

plot_df = pd.DataFrame({
    'Shared smORFs' : shared_per_sample[order],
    'Unique smORFs' : unique_per_sample[order]
})

# Plot
ax = plot_df.plot(kind='bar', stacked=True, width=0.8,
                  color=[sns.color_palette("Dark2")[1], sns.color_palette("Dark2")[2]],
                  edgecolor='white', linewidth=0.6)

plt.xlabel("Sample ID", fontsize=14, labelpad=12)
plt.ylabel("Number of predicted ≤100 aa smORFs", fontsize=14)
plt.title("Total predicted smORFs per sample", fontsize=18, fontweight='bold', pad=20)

plt.legend(fontsize=14, frameon=True, fancybox=False, edgecolor='black')

plt.xticks(rotation=90, fontsize=9)
plt.tick_params(axis='x', length=4)

plt.tight_layout()
plt.savefig("Fig_smORFs_per_sample_NO_GRID_FINAL.png", dpi=600, bbox_inches='tight')
plt.show()
print(f"Average total detections per sample : {total_per_sample.mean():,.0f}")
print(f"Average sample-specific smORFs       : {unique_per_sample.mean():,.0f}")
print(f"Percentage sample-specific           : {100*unique_per_sample.mean()/total_per_sample.mean():.1f}%")
