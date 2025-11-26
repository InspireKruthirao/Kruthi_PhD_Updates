import pandas as pd

path = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/SHD_SMORF_resource/100AA_SmORFs_origins.tsv.gz"

df = pd.read_csv(
    path,
    sep="\t",
    header=0,          
    usecols=[0, 1],   
    compression="gzip"
)

df.columns = ["smorf", "sample"]

# Counts
total_rows = len(df)
unique_smorfs = df["smorf"].nunique()
unique_samples = df["sample"].nunique()

prevalence = df.groupby("smorf")["sample"].nunique()

print("SUMMARY")
print("───────")
print(f"Total smORF occurrences (rows): {total_rows:,}")
print(f"Unique smORFs                  : {unique_smorfs:,}")
print(f"Unique dog samples             : {unique_samples:,}")
print(f"Most prevalent smORF in        : {prevalence.max()} dogs")

print("\nPREVALENCE THRESHOLDS")
print("─────────────────────")
print(f"Present in ≥40 dogs            : {(prevalence >= 40).sum():,}")
print(f"Present in ≥47 dogs            : {(prevalence >= 47).sum():,}")
print(f"Present in all 52 dogs         : {(prevalence == 52).sum():,}")
