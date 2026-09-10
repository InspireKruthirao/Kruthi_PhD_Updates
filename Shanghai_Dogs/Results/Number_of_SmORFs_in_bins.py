import gzip
from collections import Counter

INPUT = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
BIN_WIDTH = 5

def bin_label(n: int) -> str:
    start = ((n - 1) // BIN_WIDTH) * BIN_WIDTH + 1
    end = start + BIN_WIDTH - 1
    return f"{start}-{end}"

def main():
    smorf_counts = Counter()

    with gzip.open(INPUT, "rt") as f:
        next(f)
        for line in f:
            smorf_id = line.split("\t", 1)[0]
            smorf_counts[smorf_id] += 1

    max_count = max(smorf_counts.values())
    max_bin_end = ((max_count - 1) // BIN_WIDTH + 1) * BIN_WIDTH

    binned = Counter()
    for c in smorf_counts.values():
        binned[bin_label(c)] += 1

    for start in range(1, max_bin_end + 1, BIN_WIDTH):
        end = start + BIN_WIDTH - 1
        lbl = f"{start}-{end}"
        print(f"{lbl}\t{binned.get(lbl, 0)}")

if __name__ == "__main__":
    main()
