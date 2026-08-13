import gzip
from collections import Counter

INPUT = "/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil/UrbanSoil_SMORF_resource/UrbanSoil_100AA_SMORFs_origins.tsv.gz"
BIN_WIDTH = 5


def bin_label(n: int) -> str:
    start = ((n - 1) // BIN_WIDTH) * BIN_WIDTH + 1
    end = start + BIN_WIDTH - 1
    return f"{start}-{end}"


def main():
    smorf_counts = Counter()
    total_occurrences = 0

    # Count occurrences of each smORF
    with gzip.open(INPUT, "rt") as f:
        next(f)  # skip header

        for line in f:
            smorf_id = line.split("\t", 1)[0]

            smorf_counts[smorf_id] += 1
            total_occurrences += 1

    # Summary statistics
    unique_smorfs = len(smorf_counts)

    singleton_smorfs = sum(
        1 for count in smorf_counts.values()
        if count == 1
    )

    repeated_smorfs = sum(
        1 for count in smorf_counts.values()
        if count > 1
    )

    max_count = max(smorf_counts.values())

    print("\n===== Urban Soil smORF Summary =====")

    print(f"Total smORF occurrences (redundant): {total_occurrences:,}")
    print(f"Unique smORFs (non-redundant):       {unique_smorfs:,}")
    print(f"Singleton smORFs (occur once):       {singleton_smorfs:,}")
    print(f"Repeated smORFs (occur >1 time):     {repeated_smorfs:,}")
    print(f"Maximum occurrences of one smORF:    {max_count:,}")

    # Bin occurrence frequencies
    max_bin_end = (
        ((max_count - 1) // BIN_WIDTH + 1)
        * BIN_WIDTH
    )

    binned = Counter()

    for count in smorf_counts.values():
        binned[bin_label(count)] += 1

    print("\n===== smORF Occurrence Distribution =====")

    print(
        f"{'Occurrence_Range':<18}"
        f"{'Number_of_Unique_smORFs'}"
    )

    for start in range(1, max_bin_end + 1, BIN_WIDTH):
        end = start + BIN_WIDTH - 1
        label = f"{start}-{end}"

        print(
            f"{label:<18}"
            f"{binned.get(label, 0)}"
        )


if __name__ == "__main__":
    main()
