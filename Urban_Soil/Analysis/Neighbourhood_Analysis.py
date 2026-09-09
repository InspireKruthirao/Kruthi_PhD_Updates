cat > /work/microbiome/users/kruthi/intermediate_results/Kruthi_PhD_Updates/Urban_Soil/Neighbourhood_Analysis/Prepare_neighbourhood_inputs.py <<'PY'
#!/usr/bin/env python3

import csv
import gzip
import lzma
from collections import Counter
from pathlib import Path

SMORF_FILE = Path(
    "/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil/"
    "UrbanSoil_SMORF_resource/UrbanSoil_100AA_SMORFs_origins.tsv.gz"
)

GENE_CATALOG = Path(
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil/"
    "Gene_Catalog_Urban_Soil/UrbanSoil.ORF.orig.tsv.xz"
)

EGGNOG_DIR = Path(
    "/work/microbiome/users/kruthi/intermediate_results/egg_nog/results"
)

OUT_ROOT = Path(
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil/"
    "Neighbourhood_Analysis"
)

FINAL_DIR = OUT_ROOT / "prepared_inputs"
PART_DIR = OUT_ROOT / "prepared_inputs.part"

BIN_WIDTH = 5
PROGRESS_EVERY = 5_000_000


def get_samples():
    samples = []

    for sample_dir in sorted(EGGNOG_DIR.iterdir()):
        if not sample_dir.is_dir():
            continue

        sample = sample_dir.name

        plain_file = sample_dir / f"{sample}.emapper.annotations"
        gzip_file = sample_dir / f"{sample}.emapper.annotations.gz"

        if plain_file.exists() or gzip_file.exists():
            samples.append(sample)

    return samples


def require_columns(header, required):
    missing = [column for column in required if column not in header]

    if missing:
        raise ValueError(
            "Missing required columns: " + ", ".join(missing)
        )

    return {
        column: header.index(column)
        for column in required
    }


def open_sample_files(directory, samples, header):
    directory.mkdir(parents=True, exist_ok=False)

    handles = {}

    for sample in samples:
        path = directory / f"{sample}.tsv.gz"

        handle = gzip.open(
            path,
            "wt",
            compresslevel=1,
            newline=""
        )

        handle.write(header.rstrip("\r\n") + "\n")
        handles[sample] = handle

    return handles


def close_files(handles):
    for handle in handles.values():
        handle.close()


def main():

    for input_path in (SMORF_FILE, GENE_CATALOG, EGGNOG_DIR):
        if not input_path.exists():
            raise FileNotFoundError(input_path)

    samples = get_samples()

    if len(samples) != 58:
        raise RuntimeError(
            f"Expected 58 samples, but found {len(samples)}"
        )

    sample_set = set(samples)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    if FINAL_DIR.exists():
        raise FileExistsError(
            f"{FINAL_DIR} already exists. Nothing was overwritten."
        )

    if PART_DIR.exists():
        raise File FileExistsError(
            f"{PART_DIR} already exists. "
            "This may be from an incomplete previous run."
        )

    PART_DIR.mkdir()

    # ---------------------------------------------------------
    # PASS 1: count global occurrences of each smORF
    # ---------------------------------------------------------

    print("PASS 1/3: Counting smORF occurrences...", flush=True)

    smorf_counts = Counter()
    total_origin_rows = 0

    with gzip.open(SMORF_FILE, "rt", newline="") as source:

        reader = csv.reader(source, delimiter="\t")
        header = next(reader)

        indexes = require_columns(
            header,
            ["SmORF ID", "Sample ID"]
        )

        smorf_index = indexes["SmORF ID"]

        for row in reader:

            if not row:
                continue

            smorf_counts[row[smorf_index]] += 1
            total_origin_rows += 1

            if total_origin_rows % PROGRESS_EVERY == 0:
                print(
                    f"Counted {total_origin_rows:,} origin rows",
                    flush=True
                )

    if not smorf_counts:
        raise RuntimeError("No smORF records were found")

    counts_file = PART_DIR / "smorf_counts.tsv.gz"

    with gzip.open(
        counts_file,
        "wt",
        compresslevel=1,
        newline=""
    ) as output:

        output.write(
            "SmORF_ID\tTotalOccurrences\t"
            "BinStart\tBinEnd\tBinLabel\n"
        )

        for smorf_id, count in smorf_counts.items():

            bin_start = (
                ((count - 1) // BIN_WIDTH) * BIN_WIDTH + 1
            )

            bin_end = bin_start + BIN_WIDTH - 1

            output.write(
                f"{smorf_id}\t{count}\t"
                f"{bin_start}\t{bin_end}\t"
                f"{bin_start}_{bin_end}\n"
            )

    print(
        f"Unique smORFs: {len(smorf_counts):,}",
        flush=True
    )

    # ---------------------------------------------------------
    # PASS 2: split smORF origins by sample
    # ---------------------------------------------------------

    print(
        "PASS 2/3: Splitting smORF origins by sample...",
        flush=True
    )

    origins_directory = PART_DIR / "origins_by_sample"

    origins_header = (
        "SmORF ID\tSample ID\tContig\tCoordinates\tStrand\t"
        "TotalOccurrences\tBinStart\tBinEnd\tBinLabel"
    )

    origin_outputs = open_sample_files(
        origins_directory,
        samples,
        origins_header
    )

    origin_counts = Counter()
    written_origins = 0
    unexpected_samples = Counter()

    try:

        with gzip.open(SMORF_FILE, "rt", newline="") as source:

            reader = csv.reader(source, delimiter="\t")
            header = next(reader)

            indexes = require_columns(
                header,
                ["SmORF ID", "Sample ID"]
            )

            smorf_index = indexes["SmORF ID"]
            sample_index = indexes["Sample ID"]

            for row in reader:

                if not row:
                    continue

                smorf_id = row[smorf_index]
                sample = row[sample_index]

                if sample not in sample_set:
                    unexpected_samples[sample] += 1
                    continue

                count = smorf_counts[smorf_id]

                bin_start = (
                    ((count - 1) // BIN_WIDTH) * BIN_WIDTH + 1
                )

                bin_end = bin_start + BIN_WIDTH - 1

                origin_outputs[sample].write(
                    "\t".join(row)
                    + f"\t{count}"
                    + f"\t{bin_start}"
                    + f"\t{bin_end}"
                    + f"\t{bin_start}_{bin_end}\n"
                )

                origin_counts[sample] += 1
                written_origins += 1

                if written_origins % PROGRESS_EVERY == 0:
                    print(
                        f"Wrote {written_origins:,} origin rows",
                        flush=True
                    )

    finally:
        close_files(origin_outputs)

    if unexpected_samples:
        raise RuntimeError(
            "Unexpected origin samples: "
            + ", ".join(sorted(unexpected_samples))
        )

    if written_origins != total_origin_rows:
        raise RuntimeError(
            f"Origin row mismatch: counted {total_origin_rows:,}, "
            f"wrote {written_origins:,}"
        )

    # ---------------------------------------------------------
    # PASS 3: split gene catalogue by sample
    # ---------------------------------------------------------

    print(
        "PASS 3/3: Splitting gene catalogue by sample...",
        flush=True
    )

    genes_directory = PART_DIR / "genes_by_sample"

    gene_outputs = None
    gene_counts = Counter()
    unexpected_gene_samples = Counter()
    total_gene_rows = 0

    try:

        with lzma.open(
            GENE_CATALOG,
            "rt",
            newline=""
        ) as source:

            gene_header_line = source.readline()

            if not gene_header_line:
                raise RuntimeError("Gene catalogue is empty")

            gene_header = (
                gene_header_line.rstrip("\r\n").split("\t")
            )

            indexes = require_columns(
                gene_header,
                [
                    "ORF",
                    "Sample",
                    "Original_ID",
                    "Start",
                    "End",
                    "Strand",
                    "Partial"
                ]
            )

            sample_index = indexes["Sample"]

            gene_outputs = open_sample_files(
                genes_directory,
                samples,
                gene_header_line
            )

            for line in source:

                if not line.strip():
                    continue

                fields = line.rstrip("\r\n").split("\t")

                if len(fields) <= sample_index:
                    raise ValueError(
                        "Malformed gene-catalogue row near "
                        f"row {total_gene_rows + 2:,}"
                    )

                sample = fields[sample_index]

                if sample not in sample_set:
                    unexpected_gene_samples[sample] += 1
                    continue

                gene_outputs[sample].write(
                    "\t".join(fields) + "\n"
                )

                gene_counts[sample] += 1
                total_gene_rows += 1

                if total_gene_rows % PROGRESS_EVERY == 0:
                    print(
                        f"Wrote {total_gene_rows:,} gene rows",
                        flush=True
                    )

    finally:

        if gene_outputs is not None:
            close_files(gene_outputs)

    if unexpected_gene_samples:
        raise RuntimeError(
            "Unexpected gene samples: "
            + ", ".join(sorted(unexpected_gene_samples))
        )

    missing_origins = [
        sample
        for sample in samples
        if origin_counts[sample] == 0
    ]

    missing_genes = [
        sample
        for sample in samples
        if gene_counts[sample] == 0
    ]

    if missing_origins or missing_genes:
        raise RuntimeError(
            f"Samples without origins: {missing_origins}; "
            f"samples without genes: {missing_genes}"
        )

    samples_file = PART_DIR / "samples.tsv"

    with samples_file.open("w") as output:

        output.write(
            "ArrayIndex\tSample\tOriginRows\tGeneRows\n"
        )

        for array_index, sample in enumerate(
            samples,
            start=1
        ):

            output.write(
                f"{array_index}\t{sample}\t"
                f"{origin_counts[sample]}\t"
                f"{gene_counts[sample]}\n"
            )

    summary_file = PART_DIR / "prepare_summary.txt"

    with summary_file.open("w") as output:

        output.write(f"Samples: {len(samples)}\n")
        output.write(
            f"Unique smORFs: {len(smorf_counts)}\n"
        )
        output.write(
            f"Origin rows: {total_origin_rows}\n"
        )
        output.write(
            f"Gene rows: {total_gene_rows}\n"
        )
        output.write(f"Bin width: {BIN_WIDTH}\n")

    del smorf_counts

    PART_DIR.rename(FINAL_DIR)

    print("PREPARATION COMPLETE", flush=True)
    print(f"Prepared inputs: {FINAL_DIR}", flush=True)
    print(
        f"Origin rows: {total_origin_rows:,}",
        flush=True
    )
    print(
        f"Gene rows: {total_gene_rows:,}",
        flush=True
    )


if __name__ == "__main__":
    main()
PY
