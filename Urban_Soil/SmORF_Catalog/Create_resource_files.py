#!/usr/bin/env python3

import csv
import gzip
import importlib.util
import os
from collections import Counter, defaultdict


FASTA_PATH = (
    "/work/microbiome/shanghai_dogs/resource_generation/fasta.py"
)

LIB_PATH = (
    "/work/microbiome/shanghai_dogs/resource_generation/lib.py"
)

BASE_DIR = (
    "/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil"
)

RESOURCE_DIR = os.path.join(
    BASE_DIR,
    "UrbanSoil_SMORF_resource"
)

SAMPLE_SUFFIX = "_medaka_polypolish"
SMORF_PREFIX = "US_SM.100AA"

OUTPUT_FASTA = os.path.join(
    RESOURCE_DIR,
    "UrbanSoil_100AA_SMORFs_sequences.faa.gz"
)

ORIGINS_TSV = os.path.join(
    RESOURCE_DIR,
    "UrbanSoil_100AA_SMORFs_origins.tsv.gz"
)

HABITAT_TAXONOMY_TSV = os.path.join(
    RESOURCE_DIR,
    "UrbanSoil_100AA_SMORFs_habitat_taxonomy.tsv"
)

LOG_FILE = os.path.join(
    RESOURCE_DIR,
    "merge_log.txt"
)


def load_module(module_name, module_path):
    spec = importlib.util.spec_from_file_location(
        module_name,
        module_path
    )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def clean_sequence(sequence):
    return (
        sequence
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .upper()
    )


def get_sample_directories():
    sample_directories = []

    for sample in sorted(os.listdir(BASE_DIR)):
        sample_path = os.path.join(BASE_DIR, sample)

        if not os.path.isdir(sample_path):
            continue

        if not sample.endswith(SAMPLE_SUFFIX):
            continue

        sample_directories.append((sample, sample_path))

    return sample_directories


def check_required_files(sample_directories):
    missing_files = []

    for sample, sample_path in sample_directories:
        mapped_file = os.path.join(
            sample_path,
            "mapped.smorfs.faa"
        )

        predicted_file = os.path.join(
            sample_path,
            "predicted.filtered.smorf.faa"
        )

        if not os.path.isfile(mapped_file):
            missing_files.append(mapped_file)

        if not os.path.isfile(predicted_file):
            missing_files.append(predicted_file)

    if missing_files:
        print("Required input files are missing:")

        for path in missing_files:
            print(path)

        raise FileNotFoundError(
            f"{len(missing_files)} required files are missing"
        )


def remove_old_part_files(paths):
    for path in paths:
        if os.path.exists(path):
            os.remove(path)


def main():
    fasta = load_module("fasta", FASTA_PATH)
    lib = load_module("lib", LIB_PATH)

    print("Using:")
    print(fasta.__file__)
    print(lib.__file__)
    print("First ID:", lib.pad9(SMORF_PREFIX, 0))

    os.makedirs(RESOURCE_DIR, exist_ok=True)

    sample_directories = get_sample_directories()

    if not sample_directories:
        raise RuntimeError(
            "No Urban Soil sample directories were found"
        )

    print(
        f"Urban Soil sample directories: "
        f"{len(sample_directories)}"
    )

    check_required_files(sample_directories)

    output_fasta_part = OUTPUT_FASTA + ".part"
    origins_tsv_part = ORIGINS_TSV + ".part"
    habitat_taxonomy_part = HABITAT_TAXONOMY_TSV + ".part"
    log_file_part = LOG_FILE + ".part"

    part_files = [
        output_fasta_part,
        origins_tsv_part,
        habitat_taxonomy_part,
        log_file_part,
    ]

    remove_old_part_files(part_files)

    print("Reading mapped smORF sequences...")

    mapped_sequence_counts = Counter()
    mapped_records = 0

    for sample_number, (sample, sample_path) in enumerate(
        sample_directories,
        start=1
    ):
        mapped_file = os.path.join(
            sample_path,
            "mapped.smorfs.faa"
        )

        for header, sequence in fasta.fasta_iter(
            mapped_file,
            full_header=True
        ):
            clean_seq = clean_sequence(sequence)

            if not clean_seq:
                continue

            mapped_sequence_counts[clean_seq] += 1
            mapped_records += 1

        print(
            f"[{sample_number}/{len(sample_directories)}] "
            f"Read mapped sequences from {sample}"
        )

    if not mapped_sequence_counts:
        raise RuntimeError(
            "No mapped smORF sequences were found"
        )

    print(f"Mapped sequence records: {mapped_records:,}")
    print(
        f"Unique mapped sequences: "
        f"{len(mapped_sequence_counts):,}"
    )

    print("Sorting unique mapped sequences...")

    sorted_sequences = sorted(
        mapped_sequence_counts,
        key=lambda sequence: (
            -mapped_sequence_counts[sequence],
            sequence
        )
    )

    sequence_to_id = {}

    print("Writing pad9 smORF FASTA...")

    with gzip.open(output_fasta_part, "wt") as output_handle:
        for index, sequence in enumerate(sorted_sequences):
            smorf_id = lib.pad9(SMORF_PREFIX, index)

            output_handle.write(
                f">{smorf_id}\n{sequence}\n"
            )

            sequence_to_id[sequence] = smorf_id

    first_id = lib.pad9(SMORF_PREFIX, 0)

    last_id = lib.pad9(
        SMORF_PREFIX,
        len(sequence_to_id) - 1
    )

    print(f"First smORF ID: {first_id}")
    print(f"Last smORF ID:  {last_id}")

    del sorted_sequences
    del mapped_sequence_counts

    sequence_metadata = defaultdict(list)
    sequence_counts = Counter()

    smorf_to_habitat = defaultdict(set)
    smorf_to_taxonomy = {}

    predicted_records = 0
    matched_records = 0
    unmatched_records = 0
    invalid_headers = 0
    missing_habitat_files = 0
    missing_taxonomy_files = 0

    print("Reading predicted smORF metadata...")

    for sample_number, (sample, sample_path) in enumerate(
        sample_directories,
        start=1
    ):
        sample_id = sample.removesuffix(SAMPLE_SUFFIX)

        predicted_file = os.path.join(
            sample_path,
            "predicted.filtered.smorf.faa"
        )

        for header, sequence in fasta.fasta_iter(
            predicted_file,
            full_header=True
        ):
            predicted_records += 1
            clean_seq = clean_sequence(sequence)

            if clean_seq not in sequence_to_id:
                unmatched_records += 1
                continue

            parts = header.split("#")

            if len(parts) < 5:
                invalid_headers += 1
                continue

            smorf_id_raw = parts[0].strip()
            contig = parts[1].strip()
            start = parts[2].strip()
            end = parts[3].strip()
            strand_value = parts[4].strip()

            if strand_value == "1":
                strand = "+"
            elif strand_value == "-1":
                strand = "-"
            else:
                invalid_headers += 1
                continue

            coordinates = f"{start}-{end}"

            sequence_counts[clean_seq] += 1
            matched_records += 1

            sequence_metadata[clean_seq].append(
                (
                    sample_id,
                    contig,
                    coordinates,
                    strand,
                    smorf_id_raw,
                )
            )

        habitat_file = os.path.join(
            sample_path,
            "habitat.out.smorfs.tsv"
        )

        if os.path.isfile(habitat_file):
            with open(habitat_file, newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                next(reader, None)

                for row in reader:
                    if len(row) < 2:
                        continue

                    smorf_id_raw = row[0].strip()

                    habitats = [
                        habitat.strip()
                        for habitat in row[1].split(",")
                        if habitat.strip()
                    ]

                    smorf_to_habitat[
                        (smorf_id_raw, sample_id)
                    ].update(habitats)
        else:
            missing_habitat_files += 1

        taxonomy_file = os.path.join(
            sample_path,
            "taxonomy.out.smorfs.tsv"
        )

        if os.path.isfile(taxonomy_file):
            with open(taxonomy_file, newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                next(reader, None)

                for row in reader:
                    if len(row) < 2:
                        continue

                    smorf_id_raw = row[0].strip()
                    taxonomy = row[1].strip()

                    smorf_to_taxonomy[
                        (smorf_id_raw, sample_id)
                    ] = taxonomy
        else:
            missing_taxonomy_files += 1

        print(
            f"[{sample_number}/{len(sample_directories)}] "
            f"Processed metadata from {sample}"
        )

    sequences_with_origins = len(sequence_metadata)

    missing_origin_sequences = (
        len(sequence_to_id) - sequences_with_origins
    )

    print(f"Predicted records read: {predicted_records:,}")
    print(f"Matched origin records: {matched_records:,}")
    print(f"Unmatched records:      {unmatched_records:,}")
    print(f"Invalid headers:        {invalid_headers:,}")
    print(
        f"Sequences with origins: "
        f"{sequences_with_origins:,}"
    )
    print(
        f"Sequences without origins: "
        f"{missing_origin_sequences:,}"
    )

    if missing_origin_sequences != 0:
        raise RuntimeError(
            f"{missing_origin_sequences:,} mapped sequences "
            "have no valid origin metadata. "
            "The final files were not replaced."
        )

    smorf_order = sorted(
        sequence_to_id,
        key=lambda sequence: (
            -sequence_counts[sequence],
            sequence
        )
    )

    print("Writing origins file...")

    origin_rows_written = 0

    with gzip.open(
        origins_tsv_part,
        "wt",
        newline=""
    ) as output_handle:
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(
            [
                "SmORF ID",
                "Sample ID",
                "Contig",
                "Coordinates",
                "Strand",
            ]
        )

        for sequence in smorf_order:
            smorf_id = sequence_to_id[sequence]

            for entry in sequence_metadata[sequence]:
                (
                    sample_id,
                    contig,
                    coordinates,
                    strand,
                    smorf_id_raw,
                ) = entry

                clean_contig = contig.rsplit("_", maxsplit=1)[0]

                writer.writerow(
                    [
                        smorf_id,
                        sample_id,
                        clean_contig,
                        coordinates,
                        strand,
                    ]
                )

                origin_rows_written += 1

    print("Writing habitat and taxonomy file...")

    habitat_rows_written = 0

    with open(
        habitat_taxonomy_part,
        "w",
        newline=""
    ) as output_handle:
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(
            [
                "SmORF ID",
                "Sample ID",
                "Habitat",
                "Taxonomy",
            ]
        )

        for sequence in smorf_order:
            smorf_id = sequence_to_id[sequence]

            for entry in sequence_metadata[sequence]:
                (
                    sample_id,
                    contig,
                    coordinates,
                    strand,
                    smorf_id_raw,
                ) = entry

                habitat = ",".join(
                    sorted(
                        smorf_to_habitat.get(
                            (smorf_id_raw, sample_id),
                            set()
                        )
                    )
                )

                taxonomy = smorf_to_taxonomy.get(
                    (smorf_id_raw, sample_id),
                    ""
                )

                writer.writerow(
                    [
                        smorf_id,
                        sample_id,
                        habitat,
                        taxonomy,
                    ]
                )

                habitat_rows_written += 1

    print("Writing log file...")

    duplicate_count = sum(
        1
        for sequence in sequence_to_id
        if sequence_counts[sequence] > 1
    )

    with open(log_file_part, "w") as log_handle:
        log_handle.write(
            f"Samples processed: "
            f"{len(sample_directories)}\n"
        )

        log_handle.write(
            f"Total mapped sequence records: "
            f"{mapped_records}\n"
        )

        log_handle.write(
            f"Total unique sequences: "
            f"{len(sequence_to_id)}\n"
        )

        log_handle.write(
            f"Predicted records read: "
            f"{predicted_records}\n"
        )

        log_handle.write(
            f"Matched origin records: "
            f"{matched_records}\n"
        )

        log_handle.write(
            f"Unmatched predicted records: "
            f"{unmatched_records}\n"
        )

        log_handle.write(
            f"Invalid FASTA headers: "
            f"{invalid_headers}\n"
        )

        log_handle.write(
            f"Sequences with origins: "
            f"{sequences_with_origins}\n"
        )

        log_handle.write(
            f"Sequences without origins: "
            f"{missing_origin_sequences}\n"
        )

        log_handle.write(
            f"Origin rows written: "
            f"{origin_rows_written}\n"
        )

        log_handle.write(
            f"Habitat/taxonomy rows written: "
            f"{habitat_rows_written}\n"
        )

        log_handle.write(
            f"Sequences occurring more than once: "
            f"{duplicate_count}\n"
        )

        log_handle.write(
            f"Missing habitat files: "
            f"{missing_habitat_files}\n"
        )

        log_handle.write(
            f"Missing taxonomy files: "
            f"{missing_taxonomy_files}\n"
        )

        log_handle.write(
            f"First smORF ID: {first_id}\n"
        )

        log_handle.write(
            f"Last smORF ID: {last_id}\n"
        )

        log_handle.write(
            "Detailed duplicate counts:\n"
        )

        for sequence in smorf_order:
            count = sequence_counts[sequence]

            if count <= 1:
                continue

            smorf_id = sequence_to_id[sequence]

            raw_smorf_ids = sorted(
                {
                    entry[4]
                    for entry in sequence_metadata[sequence]
                }
            )

            log_handle.write(
                f"Sequence with ID {smorf_id} "
                f"appeared {count} times "
                f"(smORFs: {', '.join(raw_smorf_ids)})\n"
            )

    print("Replacing final output files...")

    os.replace(output_fasta_part, OUTPUT_FASTA)
    os.replace(origins_tsv_part, ORIGINS_TSV)
    os.replace(
        habitat_taxonomy_part,
        HABITAT_TAXONOMY_TSV
    )
    os.replace(log_file_part, LOG_FILE)

    print("Urban Soil smORF resource completed")
    print(f"Unique smORFs:          {len(sequence_to_id):,}")
    print(f"Origin rows:            {origin_rows_written:,}")
    print(f"FASTA:                  {OUTPUT_FASTA}")
    print(f"Origins TSV:            {ORIGINS_TSV}")
    print(
        f"Habitat/taxonomy TSV:   "
        f"{HABITAT_TAXONOMY_TSV}"
    )
    print(f"Log:                    {LOG_FILE}")


if __name__ == "__main__":
    main()
