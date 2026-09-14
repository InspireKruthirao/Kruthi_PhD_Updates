#!/usr/bin/env python3

import csv
import gzip
import lzma
import os
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

ROOT = Path(
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil/"
    "Neighbourhood_Analysis_SQL_ge6"
)

GENES = Path(
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil/"
    "Gene_Catalog_Urban_Soil/UrbanSoil.ORF.orig.tsv.xz"
)

EGGNOG = Path(
    "/work/microbiome/users/kruthi/intermediate_results/egg_nog/results"
)

PROGRESS = 10000


def gene_key(row):
    return (
        row["Sample"].strip(),
        row["ORF"].strip(),
        row["Contig"].strip(),
        str(int(row["Start"])),
        str(int(row["End"])),
    )


def read_neighbors(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "ORF",
            "Sample",
            "Contig",
            "Start",
            "End",
            "Strand",
            "Partial",
            "Note",
        }

        if (
            reader.fieldnames is None
            or not required.issubset(reader.fieldnames)
        ):
            raise ValueError(
                f"Bad header in {path}: {reader.fieldnames}"
            )

        return list(reader)


def neighbor_files():
    for directory, _, files in os.walk(ROOT):
        if "neighbors.tsv" in files:
            yield Path(directory) / "neighbors.tsv"


def get_cog(line):
    fields = line.rstrip("\n").split("\t")

    if len(fields) < 7:
        return "Unknown"

    category = fields[6].strip()

    for item in fields[4].split(","):
        cog = item.split("@", 1)[0].strip()

        if cog.startswith("COG"):
            if category and category != "-":
                return f"{cog}-{category}"

            return cog

    return "Unknown"


def strand(value):
    value = str(value).strip()

    if value in {"1", "+", "1.0"}:
        return "+"

    if value in {"-1", "-", "-1.0"}:
        return "-"

    return "."


def escape_gff(value):
    return quote(
        str(value),
        safe="._:-"
    )


def write_gff(
    neighbors_path,
    rows,
    original_ids,
    annotations,
):
    contig = rows[0]["Contig"].strip()

    output = (
        neighbors_path.parent
        / f"{contig}.gff"
    )

    # Allows a failed job to be resumed safely.
    if (
        output.is_file()
        and output.stat().st_size > 0
    ):
        return "existing", 0, 0

    annotated = 0
    missing_original = 0

    # Only one temporary file exists while each GFF is written.
    part = output.with_suffix(".gff.part")

    with part.open("w") as gff:
        gff.write("##gff-version 3\n")

        for row in rows:
            key = gene_key(row)
            original = original_ids.get(key)

            if original is None:
                missing_original += 1
                name = "Unknown"
            else:
                name = annotations.get(
                    original,
                    "Unknown",
                )

            if name != "Unknown":
                annotated += 1

            target = (
                "TARGET"
                in row.get("Note", "").upper()
            )

            attributes = [
                f"ID={escape_gff(row['ORF'].strip())}",
                f"Name={escape_gff(name)}",
                f"target={int(target)}",
            ]

            if original is not None:
                attributes.append(
                    f"Original_ID={escape_gff(original)}"
                )

            if target:
                attributes.extend([
                    "Note=TARGET_smORF",
                    "colour=1",
                ])

            score = (
                "1"
                if name != "Unknown"
                else "0.0"
            )

            gff.write(
                "\t".join([
                    contig,
                    "Prodigal_v2.6.3",
                    "CDS",
                    str(int(row["Start"])),
                    str(int(row["End"])),
                    score,
                    strand(row["Strand"]),
                    "0",
                    ";".join(attributes),
                ])
                + "\n"
            )

    part.replace(output)

    return (
        "written",
        annotated,
        missing_original,
    )


def validate_inputs():
    for path in (
        ROOT,
        GENES,
        EGGNOG,
    ):
        if not path.exists():
            raise FileNotFoundError(path)


def main():
    validate_inputs()

    print(
        "Pass 1: scanning neighbors.tsv files...",
        flush=True,
    )

    needed = set()
    paths_by_sample = defaultdict(list)

    file_count = 0
    row_count = 0

    for path in neighbor_files():
        rows = read_neighbors(path)

        if not rows:
            continue

        samples = {
            row["Sample"].strip()
            for row in rows
        }

        if len(samples) != 1:
            raise ValueError(
                f"Multiple samples in {path}: "
                f"{samples}"
            )

        sample = next(iter(samples))

        paths_by_sample[sample].append(path)
        file_count += 1

        for row in rows:
            needed.add(gene_key(row))
            row_count += 1

        if file_count % PROGRESS == 0:
            print(
                f"  files={file_count:,}; "
                f"unique genes={len(needed):,}",
                flush=True,
            )

    if file_count == 0:
        raise RuntimeError(
            f"No neighbors.tsv files found below {ROOT}"
        )

    print(
        f"Neighborhood files: {file_count:,}",
        flush=True,
    )

    print(
        f"Neighbor rows: {row_count:,}",
        flush=True,
    )

    print(
        f"Unique gene occurrences needed: "
        f"{len(needed):,}",
        flush=True,
    )

    print(
        "Pass 2: recovering Original_ID values...",
        flush=True,
    )

    original_ids = {}

    with lzma.open(
        GENES,
        "rt",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required = {
            "ORF",
            "Sample",
            "Original_ID",
            "Start",
            "End",
        }

        if (
            reader.fieldnames is None
            or not required.issubset(
                reader.fieldnames
            )
        ):
            raise ValueError(
                "Bad gene-catalogue header: "
                f"{reader.fieldnames}"
            )

        for number, row in enumerate(
            reader,
            1,
        ):
            original = (
                row["Original_ID"].strip()
            )

            base_contig = original.split(
                "_polypolish",
                1,
            )[0]

            key = (
                row["Sample"].strip(),
                row["ORF"].strip(),
                base_contig,
                str(int(row["Start"])),
                str(int(row["End"])),
            )

            if key in needed:
                original_ids[key] = original
                needed.remove(key)

            if number % 10000000 == 0:
                print(
                    f"  catalogue={number:,}; "
                    f"matched={len(original_ids):,}; "
                    f"missing={len(needed):,}",
                    flush=True,
                )

    print(
        f"Original IDs recovered: "
        f"{len(original_ids):,}",
        flush=True,
    )

    print(
        f"Unmatched gene occurrences: "
        f"{len(needed):,}",
        flush=True,
    )

    queries_by_sample = defaultdict(set)

    for key, original in original_ids.items():
        sample = key[0]
        queries_by_sample[sample].add(
            original
        )

    print(
        "Pass 3: reading eggNOG and writing GFFs...",
        flush=True,
    )

    written = 0
    existing = 0
    annotated_rows = 0
    missing_rows = 0

    samples = sorted(paths_by_sample)

    for sample_number, sample in enumerate(
        samples,
        1,
    ):
        required_queries = (
            queries_by_sample.get(
                sample,
                set(),
            )
        )

        eggnog_file = (
            EGGNOG
            / sample
            / f"{sample}.emapper.annotations.gz"
        )

        annotations = {}

        if eggnog_file.is_file():
            with gzip.open(
                eggnog_file,
                "rt",
            ) as handle:

                for line in handle:
                    if (
                        not line.strip()
                        or line.startswith("#")
                    ):
                        continue

                    query_id = line.split(
                        "\t",
                        1,
                    )[0]

                    if query_id in required_queries:
                        annotations[query_id] = (
                            get_cog(line)
                        )

        print(
            f"[{sample_number}/{len(samples)}] "
            f"{sample}: "
            f"neighborhoods="
            f"{len(paths_by_sample[sample]):,}; "
            f"queries="
            f"{len(required_queries):,}; "
            f"eggNOG matched="
            f"{len(annotations):,}",
            flush=True,
        )

        for number, path in enumerate(
            paths_by_sample[sample],
            1,
        ):
            rows = read_neighbors(path)

            status, annotated, missing = (
                write_gff(
                    path,
                    rows,
                    original_ids,
                    annotations,
                )
            )

            if status == "written":
                written += 1
            else:
                existing += 1

            annotated_rows += annotated
            missing_rows += missing

            if number % PROGRESS == 0:
                print(
                    f"  {sample}: "
                    f"{number:,}/"
                    f"{len(paths_by_sample[sample]):,}",
                    flush=True,
                )

        del annotations

    print(
        "\nGFF CREATION COMPLETE",
        flush=True,
    )

    print(
        f"GFF files written: {written:,}",
        flush=True,
    )

    print(
        f"Existing GFF files skipped: "
        f"{existing:,}",
        flush=True,
    )

    print(
        f"Annotated GFF rows: "
        f"{annotated_rows:,}",
        flush=True,
    )

    print(
        f"Rows missing Original_ID: "
        f"{missing_rows:,}",
        flush=True,
    )


if __name__ == "__main__":
    main()
