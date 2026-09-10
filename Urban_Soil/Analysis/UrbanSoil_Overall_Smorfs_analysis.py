
#!/usr/bin/env python3

import gzip
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SMORF_FILE = Path(
    "/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil/"
    "UrbanSoil_SMORF_resource/"
    "UrbanSoil_100AA_SMORFs_origins.tsv.gz"
)

GENE_CATALOG = Path(
    "/work/microbiome/users/kruthi/intermediate_results/"
    "urban_soil/Gene_Catalog_Urban_Soil/"
    "UrbanSoil.ORF.orig.tsv.xz"
)

EGGNOG_DIR = Path(
    "/work/microbiome/users/kruthi/intermediate_results/"
    "egg_nog/results"
)

BASE_OUT = Path(
    "/work/microbiome/users/kruthi/intermediate_results/"
    "urban_soil/Neighbourhood_Analysis"
)

BIN_WIDTH = 5
WINDOW = 5
PROGRESS_EVERY = 25


def parse_coords(coord_string):
    start, end = coord_string.split("-", 1)
    return int(start), int(end)


def base_contig_name(contig):
    return contig.split("_polypolish", 1)[0]


def bin_range(n):
    start = ((n - 1) // BIN_WIDTH) * BIN_WIDTH + 1
    end = start + BIN_WIDTH - 1
    return start, end


def bin_label(start, end):
    return f"{start}_{end}"


def load_gene_index():

    print("Loading gene catalog...", flush=True)

    gene_df = pd.read_csv(
        GENE_CATALOG,
        sep="\t",
        compression="xz",
        header=0,
        dtype=str,
        low_memory=False,
    )

    gene_df = gene_df[
        [
            "ORF",
            "Sample",
            "Original_ID",
            "Start",
            "End",
            "Strand",
            "Partial",
        ]
    ].copy()

    gene_df = gene_df.rename(
        columns={"Original_ID": "Contig"}
    )

    gene_df["Start"] = pd.to_numeric(
        gene_df["Start"],
        errors="coerce",
    )

    gene_df["End"] = pd.to_numeric(
        gene_df["End"],
        errors="coerce",
    )

    gene_df = gene_df.dropna(
        subset=["Start", "End"]
    )

    gene_df["Start"] = gene_df["Start"].astype(int)
    gene_df["End"] = gene_df["End"].astype(int)

    gene_df["BaseContig"] = (
        gene_df["Contig"]
        .str.split("_polypolish", n=1)
        .str[0]
    )

    gene_df = gene_df.sort_values(
        ["Sample", "BaseContig", "Start"]
    ).reset_index(drop=True)

    gene_index = {}

    for (sample, base_contig), sub in gene_df.groupby(
        ["Sample", "BaseContig"],
        sort=False,
    ):
        gene_index[(sample, base_contig)] = (
            sub.reset_index(drop=True)
        )

    print(
        f"Gene index built: {len(gene_index)} contig groups",
        flush=True,
    )

    return gene_index


def get_neighbours(
    gene_sub_df,
    smorf_start,
    smorf_end,
    smorf_id,
):

    if gene_sub_df is None or gene_sub_df.empty:
        return None

    overlap = (
        np.minimum(
            gene_sub_df["End"].values,
            smorf_end,
        )
        - np.maximum(
            gene_sub_df["Start"].values,
            smorf_start,
        )
    )

    overlap = np.clip(overlap, 0, None)

    # Case 1: the smORF overlaps a predicted gene.
    if overlap.max() > 0:

        anchor_idx = int(overlap.argmax())

        note = (
            f"TARGET ({smorf_id}) - overlaps by "
            f"{int(overlap.max())} bp"
        )

    # Case 2: the smORF is intergenic.
    else:

        anchor_idx = int(
            np.searchsorted(
                gene_sub_df["Start"].values,
                smorf_start,
            )
        )

        anchor_idx = max(
            0,
            min(
                anchor_idx,
                len(gene_sub_df) - 1,
            ),
        )

        note = f"TARGET ({smorf_id}) - intergenic"

    left = max(
        0,
        anchor_idx - WINDOW,
    )

    right = min(
        len(gene_sub_df),
        anchor_idx + WINDOW + 1,
    )

    neighbours = gene_sub_df.iloc[left:right].copy()

    neighbours["Note"] = ""

    neighbours.iloc[
        anchor_idx - left,
        neighbours.columns.get_loc("Note"),
    ] = note

    return neighbours


def load_eggnog_cache_for_sample(sample):

    annotation_file = (
        EGGNOG_DIR
        / sample
        / f"{sample}.emapper.annotations.gz"
    )

    if not annotation_file.exists():
        return None

    query_to_lines = defaultdict(list)

    with gzip.open(
        annotation_file,
        "rt",
    ) as handle:

        for line in handle:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            query_id = line.split("\t", 1)[0]

            query_to_lines[query_id].append(line)

    return query_to_lines


def validate_inputs():

    missing = []

    if not SMORF_FILE.is_file():
        missing.append(str(SMORF_FILE))

    if not GENE_CATALOG.is_file():
        missing.append(str(GENE_CATALOG))

    if not EGGNOG_DIR.is_dir():
        missing.append(str(EGGNOG_DIR))

    if missing:

        print("ERROR: required inputs are missing:")

        for path in missing:
            print(f"  {path}")

        raise FileNotFoundError(
            "One or more required inputs are missing."
        )

    annotation_files = list(
        EGGNOG_DIR.glob(
            "*/*.emapper.annotations.gz"
        )
    )

    if not annotation_files:
        raise FileNotFoundError(
            f"No compressed eggNOG annotation files found in "
            f"{EGGNOG_DIR}"
        )

    print(
        f"eggNOG annotation files found: "
        f"{len(annotation_files)}",
        flush=True,
    )


def main():

    validate_inputs()

    BASE_OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "Pass 1: counting SmORF occurrences...",
        flush=True,
    )

    smorf_counts = Counter()

    with gzip.open(
        SMORF_FILE,
        "rt",
    ) as handle:

        header = next(handle, None)

        if header is None:
            raise ValueError(
                f"SmORF origins file is empty: {SMORF_FILE}"
            )

        for line in handle:

            if not line.strip():
                continue

            smorf_id = line.split("\t", 1)[0]

            smorf_counts[smorf_id] += 1

    if not smorf_counts:
        raise ValueError(
            "No SmORF records were read from the origins file."
        )

    print(
        f"Unique SmORFs: {len(smorf_counts)}",
        flush=True,
    )

    max_count = max(smorf_counts.values())

    all_bins = []

    for start in range(
        1,
        max_count + 1,
        BIN_WIDTH,
    ):

        end = start + BIN_WIDTH - 1

        all_bins.append((start, end))

    bins_to_smorfs = defaultdict(list)

    for smorf_id, count in smorf_counts.items():

        start, end = bin_range(count)

        bins_to_smorfs[(start, end)].append(smorf_id)

    log_path = BASE_OUT / "bin_creation_log.tsv"

    with log_path.open("w") as log:

        log.write(
            "BinStart\t"
            "BinEnd\t"
            "Label\t"
            "SmORF_Count\t"
            "Status\t"
            "Folder\n"
        )

        for start, end in all_bins:

            label = bin_label(start, end)

            number_of_smorfs = len(
                bins_to_smorfs.get(
                    (start, end),
                    [],
                )
            )

            if number_of_smorfs == 0:

                log.write(
                    f"{start}\t"
                    f"{end}\t"
                    f"{label}\t"
                    f"0\t"
                    f"SKIPPED_EMPTY\t"
                    f"NA\n"
                )

            else:

                folder = (
                    f"SmORF_neighbourhoods_{label}"
                )

                log.write(
                    f"{start}\t"
                    f"{end}\t"
                    f"{label}\t"
                    f"{number_of_smorfs}\t"
                    f"CREATED\t"
                    f"{folder}\n"
                )

    print(
        f"Bin log written: {log_path}",
        flush=True,
    )

    gene_index = load_gene_index()

    print(
        "Pass 2: loading SmORF rows...",
        flush=True,
    )

    smorf_rows = defaultdict(list)

    with gzip.open(
        SMORF_FILE,
        "rt",
    ) as handle:

        next(handle, None)

        for line in handle:

            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 4:
                continue

            smorf_id = parts[0]
            sample = parts[1]
            contig = parts[2]
            coordinates = parts[3]

            smorf_rows[smorf_id].append(
                (
                    sample,
                    contig,
                    coordinates,
                )
            )

    eggnog_cache = {}

    for start, end in all_bins:

        selected_list = sorted(
            bins_to_smorfs.get(
                (start, end),
                [],
            )
        )

        if not selected_list:
            continue

        label = bin_label(start, end)

        output_root = (
            BASE_OUT
            / f"SmORF_neighbourhoods_{label}"
        )

        output_root.mkdir(
            parents=True,
            exist_ok=True,
        )

        print(
            "\n" + "=" * 80,
            flush=True,
        )

        print(
            f"PROCESSING BIN: {label}",
            flush=True,
        )

        print(
            f"SmORFs in bin: {len(selected_list)}",
            flush=True,
        )

        print(
            "=" * 80,
            flush=True,
        )

        summary_path = output_root / "summary.tsv"

        with summary_path.open("w") as summary:

            summary.write(
                "SmORF_ID\t"
                "TotalOccurrences\t"
                "OccurrencesProcessed\t"
                "NeighborsWritten\t"
                "SkippedNoGenes\t"
                "EggNOGMissing\n"
            )

        completed_smorfs = 0
        total_occurrences_processed = 0
        total_neighbours_written = 0

        for smorf_id in selected_list:

            occurrences = smorf_rows.get(
                smorf_id,
                [],
            )

            total_occurrences = smorf_counts[smorf_id]

            processed = 0
            written = 0
            skipped_no_genes = 0
            eggnog_missing = 0

            for sample, contig, coordinates in occurrences:

                processed += 1
                total_occurrences_processed += 1

                base_contig = base_contig_name(contig)

                gene_subset = gene_index.get(
                    (
                        sample,
                        base_contig,
                    )
                )

                if (
                    gene_subset is None
                    or gene_subset.empty
                ):

                    skipped_no_genes += 1
                    continue

                smorf_start, smorf_end = parse_coords(
                    coordinates
                )

                neighbours = get_neighbours(
                    gene_subset,
                    smorf_start,
                    smorf_end,
                    smorf_id,
                )

                if (
                    neighbours is None
                    or neighbours.empty
                ):

                    skipped_no_genes += 1
                    continue

                output_directory = (
                    output_root
                    / smorf_id
                    / f"{sample}_{base_contig}"
                )

                output_directory.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                neighbours[
                    [
                        "ORF",
                        "Sample",
                        "Contig",
                        "Start",
                        "End",
                        "Strand",
                        "Partial",
                        "Note",
                    ]
                ].to_csv(
                    output_directory
                    / "neighbors.tsv",
                    sep="\t",
                    index=False,
                )

                if sample not in eggnog_cache:

                    eggnog_cache[sample] = (
                        load_eggnog_cache_for_sample(
                            sample
                        )
                    )

                query_map = eggnog_cache.get(sample)

                eggnog_output = (
                    output_directory
                    / "EggNOG_annotations.tsv"
                )

                with eggnog_output.open("w") as output:

                    output.write(
                        "#query\t"
                        "ORF\t"
                        "Sample\t"
                        "Contig\t"
                        "seed_ortholog\t"
                        "evalue\t"
                        "score\t"
                        "eggNOG_OGs\t"
                        "max_annot_lvl\t"
                        "COG_category\t"
                        "Description\t"
                        "Preferred_name\t"
                        "GOs\t"
                        "EC\t"
                        "KEGG_ko\t"
                        "KEGG_Pathway\t"
                        "KEGG_Module\t"
                        "KEGG_Reaction\t"
                        "KEGG_rclass\t"
                        "BRITE\t"
                        "KEGG_TC\t"
                        "CAZy\t"
                        "BiGG_Reaction\t"
                        "PFAMs\n"
                    )

                    if query_map is None:

                        eggnog_missing += 1

                    else:

                        for _, row in neighbours.iterrows():

                            orf = row["ORF"]
                            query_id = row["Contig"]

                            annotation_lines = (
                                query_map.get(
                                    query_id,
                                    [],
                                )
                            )

                            if annotation_lines:

                                output.write(
                                    f"{smorf_id}\t"
                                    f"{orf}\t"
                                    f"{sample}\t"
                                    f"{query_id}\t"
                                    f"{annotation_lines[0]}"
                                )

                written += 1
                total_neighbours_written += 1

            with summary_path.open("a") as summary:

                summary.write(
                    f"{smorf_id}\t"
                    f"{total_occurrences}\t"
                    f"{processed}\t"
                    f"{written}\t"
                    f"{skipped_no_genes}\t"
                    f"{eggnog_missing}\n"
                )

            completed_smorfs += 1

            if (
                completed_smorfs % PROGRESS_EVERY == 0
                or completed_smorfs
                == len(selected_list)
            ):

                print(
                    f"[{label}] "
                    f"{completed_smorfs}/"
                    f"{len(selected_list)} SmORFs done | "
                    f"occurrences processed="
                    f"{total_occurrences_processed} | "
                    f"neighbors written="
                    f"{total_neighbours_written}",
                    flush=True,
                )

        print(
            f"Completed bin: {label}",
            flush=True,
        )

    print(
        "\nALL DONE",
        flush=True,
    )


if __name__ == "__main__":
    main()
