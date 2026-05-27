#!/usr/bin/env python3

import gzip
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SMORF_FILE = Path(
    "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
)

GENE_CATALOG = Path(
    "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
)

EGGNOG_DIR = Path(
    "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"
)

BASE_OUT = Path("/work/microbiome/users/kruthi")

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

    print("Loading gene catalog...")

    gene_df = pd.read_csv(
        GENE_CATALOG,
        sep="\t",
        compression="xz",
        header=0,
        dtype=str,
        low_memory=False,
    )

    gene_df = gene_df[
        ["ORF", "Sample", "Original_ID", "Start", "End", "Strand", "Partial"]
    ].copy()

    gene_df = gene_df.rename(columns={"Original_ID": "Contig"})

    gene_df["Start"] = pd.to_numeric(gene_df["Start"], errors="coerce")
    gene_df["End"] = pd.to_numeric(gene_df["End"], errors="coerce")

    gene_df = gene_df.dropna(subset=["Start", "End"])

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

    for (sample, bc), sub in gene_df.groupby(
        ["Sample", "BaseContig"],
        sort=False
    ):
        gene_index[(sample, bc)] = sub.reset_index(drop=True)

    print(f"Gene index built: {len(gene_index)} contig groups")

    return gene_index


def get_neighbours(gene_sub_df, smorf_start, smorf_end, smorf_id):

    if gene_sub_df is None or gene_sub_df.empty:
        return None

    overlap = (
        np.minimum(gene_sub_df["End"].values, smorf_end)
        - np.maximum(gene_sub_df["Start"].values, smorf_start)
    )

    overlap = np.clip(overlap, 0, None)

    # CASE 1: overlapping gene exists
    if overlap.max() > 0:

        anchor_idx = int(overlap.argmax())

        note = (
            f"TARGET ({smorf_id}) - overlaps by "
            f"{int(overlap.max())} bp"
        )

    # CASE 2: intergenic smORF
    else:

        anchor_idx = int(
            np.searchsorted(
                gene_sub_df["Start"].values,
                smorf_start
            )
        )

        anchor_idx = max(
            0,
            min(anchor_idx, len(gene_sub_df) - 1)
        )

        note = f"TARGET ({smorf_id}) - intergenic"

    left = max(0, anchor_idx - WINDOW)

    right = min(
        len(gene_sub_df),
        anchor_idx + WINDOW + 1
    )

    neigh = gene_sub_df.iloc[left:right].copy()

    neigh["Note"] = ""

    neigh.iloc[
        anchor_idx - left,
        neigh.columns.get_loc("Note")
    ] = note

    return neigh


def load_eggnog_cache_for_sample(sample):

    annot_file = EGGNOG_DIR / sample / f"{sample}.emapper.annotations"

    if not annot_file.exists():
        return None

    query_to_lines = defaultdict(list)

    with annot_file.open("r") as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            query_id = line.split("\t", 1)[0]

            query_to_lines[query_id].append(line)

    return query_to_lines


def main():

    BASE_OUT.mkdir(parents=True, exist_ok=True)

    print("Pass 1: counting SmORF occurrences...")

    smorf_counts = Counter()

    with gzip.open(SMORF_FILE, "rt") as f:

        next(f)

        for line in f:

            smorf_id = line.split("\t", 1)[0]

            smorf_counts[smorf_id] += 1

    print(f"Unique SmORFs: {len(smorf_counts)}")

    max_count = max(smorf_counts.values())

    all_bins = []

    for start in range(1, max_count + 1, BIN_WIDTH):

        end = start + BIN_WIDTH - 1

        all_bins.append((start, end))

    bins_to_smorfs = defaultdict(list)

    for smorf_id, count in smorf_counts.items():

        start, end = bin_range(count)

        bins_to_smorfs[(start, end)].append(smorf_id)

    log_path = BASE_OUT / "bin_creation_log.tsv"

    with log_path.open("w") as log:

        log.write(
            "BinStart\tBinEnd\tLabel\tSmORF_Count\tStatus\tFolder\n"
        )

        for start, end in all_bins:

            label = bin_label(start, end)

            n = len(bins_to_smorfs.get((start, end), []))

            if n == 0:

                log.write(
                    f"{start}\t{end}\t{label}\t0\tSKIPPED_EMPTY\tNA\n"
                )

            else:

                folder = f"SmORF_neighbourhoods_{label}"

                log.write(
                    f"{start}\t{end}\t{label}\t{n}\tCREATED\t{folder}\n"
                )

    print(f"Bin log written: {log_path}")

    gene_index = load_gene_index()

    print("Pass 2: loading SmORF rows...")

    smorf_rows = defaultdict(list)

    with gzip.open(SMORF_FILE, "rt") as f:

        next(f)

        for line in f:

            parts = line.rstrip("\n").split("\t")

            smorf_id = parts[0]
            sample = parts[1]
            contig = parts[2]
            coords = parts[3]

            smorf_rows[smorf_id].append(
                (sample, contig, coords)
            )

    eggnog_cache = {}

    for start, end in all_bins:

        selected_list = sorted(
            bins_to_smorfs.get((start, end), [])
        )

        if not selected_list:
            continue

        label = bin_label(start, end)

        out_root = BASE_OUT / f"SmORF_neighbourhoods_{label}"

        out_root.mkdir(parents=True, exist_ok=True)

        print("\n" + "=" * 80)
        print(f"PROCESSING BIN: {label}")
        print(f"SmORFs in bin: {len(selected_list)}")
        print("=" * 80)

        summary_path = out_root / "summary.tsv"

        with summary_path.open("w") as s:

            s.write(
                "SmORF_ID\t"
                "TotalOccurrences\t"
                "OccurrencesProcessed\t"
                "NeighborsWritten\t"
                "SkippedNoGenes\t"
                "EggNOGMissing\n"
            )

        done = 0

        total_occ_processed = 0
        total_neighbors_written = 0

        for smorf_id in selected_list:

            occs = smorf_rows.get(smorf_id, [])

            total_occ = smorf_counts[smorf_id]

            processed = 0
            written = 0
            skipped_no_genes = 0
            eggnog_missing = 0

            for sample, contig, coords in occs:

                processed += 1
                total_occ_processed += 1

                bc = base_contig_name(contig)

                gene_sub = gene_index.get((sample, bc))

                if gene_sub is None or gene_sub.empty:

                    skipped_no_genes += 1
                    continue

                smorf_start, smorf_end = parse_coords(coords)

                neigh = get_neighbours(
                    gene_sub,
                    smorf_start,
                    smorf_end,
                    smorf_id
                )

                if neigh is None or neigh.empty:

                    skipped_no_genes += 1
                    continue

                out_dir = out_root / smorf_id / f"{sample}_{bc}"

                out_dir.mkdir(parents=True, exist_ok=True)

                neigh[
                    [
                        "ORF",
                        "Sample",
                        "Contig",
                        "Start",
                        "End",
                        "Strand",
                        "Partial",
                        "Note"
                    ]
                ].to_csv(
                    out_dir / "neighbors.tsv",
                    sep="\t",
                    index=False
                )

                if sample not in eggnog_cache:

                    eggnog_cache[sample] = (
                        load_eggnog_cache_for_sample(sample)
                    )

                query_map = eggnog_cache.get(sample)

                egg_out = out_dir / "EggNOG_annotations.tsv"

                with egg_out.open("w") as out:

                    out.write(
                        "#query\tORF\tSample\tContig\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
                    )

                    if query_map is None:

                        eggnog_missing += 1

                    else:

                        for _, r in neigh.iterrows():

                            orf = r["ORF"]
                            query_id = r["Contig"]

                            lines = query_map.get(query_id, [])

                            if lines:

                                out.write(
                                    f"{smorf_id}\t"
                                    f"{orf}\t"
                                    f"{sample}\t"
                                    f"{query_id}\t"
                                    f"{lines[0]}"
                                )

                written += 1
                total_neighbors_written += 1

            with summary_path.open("a") as s:

                s.write(
                    f"{smorf_id}\t"
                    f"{total_occ}\t"
                    f"{processed}\t"
                    f"{written}\t"
                    f"{skipped_no_genes}\t"
                    f"{eggnog_missing}\n"
                )

            done += 1

            if done % PROGRESS_EVERY == 0 or done == len(selected_list):

                print(
                    f"[{label}] "
                    f"{done}/{len(selected_list)} SmORFs done | "
                    f"occurrences processed={total_occ_processed} | "
                    f"neighbors written={total_neighbors_written}"
                )

        print(f"Completed bin: {label}")

    print("\nALL DONE")


if __name__ == "__main__":
    main()
