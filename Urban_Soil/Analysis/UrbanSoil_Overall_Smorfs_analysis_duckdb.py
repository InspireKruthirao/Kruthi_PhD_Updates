SCRIPT=/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_DuckDB/UrbanSoil_Overall_Smorfs_analysis_duckdb.py

cat > "$SCRIPT" <<'PY'
#!/usr/bin/env python3

import gc
import gzip
from collections import Counter
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd


SMORF_FILE = Path("/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil/UrbanSoil_SMORF_resource/UrbanSoil_100AA_SMORFs_origins.tsv.gz")
GENE_CATALOG = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Gene_Catalog_Urban_Soil/UrbanSoil.ORF.orig.tsv.xz")
EGGNOG_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/egg_nog/results")
BASE_OUT = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_DuckDB")
DB_PATH = BASE_OUT / "urban_soil.duckdb"

BIN_WIDTH = 5
WINDOW = 5
PROGRESS_EVERY = 25
EGGNOG_BATCH_SIZE = 250_000


def validate_inputs():
    missing = [
        str(path)
        for path in (SMORF_FILE, GENE_CATALOG)
        if not path.is_file()
    ]

    if not EGGNOG_DIR.is_dir():
        missing.append(str(EGGNOG_DIR))

    if missing:
        raise FileNotFoundError(
            "Missing required input(s):\n" + "\n".join(missing)
        )


def create_table_from_dataframe(
    conn,
    table_name,
    view_name,
    dataframe,
):
    conn.register(view_name, dataframe)

    try:
        conn.execute(
            f"CREATE TABLE {table_name} "
            f"AS SELECT * FROM {view_name}"
        )
    finally:
        conn.unregister(view_name)


def insert_eggnog_batch(conn, rows):
    if not rows:
        return

    batch = pd.DataFrame(
        rows,
        columns=[
            "sample",
            "query_id",
            "annotation_line",
        ],
    )

    conn.register("eggnog_batch", batch)

    try:
        conn.execute(
            "INSERT INTO eggnog "
            "SELECT * FROM eggnog_batch"
        )
    finally:
        conn.unregister("eggnog_batch")


def database_is_complete():
    if not DB_PATH.is_file():
        return False

    try:
        conn = duckdb.connect(
            str(DB_PATH),
            read_only=True,
        )

        try:
            row = conn.execute(
                "SELECT completed "
                "FROM build_status "
                "LIMIT 1"
            ).fetchone()

            return bool(row and row[0])
        finally:
            conn.close()

    except Exception:
        return False


def load_data():
    print(
        "Initializing DuckDB database...",
        flush=True,
    )

    conn = duckdb.connect(str(DB_PATH))

    print(
        "Pass 1: Loading smORFs...",
        flush=True,
    )

    smorf_counts = Counter()
    smorf_data = []

    with gzip.open(SMORF_FILE, "rt") as handle:
        next(handle, None)

        for line in handle:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 4:
                continue

            smorf_id, sample, contig, coordinates = parts[:4]

            try:
                smorf_start, smorf_end = map(
                    int,
                    coordinates.split("-", 1),
                )
            except ValueError:
                continue

            base_contig = contig.split(
                "_polypolish",
                1,
            )[0]

            smorf_counts[smorf_id] += 1

            smorf_data.append(
                {
                    "smorf_id": smorf_id,
                    "sample": sample,
                    "base_contig": base_contig,
                    "smorf_start": smorf_start,
                    "smorf_end": smorf_end,
                }
            )

    print(
        f"Unique SmORFs: {len(smorf_counts)}",
        flush=True,
    )

    smorf_df = pd.DataFrame(smorf_data)

    create_table_from_dataframe(
        conn,
        "smorf_occurrences",
        "smorf_df",
        smorf_df,
    )

    conn.execute(
        "CREATE INDEX idx_smorf_id "
        "ON smorf_occurrences(smorf_id)"
    )

    conn.execute(
        "CREATE INDEX idx_smorf_sample_contig "
        "ON smorf_occurrences(sample, base_contig)"
    )

    conn.execute(
        """
        CREATE TABLE smorf_summary AS
        SELECT
            smorf_id,
            COUNT(*) AS occurrence_count
        FROM smorf_occurrences
        GROUP BY smorf_id
        """
    )

    del smorf_counts
    del smorf_data
    del smorf_df
    gc.collect()

    print(
        "Loading genes...",
        flush=True,
    )

    gene_df = pd.read_csv(
        GENE_CATALOG,
        sep="\t",
        compression="xz",
        dtype=str,
        low_memory=False,
    )

    required_columns = [
        "ORF",
        "Sample",
        "Original_ID",
        "Start",
        "End",
        "Strand",
        "Partial",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in gene_df.columns
    ]

    if missing_columns:
        raise ValueError(
            "Gene catalogue is missing columns: "
            + ", ".join(missing_columns)
        )

    gene_df = gene_df[required_columns].copy()

    gene_df.columns = [
        "orf",
        "sample",
        "annotation_query",
        "gene_start",
        "gene_end",
        "strand",
        "partial",
    ]

    gene_df["gene_start"] = pd.to_numeric(
        gene_df["gene_start"],
        errors="coerce",
    )

    gene_df["gene_end"] = pd.to_numeric(
        gene_df["gene_end"],
        errors="coerce",
    )

    gene_df = gene_df.dropna(
        subset=[
            "gene_start",
            "gene_end",
        ]
    )

    gene_df["gene_start"] = gene_df[
        "gene_start"
    ].astype("int64")

    gene_df["gene_end"] = gene_df[
        "gene_end"
    ].astype("int64")

    gene_df["base_contig"] = gene_df[
        "annotation_query"
    ].str.split(
        "_polypolish",
        n=1,
    ).str[0]

    genes_df = gene_df[
        [
            "orf",
            "sample",
            "base_contig",
            "annotation_query",
            "gene_start",
            "gene_end",
            "strand",
            "partial",
        ]
    ].copy()

    create_table_from_dataframe(
        conn,
        "genes",
        "genes_df",
        genes_df,
    )

    conn.execute(
        "CREATE INDEX idx_gene_sample_contig "
        "ON genes("
        "sample, base_contig, gene_start, gene_end"
        ")"
    )

    del gene_df
    del genes_df
    gc.collect()

    print(
        "Loading eggNOG files...",
        flush=True,
    )

    annotation_files = sorted(
        EGGNOG_DIR.glob(
            "*/*.emapper.annotations.gz"
        )
    )

    print(
        f"Found {len(annotation_files)} annotation files",
        flush=True,
    )

    if not annotation_files:
        raise FileNotFoundError(
            "No eggNOG annotation files found in "
            f"{EGGNOG_DIR}"
        )

    conn.execute(
        """
        CREATE TABLE eggnog (
            sample VARCHAR,
            query_id VARCHAR,
            annotation_line VARCHAR
        )
        """
    )

    annotation_header = None
    annotation_batch = []
    annotation_rows = 0

    for file_number, annotation_file in enumerate(
        annotation_files,
        1,
    ):
        sample = annotation_file.parent.name

        print(
            "  Loading eggNOG file "
            f"{file_number}/{len(annotation_files)}: "
            f"{sample}",
            flush=True,
        )

        with gzip.open(annotation_file, "rt") as handle:
            for line in handle:
                if line.startswith("#query"):
                    if annotation_header is None:
                        annotation_header = line.rstrip("\n")
                    continue

                if line.startswith("#") or not line.strip():
                    continue

                query_id = line.split("\t", 1)[0]

                annotation_batch.append(
                    (
                        sample,
                        query_id,
                        line.rstrip("\n"),
                    )
                )

                if (
                    len(annotation_batch)
                    >= EGGNOG_BATCH_SIZE
                ):
                    insert_eggnog_batch(
                        conn,
                        annotation_batch,
                    )

                    annotation_rows += len(
                        annotation_batch
                    )

                    annotation_batch.clear()

    if annotation_batch:
        insert_eggnog_batch(
            conn,
            annotation_batch,
        )

        annotation_rows += len(
            annotation_batch
        )

        annotation_batch.clear()

    conn.execute(
        "CREATE INDEX idx_eggnog_sample_query "
        "ON eggnog(sample, query_id)"
    )

    conn.execute(
        "CREATE TABLE eggnog_metadata("
        "annotation_header VARCHAR"
        ")"
    )

    conn.execute(
        "INSERT INTO eggnog_metadata VALUES (?)",
        [annotation_header or "#query"],
    )

    conn.execute(
        "CREATE TABLE build_status("
        "completed BOOLEAN"
        ")"
    )

    conn.execute(
        "INSERT INTO build_status VALUES (TRUE)"
    )

    print(
        f"Loaded {annotation_rows} eggNOG annotation rows",
        flush=True,
    )

    print(
        "DuckDB database ready!",
        flush=True,
    )

    return conn


def get_neighbours(
    conn,
    sample,
    base_contig,
    smorf_start,
    smorf_end,
    smorf_id,
):
    result = conn.execute(
        """
        SELECT
            orf,
            annotation_query,
            gene_start,
            gene_end,
            strand,
            partial
        FROM genes
        WHERE sample = ?
          AND base_contig = ?
        ORDER BY gene_start
        """,
        [
            sample,
            base_contig,
        ],
    ).fetchall()

    if not result:
        return None

    genes = [
        {
            "ORF": row[0],
            "AnnotationQuery": row[1],
            "Start": int(row[2]),
            "End": int(row[3]),
            "Strand": row[4],
            "Partial": row[5],
        }
        for row in result
    ]

    overlap = np.array(
        [
            min(
                gene["End"],
                smorf_end,
            )
            - max(
                gene["Start"],
                smorf_start,
            )
            for gene in genes
        ]
    )

    overlap = np.clip(
        overlap,
        0,
        None,
    )

    if overlap.max() > 0:
        anchor_index = int(overlap.argmax())

        note = (
            f"TARGET ({smorf_id}) - overlaps by "
            f"{int(overlap.max())} bp"
        )

    else:
        anchor_index = int(
            np.searchsorted(
                [
                    gene["Start"]
                    for gene in genes
                ],
                smorf_start,
            )
        )

        anchor_index = max(
            0,
            min(
                anchor_index,
                len(genes) - 1,
            ),
        )

        note = (
            f"TARGET ({smorf_id}) - intergenic"
        )

    left = max(
        0,
        anchor_index - WINDOW,
    )

    right = min(
        len(genes),
        anchor_index + WINDOW + 1,
    )

    neighbours = genes[left:right]

    for neighbour in neighbours:
        neighbour["Note"] = ""

    neighbours[
        anchor_index - left
    ]["Note"] = note

    return neighbours


def get_eggnog_annotations(
    conn,
    sample,
    neighbours,
):
    query_ids = [
        neighbour["AnnotationQuery"]
        for neighbour in neighbours
    ]

    if not query_ids:
        return {}, 0

    placeholders = ",".join(
        "?"
        for _ in query_ids
    )

    rows = conn.execute(
        "SELECT query_id, annotation_line "
        "FROM eggnog "
        "WHERE sample = ? "
        f"AND query_id IN ({placeholders})",
        [
            sample,
            *query_ids,
        ],
    ).fetchall()

    annotation_map = {
        query_id: annotation_line
        for query_id, annotation_line in rows
    }

    missing = sum(
        query_id not in annotation_map
        for query_id in query_ids
    )

    return annotation_map, missing


def main():
    validate_inputs()

    BASE_OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        BASE_OUT / "logs"
    ).mkdir(
        parents=True,
        exist_ok=True,
    )

    if (
        DB_PATH.exists()
        and not database_is_complete()
    ):
        raise RuntimeError(
            "Incomplete DuckDB database exists: "
            f"{DB_PATH}\n"
            "Move it aside before rerunning this script."
        )

    if database_is_complete():
        print(
            "Complete database already exists; "
            "skipping data load",
            flush=True,
        )

        conn = duckdb.connect(
            str(DB_PATH)
        )

    else:
        conn = load_data()

    annotation_header = conn.execute(
        "SELECT annotation_header "
        "FROM eggnog_metadata "
        "LIMIT 1"
    ).fetchone()[0]

    print(
        "Counting SmORF occurrences "
        "from the database...",
        flush=True,
    )

    summary_rows = conn.execute(
        "SELECT smorf_id, occurrence_count "
        "FROM smorf_summary"
    ).fetchall()

    smorf_counts = dict(summary_rows)

    max_count = max(
        smorf_counts.values()
    )

    all_bins = [
        (
            bin_start,
            bin_start + BIN_WIDTH - 1,
        )
        for bin_start in range(
            1,
            max_count + 1,
            BIN_WIDTH,
        )
    ]

    print(
        f"Unique SmORFs: {len(smorf_counts)}",
        flush=True,
    )

    log_path = (
        BASE_OUT / "bin_creation_log.tsv"
    )

    with log_path.open("w") as log:
        log.write(
            "BinStart\tBinEnd\tLabel\t"
            "SmORF_Count\tStatus\tFolder\n"
        )

        for bin_start, bin_end in all_bins:
            count = sum(
                bin_start
                <= occurrence_count
                <= bin_end
                for occurrence_count
                in smorf_counts.values()
            )

            label = (
                f"{bin_start}_{bin_end}"
            )

            status = (
                "CREATED"
                if count
                else "SKIPPED_EMPTY"
            )

            folder = (
                f"SmORF_neighbourhoods_{label}"
                if count
                else "NA"
            )

            log.write(
                f"{bin_start}\t{bin_end}\t"
                f"{label}\t{count}\t"
                f"{status}\t{folder}\n"
            )

    print(
        f"Bin log written: {log_path}",
        flush=True,
    )

    for bin_start, bin_end in all_bins:
        label = (
            f"{bin_start}_{bin_end}"
        )

        smorfs_in_bin = sorted(
            smorf_id
            for smorf_id, occurrence_count
            in smorf_counts.items()
            if (
                bin_start
                <= occurrence_count
                <= bin_end
            )
        )

        if not smorfs_in_bin:
            continue

        output_root = (
            BASE_OUT
            / f"SmORF_neighbourhoods_{label}"
        )

        output_root.mkdir(
            parents=True,
            exist_ok=True,
        )

        print(
            f"\nProcessing bin: {label} "
            f"({len(smorfs_in_bin)} smORFs)",
            flush=True,
        )

        summary_path = (
            output_root / "summary.tsv"
        )

        with summary_path.open("w") as summary:
            summary.write(
                "SmORF_ID\tTotalOccurrences\t"
                "OccurrencesProcessed\t"
                "NeighborsWritten\t"
                "SkippedNoGenes\t"
                "EggNOGMissing\n"
            )

        for index, smorf_id in enumerate(
            smorfs_in_bin,
            1,
        ):
            occurrences = conn.execute(
                """
                SELECT
                    sample,
                    base_contig,
                    smorf_start,
                    smorf_end
                FROM smorf_occurrences
                WHERE smorf_id = ?
                """,
                [smorf_id],
            ).fetchall()

            processed = 0
            written = 0
            skipped_no_genes = 0
            eggnog_missing = 0

            for (
                sample,
                base_contig,
                smorf_start,
                smorf_end,
            ) in occurrences:
                processed += 1

                neighbours = get_neighbours(
                    conn,
                    sample,
                    base_contig,
                    smorf_start,
                    smorf_end,
                    smorf_id,
                )

                if not neighbours:
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

                neighbours_path = (
                    output_directory
                    / "neighbors.tsv"
                )

                with neighbours_path.open("w") as out:
                    out.write(
                        "ORF\tSample\tContig\t"
                        "Start\tEnd\tStrand\t"
                        "Partial\tNote\n"
                    )

                    for neighbour in neighbours:
                        out.write(
                            f"{neighbour['ORF']}\t"
                            f"{sample}\t"
                            f"{base_contig}\t"
                            f"{neighbour['Start']}\t"
                            f"{neighbour['End']}\t"
                            f"{neighbour['Strand']}\t"
                            f"{neighbour['Partial']}\t"
                            f"{neighbour['Note']}\n"
                        )

                (
                    annotation_map,
                    missing,
                ) = get_eggnog_annotations(
                    conn,
                    sample,
                    neighbours,
                )

                eggnog_missing += missing

                annotation_path = (
                    output_directory
                    / "EggNOG_annotations.tsv"
                )

                with annotation_path.open("w") as out:
                    out.write(
                        annotation_header.rstrip("\n")
                        + "\n"
                    )

                    for neighbour in neighbours:
                        annotation_line = (
                            annotation_map.get(
                                neighbour[
                                    "AnnotationQuery"
                                ]
                            )
                        )

                        if annotation_line:
                            out.write(
                                annotation_line.rstrip(
                                    "\n"
                                )
                                + "\n"
                            )

                written += 1

            with summary_path.open("a") as summary:
                summary.write(
                    f"{smorf_id}\t"
                    f"{len(occurrences)}\t"
                    f"{processed}\t"
                    f"{written}\t"
                    f"{skipped_no_genes}\t"
                    f"{eggnog_missing}\n"
                )

            if (
                index % PROGRESS_EVERY == 0
                or index == len(smorfs_in_bin)
            ):
                print(
                    f"  [{label}] "
                    f"{index}/"
                    f"{len(smorfs_in_bin)} done",
                    flush=True,
                )

        print(
            f"Completed bin: {label}",
            flush=True,
        )

    conn.close()

    print(
        "\nALL DONE",
        flush=True,
    )


if __name__ == "__main__":
    main()
PY
