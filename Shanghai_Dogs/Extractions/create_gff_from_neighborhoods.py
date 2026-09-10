#!/usr/bin/env python3

"""
Generate GFF files for smORF gene-neighbourhood analysis.

Combines:
    - neighbors.tsv
    - EggNOG_annotations.tsv

to create a GFF file containing gene coordinates,
COG annotations and the target smORF.
"""

import sys
from pathlib import Path

import pandas as pd


# Input neighbourhood folder
BASE_DIR = Path(sys.argv[1])

missing_neighbors = []
missing_eggnog = []


def create_gff_for_contig(contig_dir):
    """Create a GFF file for one contig neighbourhood."""

    contig_path = Path(contig_dir)

    contig_name = contig_path.name.split("_", 1)[1]
    gff_file = contig_path / f"{contig_name}.gff"

    neighbors_file = contig_path / "neighbors.tsv"
    eggnog_file = contig_path / "EggNOG_annotations.tsv"

    # Check input files
    if not neighbors_file.exists():
        missing_neighbors.append(str(contig_path))
        return

    if not eggnog_file.exists():
        missing_eggnog.append(str(contig_path))
        return

    # Load neighbouring genes
    neighbors = pd.read_csv(
        neighbors_file,
        sep="\t"
    )

    # Find target smORF
    target_contig = None

    if "Note" in neighbors.columns:
        target_row = neighbors[
            neighbors["Note"].str.contains(
                "TARGET",
                na=False
            )
        ]

        if not target_row.empty:
            target_contig = target_row.iloc[0]["Contig"]

    # Read EggNOG header
    with open(eggnog_file, "r") as f:
        header = (
            f.readline()
            .lstrip("#")
            .strip()
            .split("\t")
        )

    # Load EggNOG annotations
    eggnog = pd.read_csv(
        eggnog_file,
        sep="\t",
        skiprows=1,
        names=header
    )

    # Map gene IDs to COG annotations
    annotation_map = {}

    for _, row in eggnog.iterrows():

        cog_name = "Unknown"
        eggnog_ogs = row.get("eggNOG_OGs")

        if (
            pd.notna(eggnog_ogs)
            and "COG" in str(eggnog_ogs)
        ):
            cog_raw = (
                str(eggnog_ogs)
                .split("@")[0]
                .split(",")[0]
            )

            if cog_raw.startswith("COG"):

                cog_category = row.get(
                    "COG_category",
                    ""
                )

                cat = (
                    str(cog_category).strip()
                    if pd.notna(cog_category)
                    else ""
                )

                if cat and cat != "-":
                    cog_name = f"{cog_raw}-{cat}"
                else:
                    cog_name = cog_raw

        annotation_map[row["Contig"]] = cog_name

    # Write GFF
    with open(gff_file, "w") as gff:

        gff.write("##gff-version 3\n")

        for _, row in neighbors.iterrows():

            contig_id = row["Contig"]

            name = annotation_map.get(
                contig_id,
                "Unknown"
            )

            score = (
                "1"
                if name != "Unknown"
                else "0.0"
            )

            strand = (
                "+"
                if row["Strand"] == 1
                else "-"
            )

            attributes = [
                f"ID={contig_id}",
                f"Name={name}"
            ]

            if contig_id == target_contig:

                attributes.extend([
                    "target=1",
                    "Note=TARGET_sMORF",
                    "colour=1"
                ])

            else:
                attributes.append(
                    "target=0"
                )

            # GFF attributes use ASCII semicolons
            gff.write(
                f"{contig_id}\t"
                f"Prodigal_v2.6.3\t"
                f"CDS\t"
                f"{row['Start']}\t"
                f"{row['End']}\t"
                f"{score}\t"
                f"{strand}\t"
                f"0\t"
                f"{';'.join(attributes)}\n"
            )


def main():
    """Process all smORFs in one neighbourhood bin."""

    print(f"Processing: {BASE_DIR}")

    smorf_dirs = sorted(
        BASE_DIR.glob(
            "SHD1_SM.100AA.*"
        )
    )

    print(
        f"smORFs found: "
        f"{len(smorf_dirs)}"
    )

    smorf_count = 0
    contig_count = 0

    for smorf_dir in smorf_dirs:

        contig_dirs = sorted(
            smorf_dir.glob(
                "D*_contig_*"
            )
        )

        for contig_dir in contig_dirs:

            create_gff_for_contig(
                contig_dir
            )

            contig_count += 1

        smorf_count += 1

        if smorf_count % 1000 == 0:

            print(
                f"{smorf_count:,} / "
                f"{len(smorf_dirs):,} "
                f"smORFs completed",
                flush=True
            )

    # Summary
    print("\nSUMMARY")

    print(
        f"Folder: "
        f"{BASE_DIR.name}"
    )

    print(
        f"smORFs processed: "
        f"{smorf_count:,}"
    )

    print(
        f"Contigs processed: "
        f"{contig_count:,}"
    )

    print(
        f"Missing neighbors.tsv: "
        f"{len(missing_neighbors):,}"
    )

    print(
        f"Missing EggNOG_annotations.tsv: "
        f"{len(missing_eggnog):,}"
    )

    print("Done.")


if __name__ == "__main__":
    main()
