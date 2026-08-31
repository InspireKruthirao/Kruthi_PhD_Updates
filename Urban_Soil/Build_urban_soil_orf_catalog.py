#!/usr/bin/env python3

import gzip
import lzma
import re
from pathlib import Path


# Main Urban Soil directory
BASEDIR = Path(
    "/work/microbiome/users/kruthi/intermediate_results/urban_soil"
)

# Folder for combined Prodigal catalogue
OUTDIR = BASEDIR / "Prodigal"

# Final output files
ORF_FNA = OUTDIR / "UrbanSoil.ORF.fna.xz"
ORF_ORIG = OUTDIR / "UrbanSoil.ORF.orig.tsv.xz"


# Read FASTA files
def fasta_iter(path):

    path = Path(path)

    # Open gzipped or normal FASTA
    if str(path).endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    header = None
    seq_parts = []

    with opener(path, "rt") as f:

        for line in f:

            line = line.rstrip("\n")

            if not line:
                continue

            # New FASTA sequence
            if line.startswith(">"):

                if header is not None:
                    yield header, "".join(seq_parts)

                header = line[1:].strip()
                seq_parts = []

            else:
                seq_parts.append(line.strip())

        # Return final sequence
        if header is not None:
            yield header, "".join(seq_parts)


# Extract information from Prodigal FASTA header
def parse_prodigal_header(header):

    parts = [x.strip() for x in header.split("#")]

    if len(parts) < 4:
        raise ValueError(
            f"Unexpected Prodigal header: {header}"
        )

    original_id = parts[0].split()[0]
    start = int(parts[1])
    end = int(parts[2])
    strand = int(parts[3])

    # Get partial ORF information
    match = re.search(r"partial=([0-9]+)", header)
    partial = match.group(1) if match else "NA"

    return original_id, start, end, strand, partial


# Create Urban Soil ORF IDs
def make_orf_id(index):

    # Example:
    # US.ORF.000_000_000

    a = index // 1_000_000
    b = (index // 1_000) % 1_000
    c = index % 1_000

    return f"US.ORF.{a:03d}_{b:03d}_{c:03d}"


def main():

    # Create output directory if needed
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Find all sample Prodigal ORF files
    fnafiles = sorted(
        p for p in BASEDIR.glob("*/*_orfs.fna")
        if p.is_file()
    )

    if not fnafiles:
        raise FileNotFoundError(
            f"No *_orfs.fna files found under {BASEDIR}"
        )

    print(
        f"Found {len(fnafiles)} sample ORF files",
        flush=True
    )

    orf_index = 0

    # Open compressed output files
    with lzma.open(ORF_FNA, "wt") as out_fna, \
            lzma.open(ORF_ORIG, "wt") as out_orig:

        # Write metadata header
        out_orig.write(
            "ORF\tSample\tOriginal_ID\t"
            "Start\tEnd\tStrand\tPartial\n"
        )

        # Process each sample
        for i, fnafile in enumerate(fnafiles, 1):

            sample = fnafile.parent.name

            print(
                f"Processing {i}/{len(fnafiles)}: {sample}",
                flush=True
            )

            # Read every ORF
            for header, seq in fasta_iter(fnafile):

                original_id, start, end, strand, partial = \
                    parse_prodigal_header(header)

                # Assign new Urban Soil ORF ID
                orf_id = make_orf_id(orf_index)
                orf_index += 1

                # Write ORF origin information
                out_orig.write(
                    f"{orf_id}\t"
                    f"{sample}\t"
                    f"{original_id}\t"
                    f"{start}\t"
                    f"{end}\t"
                    f"{strand}\t"
                    f"{partial}\n"
                )

                # Write ORF sequence
                out_fna.write(
                    f">{orf_id}\n{seq}\n"
                )

    # Final summary
    print(
        f"Total ORFs: {orf_index}",
        flush=True
    )

    print(
        f"Created: {ORF_FNA}",
        flush=True
    )

    print(
        f"Created: {ORF_ORIG}",
        flush=True
    )


if __name__ == "__main__":
    main()
