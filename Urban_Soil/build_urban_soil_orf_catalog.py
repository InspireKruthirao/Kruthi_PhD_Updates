#!/usr/bin/env python3

import gzip
import lzma
import re
from pathlib import Path

BASEDIR = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil")
OUTDIR = BASEDIR / "Prodigal"

ORF_FNA = OUTDIR / "UrbanSoil.ORF.fna.xz"
ORF_ORIG = OUTDIR / "UrbanSoil.ORF.orig.tsv.xz"


def fasta_iter(path):
    path = Path(path)

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

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)

                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())

        if header is not None:
            yield header, "".join(seq_parts)


def parse_prodigal_header(header):
    parts = [x.strip() for x in header.split("#")]

    if len(parts) < 4:
        raise ValueError(f"Unexpected Prodigal header: {header}")

    original_id = parts[0].split()[0]
    start = int(parts[1])
    end = int(parts[2])
    strand = int(parts[3])

    match = re.search(r"partial=([0-9]+)", header)
    partial = match.group(1) if match else "NA"

    return original_id, start, end, strand, partial


def make_orf_id(index):
    # SHD-style ID:
    # US.ORF.000_000_000
    a = index // 1_000_000
    b = (index // 1_000) % 1_000
    c = index % 1_000

    return f"US.ORF.{a:03d}_{b:03d}_{c:03d}"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fnafiles = sorted(
        p for p in BASEDIR.glob("*/*_orfs.fna")
        if p.is_file()
    )

    if not fnafiles:
        raise FileNotFoundError(
            f"No *_orfs.fna files found under {BASEDIR}"
        )

    print(f"Found {len(fnafiles)} sample ORF files", flush=True)

    orf_index = 0

    with lzma.open(ORF_FNA, "wt") as out_fna, \
            lzma.open(ORF_ORIG, "wt") as out_orig:

        out_orig.write(
            "ORF\tSample\tOriginal_ID\tStart\tEnd\tStrand\tPartial\n"
        )

        for i, fnafile in enumerate(fnafiles, 1):
            sample = fnafile.parent.name

            print(
                f"Processing {i}/{len(fnafiles)}: {sample}",
                flush=True
            )

            for header, seq in fasta_iter(fnafile):
                original_id, start, end, strand, partial = \
                    parse_prodigal_header(header)

                orf_id = make_orf_id(orf_index)
                orf_index += 1

                out_orig.write(
                    f"{orf_id}\t"
                    f"{sample}\t"
                    f"{original_id}\t"
                    f"{start}\t"
                    f"{end}\t"
                    f"{strand}\t"
                    f"{partial}\n"
                )

                out_fna.write(
                    f">{orf_id}\n{seq}\n"
                )

    print(f"Total ORFs: {orf_index}", flush=True)
    print(f"Created: {ORF_FNA}", flush=True)
    print(f"Created: {ORF_ORIG}", flush=True)


if __name__ == "__main__":
    main()
