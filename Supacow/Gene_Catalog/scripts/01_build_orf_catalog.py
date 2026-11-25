#!/usr/bin/env python3
"""
STEP 01 - Gather the Prodigal genes of every MAG into one ORF catalogue
(as Urban_Soil/Gene_Catalog/01_Build_urban_soil_orf_catalog.py), with stable
SC.ORF.000_000_000 IDs. MAGs in sorted order, genes in Prodigal order.

    01_build_orf_catalog.py PRODIGAL_DIR N_MAGS PREFIX OUT_FNA OUT_FAA OUT_ORIG

Outputs
    OUT_FNA   ORF nucleotide sequences (plain FASTA, for MMseqs2)
    OUT_FAA   ORF protein sequences, same IDs (trailing * removed)
    OUT_ORIG  ORF -> MAG / Prodigal ID / start / end / strand / partial (.tsv.xz)
"""

import re
import sys
from pathlib import Path

from common import AtomicWriter, check, fasta_iter, grouped_id


def parse_prodigal_header(header):
    # contig_1_5 # 133 # 495 # 1 # ID=1_5;partial=00;start_type=ATG;...
    parts = [x.strip() for x in header.split("#")]
    if len(parts) < 5:
        raise ValueError(f"Unexpected Prodigal header: {header}")
    m = re.search(r"partial=([0-9]+)", header)
    return parts[0].split()[0], int(parts[1]), int(parts[2]), int(parts[3]), m.group(1) if m else "NA"


def main():
    prodigal_dir, n_mags, prefix, out_fna, out_faa, out_orig = sys.argv[1:]
    mags = sorted(p.name for p in Path(prodigal_dir).iterdir() if (p / f"{p.name}_orfs.fna").exists())
    print(f"MAGs with Prodigal output: {len(mags):,}", flush=True)
    check(len(mags) == int(n_mags), f"all {n_mags} MAGs have Prodigal output")

    n = 0
    with AtomicWriter(out_fna) as fna, AtomicWriter(out_faa) as faa, AtomicWriter(out_orig) as orig:
        orig.write("ORF\tMAG\tOriginal_ID\tStart\tEnd\tStrand\tPartial\n")
        for i, mag in enumerate(mags, 1):
            d = Path(prodigal_dir) / mag
            nuc = fasta_iter(d / f"{mag}_orfs.fna")
            prot = fasta_iter(d / f"{mag}_proteins.faa")
            for (h, seq), (hp, pseq) in zip(nuc, prot):
                original_id, start, end, strand, partial = parse_prodigal_header(h)
                if hp.split()[0] != original_id:
                    raise RuntimeError(f"{mag}: protein {hp.split()[0]} does not match gene {original_id}")
                orf = grouped_id(prefix, "ORF", n)
                n += 1
                orig.write(f"{orf}\t{mag}\t{original_id}\t{start}\t{end}\t{strand}\t{partial}\n")
                fna.write(f">{orf}\n{seq}\n")
                faa.write(f">{orf}\n{pseq.rstrip('*')}\n")
            if next(nuc, None) is not None or next(prot, None) is not None:
                raise RuntimeError(f"{mag}: different numbers of genes and proteins")
            if i % 500 == 0:
                print(f"  {i:,}/{len(mags):,} MAGs, {n:,} ORFs", flush=True)
    print(f"Total ORFs: {n:,}", flush=True)
    check(n > 0, "at least one ORF")


if __name__ == "__main__":
    main()
