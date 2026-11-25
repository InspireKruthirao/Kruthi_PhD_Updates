#!/usr/bin/env python3
"""
STEP 09 - smORF resource from the GMSC-mapper runs (as
Urban_Soil/SmORF_Catalog/Create_resource_files.py).

Unique mapped smORF sequences (mapped.smorfs.faa of every chunk), most
common first (ties: by sequence), get pad9 IDs COW_SM.100AA.000_000_000.
Every predicted smORF (predicted.filtered.smorf.faa) with one of those
sequences is an origin. GMSC-mapper ran on chunks of MAGs whose contigs
were renamed <MAG>|<contig>, so the MAG comes from the contig name.

    09_smorf_resource.py GMSC_DIR N_CHUNKS PREFIX OUT_FAA OUT_ORIGINS OUT_HABITAT_TAXONOMY OUT_LOG

Outputs (same columns as Urban Soil)
    OUT_FAA               >COW_SM.100AA.###  sequence
    OUT_ORIGINS           SmORF ID, Sample ID (MAG), Contig, Coordinates, Strand
    OUT_HABITAT_TAXONOMY  SmORF ID, Sample ID (MAG), Habitat, Taxonomy
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

from common import AtomicWriter, check, fasta_iter, grouped_id


def clean(seq):
    return seq.replace(" ", "").replace("\r", "").upper()


def read_two_columns(path):
    """GMSC-mapper table (comment line, header, then rows) -> {smORF: value}"""
    out = {}
    if not path.exists():
        return out
    with open(path) as f:
        rows = (r for r in csv.reader(f, delimiter="\t") if r and not r[0].startswith("#"))
        next(rows, None)  # header
        for r in rows:
            if len(r) >= 2:
                out[r[0].strip()] = r[1].strip()
    return out


def main():
    gmsc_dir, n_chunks, prefix, out_faa, out_origins, out_ht, out_log = sys.argv[1:]
    chunks = [Path(gmsc_dir) / f"chunk{i:02d}" for i in range(int(n_chunks))]
    for c in chunks:
        for name in ("mapped.smorfs.faa", "predicted.filtered.smorf.faa"):
            if not (c / name).exists():
                raise FileNotFoundError(c / name)

    mapped = Counter()
    n_mapped = 0
    for c in chunks:
        for _, seq in fasta_iter(c / "mapped.smorfs.faa"):
            if clean(seq):
                mapped[clean(seq)] += 1
                n_mapped += 1
    order = sorted(mapped, key=lambda s: (-mapped[s], s))
    to_id = {s: grouped_id(prefix, "100AA", i) for i, s in enumerate(order)}
    print(f"mapped smORF records: {n_mapped:,}, unique sequences: {len(to_id):,}", flush=True)
    with AtomicWriter(out_faa) as out:
        for s in order:
            out.write(f">{to_id[s]}\n{s}\n")
    del order

    origins = defaultdict(list)  # sequence -> (MAG, contig, coordinates, strand, raw id, chunk)
    counts = Counter()
    habitat, taxonomy = {}, {}
    n_pred = n_unmatched = n_bad = 0
    for k, c in enumerate(chunks):
        for raw, value in read_two_columns(c / "habitat.out.smorfs.tsv").items():
            habitat[(raw, k)] = ",".join(sorted(h.strip() for h in value.split(",") if h.strip()))
        for raw, value in read_two_columns(c / "taxonomy.out.smorfs.tsv").items():
            taxonomy[(raw, k)] = value
        for header, seq in fasta_iter(c / "predicted.filtered.smorf.faa"):
            n_pred += 1
            s = clean(seq)
            if s not in to_id:
                n_unmatched += 1
                continue
            # smORF_000000 # <MAG>|<contig>_<gene> # start # end # strand # ID=...
            parts = [p.strip() for p in header.split("#")]
            if len(parts) < 5 or parts[4] not in ("1", "-1") or "|" not in parts[1]:
                n_bad += 1
                continue
            mag, contig_gene = parts[1].split("|", 1)
            contig = contig_gene.rsplit("_", 1)[0]
            origins[s].append((mag, contig, f"{parts[2]}-{parts[3]}", "+" if parts[4] == "1" else "-",
                               parts[0].split()[0], k))
            counts[s] += 1

    print(f"predicted smORF records: {n_pred:,}, matched to a mapped sequence: {sum(counts.values()):,}, "
          f"unmatched: {n_unmatched:,}, bad headers: {n_bad:,}", flush=True)
    check(n_bad == 0, "every predicted smORF header has MAG, contig, coordinates and strand")
    check(len(origins) == len(to_id), f"each of the {len(to_id):,} unique smORFs has at least one origin")

    n_rows = 0
    with AtomicWriter(out_origins) as o, AtomicWriter(out_ht) as h:
        o.write("SmORF ID\tSample ID\tContig\tCoordinates\tStrand\n")
        h.write("SmORF ID\tSample ID\tHabitat\tTaxonomy\n")
        for s in sorted(to_id, key=lambda s: (-counts[s], s)):
            for mag, contig, coords, strand, raw, k in origins[s]:
                o.write(f"{to_id[s]}\t{mag}\t{contig}\t{coords}\t{strand}\n")
                h.write(f"{to_id[s]}\t{mag}\t{habitat.get((raw, k), '')}\t{taxonomy.get((raw, k), '')}\n")
                n_rows += 1

    mags = {e[0] for v in origins.values() for e in v}
    with AtomicWriter(out_log) as log:
        log.write(f"GMSC-mapper chunks: {len(chunks)}\n")
        log.write(f"MAGs with a mapped smORF: {len(mags)}\n")
        log.write(f"Mapped smORF records: {n_mapped}\n")
        log.write(f"Unique smORF sequences (100AA): {len(to_id)}\n")
        log.write(f"Predicted smORF records: {n_pred}\n")
        log.write(f"Origin rows: {n_rows}\n")
        log.write(f"Unique smORFs found more than once: {sum(1 for s in to_id if counts[s] > 1)}\n")
    print(f"origin rows: {n_rows:,}; MAGs with a mapped smORF: {len(mags):,}", flush=True)


if __name__ == "__main__":
    main()
