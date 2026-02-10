#!/usr/bin/env python3
import gzip
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

SMORF_FILE   = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/smORFCatalog/100AA_SmORFs_origins.tsv.gz"
GENE_CATALOG = "/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz"
EGGNOG_DIR   = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs"

BASE_OUT = "/work/microbiome/users/kruthi"
OUT_DIR  = os.path.join(BASE_OUT, "SmORF_neighbourhoods_26_30_no_min_overlap")

WINDOW = 5
PROGRESS_EVERY = 25

def parse_coords(s):
    a, b = s.split("-", 1)
    return int(a), int(b)

def base_contig_name(contig):
    return contig.split("_polypolish")[0]

def load_gene_index():
    gene_df = pd.read_csv(
        GENE_CATALOG,
        sep="\t",
        compression="xz",
        header=0,
        dtype=str,
        low_memory=False
    )

    gene_df = gene_df[["ORF", "Sample", "Original_ID", "Start", "End", "Strand", "Partial"]].copy()
    gene_df = gene_df.rename(columns={"Original_ID": "Contig"})

    gene_df["Start"] = pd.to_numeric(gene_df["Start"], errors="coerce")
    gene_df["End"]   = pd.to_numeric(gene_df["End"], errors="coerce")
    gene_df = gene_df.dropna(subset=["Start", "End"])
    gene_df["Start"] = gene_df["Start"].astype(int)
    gene_df["End"]   = gene_df["End"].astype(int)

    gene_df["BaseContig"] = gene_df["Contig"].str.split("_polypolish").str[0]
    gene_df = gene_df.sort_values(["Sample", "BaseContig", "Start"]).reset_index(drop=True)

    index = {}
    for (sample, bc), sub in gene_df.groupby(["Sample", "BaseContig"], sort=False):
        index[(sample, bc)] = sub.reset_index(drop=True)
    return index

def get_neighbours(gene_sub_df, smorf_start, smorf_end, smorf_id):
    if gene_sub_df is None or gene_sub_df.empty:
        return None

    df = gene_sub_df

    # Find the gene with the largest overlap (can be 0)
    overlap = (np.minimum(df["End"].values, smorf_end) - np.maximum(df["Start"].values, smorf_start))
    overlap = np.clip(overlap, 0, None)

    if overlap.size == 0:
        return None

    best_idx = int(overlap.argmax())  # even if max overlap is 0, take the closest/best

    # Take WINDOW genes left and right of the best match
    left = max(0, best_idx - WINDOW)
    right = min(len(df), best_idx + WINDOW + 1)

    neigh = df.iloc[left:right].copy()
    neigh["Note"] = ""
    neigh.loc[best_idx, "Note"] = f"TARGET ({smorf_id}) - overlap {int(overlap[best_idx])} bp"
    return neigh

def load_eggnog_cache_for_sample(sample):
    annot_file = os.path.join(EGGNOG_DIR, sample, f"{sample}.emapper.annotations")
    if not os.path.exists(annot_file):
        return None

    contig_to_lines = defaultdict(list)
    with open(annot_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            contig = line.split("\t", 1)[0]
            contig_to_lines[contig].append(line)
    return contig_to_lines

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Pass 1: counting SmORF occurrences...")
    counts = Counter()
    with gzip.open(SMORF_FILE, "rt") as f:
        next(f)
        for line in f:
            smorf_id = line.split("\t", 1)[0]
            counts[smorf_id] += 1

    selected = {sid for sid, c in counts.items() if 26 <= c <= 30}
    selected_list = sorted(selected)
    total_smorfs = len(selected_list)

    print(f"Selected SmORFs with counts in [26,30]: {total_smorfs}")

    print("Loading and indexing gene catalog (one-time)...")
    gene_index = load_gene_index()
    print(f"Gene index keys: {len(gene_index)}")

    print("Pass 2: collecting occurrences for selected SmORFs...")
    smorf_rows = defaultdict(list)
    with gzip.open(SMORF_FILE, "rt") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            smorf_id = parts[0]
            if smorf_id not in selected:
                continue
            sample = parts[1]
            contig = parts[2]
            coords = parts[3]
            strand = parts[4] if len(parts) > 4 else ""
            smorf_rows[smorf_id].append((sample, contig, coords, strand))

    summary_path = os.path.join(OUT_DIR, "summary.tsv")
    with open(summary_path, "w") as s:
        s.write("SmORF_ID\tTotalOccurrences\tOccurrencesFound\tOccurrencesProcessed\tNeighborsWritten\tSkippedNoGenes\tSkippedNoOverlapPossible\tEggNOGMissing\n")

    eggnog_cache = {}
    done = 0
    total_occ_processed = 0
    total_neighbors_written = 0

    for smorf_id in selected_list:
        occs = smorf_rows.get(smorf_id, [])
        total_occ = counts[smorf_id]
        occ_found = len(occs)

        processed = 0
        written = 0
        skipped_nogenes = 0
        skipped_no_overlap_possible = 0
        eggnog_missing = 0

        for (sample, contig, coords, strand) in occs:
            processed += 1
            total_occ_processed += 1

            bc = base_contig_name(contig)
            gene_sub = gene_index.get((sample, bc))
            if gene_sub is None or gene_sub.empty:
                skipped_nogenes += 1
                continue

            start, end = parse_coords(coords)
            neigh = get_neighbours(gene_sub, start, end, smorf_id)
            if neigh is None or neigh.empty:
                skipped_no_overlap_possible += 1
                continue

            out_dir = os.path.join(OUT_DIR, smorf_id, f"{sample}_{bc}")
            os.makedirs(out_dir, exist_ok=True)

            neigh[["ORF", "Sample", "Contig", "Start", "End", "Strand", "Partial", "Note"]].to_csv(
                os.path.join(out_dir, "neighbors.tsv"),
                sep="\t",
                index=False
            )

            if sample not in eggnog_cache:
                eggnog_cache[sample] = load_eggnog_cache_for_sample(sample)

            contig_map = eggnog_cache.get(sample)
            egg_out = os.path.join(out_dir, "EggNOG_annotations.tsv")
            with open(egg_out, "w") as out:
                out.write("#query\tORF\tSample\tContig\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
                if contig_map is None:
                    eggnog_missing += 1
                else:
                    for _, r in neigh.iterrows():
                        orf = r["ORF"]
                        q = r["Contig"]
                        lines = contig_map.get(q, [])
                        if lines:
                            # Write the first (usually best) annotation line for simplicity
                            out.write(f"{smorf_id}\t{orf}\t{sample}\t{q}\t{lines[0]}")

            written += 1
            total_neighbors_written += 1

        with open(summary_path, "a") as s:
            s.write(
                f"{smorf_id}\t{total_occ}\t{occ_found}\t{processed}\t{written}\t{skipped_nogenes}\t{skipped_no_overlap_possible}\t{eggnog_missing}\n"
            )

        done += 1
        if done % PROGRESS_EVERY == 0 or done == total_smorfs:
            print(
                f"Progress: {done}/{total_smorfs} SmORFs done | "
                f"occurrences processed={total_occ_processed} | neighborhoods written={total_neighbors_written}"
            )

    print("Done.")
    print(f"Output folder: {OUT_DIR}")
    print(f"Summary: {summary_path}")

if __name__ == "__main__":
    main()
