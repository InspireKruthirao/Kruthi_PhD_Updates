mkdir -p /work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_DuckDB/logs

cat > /work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_DuckDB/UrbanSoil_Overall_Smorfs_analysis_duckdb.py << 'EOFPYTHON'
#!/usr/bin/env python3

import gzip
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb

SMORF_FILE = Path("/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil/UrbanSoil_SMORF_resource/UrbanSoil_100AA_SMORFs_origins.tsv.gz")
GENE_CATALOG = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Gene_Catalog_Urban_Soil/UrbanSoil.ORF.orig.tsv.xz")
EGGNOG_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/egg_nog/results")
BASE_OUT = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_DuckDB")
DB_PATH = str(BASE_OUT / "urban_soil.duckdb")

BIN_WIDTH = 5
WINDOW = 5
PROGRESS_EVERY = 25

def load_data():
    print("Initializing DuckDB database...", flush=True)
    
    conn = duckdb.connect(DB_PATH)
    
    print("Pass 1: Loading smorfs...", flush=True)
    smorf_counts = Counter()
    smorf_data = []
    
    with gzip.open(SMORF_FILE, "rt") as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            
            smorf_id, sample, contig, coordinates = parts[:4]
            base_contig = contig.split("_polypolish", 1)[0]
            smorf_start, smorf_end = map(int, coordinates.split("-", 1))
            
            smorf_counts[smorf_id] += 1
            smorf_data.append({
                "smorf_id": smorf_id,
                "sample": sample,
                "base_contig": base_contig,
                "smorf_start": smorf_start,
                "smorf_end": smorf_end
            })
    
    print(f"Unique SmORFs: {len(smorf_counts)}", flush=True)
    
    # Create table from data
    conn.execute("CREATE TABLE smorf_occurrences AS SELECT * FROM ?", [pd.DataFrame(smorf_data)])
    conn.execute("CREATE INDEX idx_smorf_id ON smorf_occurrences(smorf_id)")
    conn.execute("CREATE INDEX idx_sample_contig ON smorf_occurrences(sample, base_contig)")
    
    # Create summary
    conn.execute("""
        CREATE TABLE smorf_summary AS
        SELECT smorf_id, COUNT(*) as occurrence_count
        FROM smorf_occurrences
        GROUP BY smorf_id
    """)
    
    print("Loading genes...", flush=True)
    
    gene_df = pd.read_csv(GENE_CATALOG, sep="\t", compression="xz", dtype=str, low_memory=False)
    gene_df = gene_df[["ORF", "Sample", "Original_ID", "Start", "End", "Strand", "Partial"]].copy()
    gene_df.columns = ["orf", "sample", "contig", "start", "end", "strand", "partial"]
    gene_df["start"] = pd.to_numeric(gene_df["start"], errors="coerce")
    gene_df["end"] = pd.to_numeric(gene_df["end"], errors="coerce")
    gene_df = gene_df.dropna(subset=["start", "end"])
    gene_df["base_contig"] = gene_df["contig"].str.split("_polypolish", n=1).str[0]
    
    conn.execute("CREATE TABLE genes AS SELECT * FROM ?", [gene_df[["orf", "sample", "base_contig", "start", "end", "strand", "partial"]]])
    conn.execute("CREATE INDEX idx_gene_sample_contig ON genes(sample, base_contig, start, end)")
    
    print("Loading eggNOG files...", flush=True)
    annot_files = list(sorted(EGGNOG_DIR.glob("*/*.emapper.annotations.gz")))
    print(f"Found {len(annot_files)} annotation files", flush=True)
    
    eggnog_data = []
    for idx, annot_file in enumerate(annot_files, 1):
        if idx % 10 == 0 or idx == len(annot_files):
            print(f"  Loading file {idx}/{len(annot_files)}", flush=True)
        
        with gzip.open(annot_file, "rt") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                query_id = line.split("\t", 1)[0]
                eggnog_data.append({"query_id": query_id, "annotation_line": line})
    
    conn.execute("CREATE TABLE eggnog AS SELECT * FROM ?", [pd.DataFrame(eggnog_data)])
    conn.execute("CREATE INDEX idx_eggnog_query ON eggnog(query_id)")
    
    print("DuckDB database ready!", flush=True)
    return conn

def get_neighbours(conn, sample, base_contig, smorf_start, smorf_end, smorf_id):
    result = conn.execute(
        "SELECT orf, start, end, strand, partial FROM genes WHERE sample = ? AND base_contig = ? ORDER BY start",
        [sample, base_contig]
    ).fetchall()
    
    if not result:
        return None
    
    gene_list = [{"ORF": r[0], "Start": r[1], "End": r[2], "Strand": r[3], "Partial": r[4]} for r in result]
    
    overlap = np.array([min(g["End"], smorf_end) - max(g["Start"], smorf_start) for g in gene_list])
    overlap = np.clip(overlap, 0, None)
    
    if overlap.max() > 0:
        anchor_idx = int(overlap.argmax())
        note = f"TARGET ({smorf_id}) - overlaps by {int(overlap.max())} bp"
    else:
        anchor_idx = int(np.searchsorted([g["Start"] for g in gene_list], smorf_start))
        anchor_idx = max(0, min(anchor_idx, len(gene_list) - 1))
        note = f"TARGET ({smorf_id}) - intergenic"
    
    left = max(0, anchor_idx - WINDOW)
    right = min(len(gene_list), anchor_idx + WINDOW + 1)
    
    neighbours = gene_list[left:right]
    for n in neighbours:
        n["Note"] = ""
    neighbours[anchor_idx - left]["Note"] = note
    
    return neighbours

def main():
    BASE_OUT.mkdir(parents=True, exist_ok=True)
    (BASE_OUT / "logs").mkdir(parents=True, exist_ok=True)
    
    if not Path(DB_PATH).exists():
        conn = load_data()
    else:
        print("Database already exists, skipping data load", flush=True)
        conn = duckdb.connect(DB_PATH)
    
    print("Pass 1: Counting SmORF occurrences (from database)...", flush=True)
    
    smorf_data = conn.execute("SELECT smorf_id, occurrence_count FROM smorf_summary").fetchall()
    smorf_data = {row[0]: row[1] for row in smorf_data}
    
    max_count = max(smorf_data.values())
    all_bins = [(start, start + BIN_WIDTH - 1) for start in range(1, max_count + 1, BIN_WIDTH)]
    
    print(f"Unique SmORFs: {len(smorf_data)}", flush=True)
    
    log_path = BASE_OUT / "bin_creation_log.tsv"
    with log_path.open("w") as log:
        log.write("BinStart\tBinEnd\tLabel\tSmORF_Count\tStatus\tFolder\n")
        for start, end in all_bins:
            label = f"{start}_{end}"
            smorfs_in_bin = [s for s, c in smorf_data.items() if start <= c <= end]
            status = "CREATED" if smorfs_in_bin else "SKIPPED_EMPTY"
            folder = f"SmORF_neighbourhoods_{label}" if smorfs_in_bin else "NA"
            log.write(f"{start}\t{end}\t{label}\t{len(smorfs_in_bin)}\t{status}\t{folder}\n")
    
    print(f"Bin log written: {log_path}", flush=True)
    
    for start, end in all_bins:
        label = f"{start}_{end}"
        smorfs_in_bin = sorted([s for s, c in smorf_data.items() if start <= c <= end])
        
        if not smorfs_in_bin:
            continue
        
        output_root = BASE_OUT / f"SmORF_neighbourhoods_{label}"
        output_root.mkdir(parents=True, exist_ok=True)
        
        print(f"\nProcessing bin: {label} ({len(smorfs_in_bin)} smorfs)", flush=True)
        
        summary_path = output_root / "summary.tsv"
        with summary_path.open("w") as summary:
            summary.write("SmORF_ID\tTotalOccurrences\tOccurrencesProcessed\tNeighborsWritten\tSkippedNoGenes\tEggNOGMissing\n")
        
        for i, smorf_id in enumerate(smorfs_in_bin):
            occurrences = conn.execute(
                "SELECT sample, base_contig, smorf_start, smorf_end FROM smorf_occurrences WHERE smorf_id = ?",
                [smorf_id]
            ).fetchall()
            
            processed = written = skipped_no_genes = eggnog_missing = 0
            
            for sample, base_contig, smorf_start, smorf_end in occurrences:
                processed += 1
                neighbours = get_neighbours(conn, sample, base_contig, smorf_start, smorf_end, smorf_id)
                
                if not neighbours:
                    skipped_no_genes += 1
                    continue
                
                output_dir = output_root / smorf_id / f"{sample}_{base_contig}"
                output_dir.mkdir(parents=True, exist_ok=True)
                
                with open(output_dir / "neighbors.tsv", "w") as f:
                    f.write("ORF\tSample\tContig\tStart\tEnd\tStrand\tPartial\tNote\n")
                    for n in neighbours:
                        f.write(f"{n['ORF']}\t{sample}\t{base_contig}\t{n['Start']}\t{n['End']}\t{n['Strand']}\t{n['Partial']}\t{n['Note']}\n")
                
                annot = conn.execute("SELECT annotation_line FROM eggnog WHERE query_id = ?", [base_contig]).fetchone()
                
                with open(output_dir / "EggNOG_annotations.tsv", "w") as f:
                    if annot:
                        f.write(annot[0])
                    else:
                        eggnog_missing += 1
                
                written += 1
            
            with summary_path.open("a") as summary:
                summary.write(f"{smorf_id}\t{len(occurrences)}\t{processed}\t{written}\t{skipped_no_genes}\t{eggnog_missing}\n")
            
            if (i + 1) % PROGRESS_EVERY == 0 or (i + 1) == len(smorfs_in_bin):
                print(f"  [{label}] {i + 1}/{len(smorfs_in_bin)} done", flush=True)
        
        print(f"Completed bin: {label}", flush=True)
    
    conn.close()
    print("\nALL DONE", flush=True)

if __name__ == "__main__":
    main()
EOFPYTHON

chmod +x /work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_DuckDB/UrbanSoil_Overall_Smorfs_analysis_duckdb.py
