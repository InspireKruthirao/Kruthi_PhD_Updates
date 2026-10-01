#!/usr/bin/env python3
"""
Build the single all-bins SQLite database -- but do the slow part (walking
every bin's millions of small files) as separate PBS jobs running AT THE
SAME TIME, instead of one script walking every bin one after another.

Three ways to run this same file:

    python3 build_all_bins_db.py
        submits one PBS job per bin that doesn't have a <bin_name>.jsonl
        yet -- every bin gets consolidated in parallel instead of in turn

    python3 build_all_bins_db.py <bin_dir>
        (this is what each PBS job runs -- you don't normally call this
        yourself) walks ONE bin, writes its <bin_name>.jsonl

    python3 build_all_bins_db.py --merge [out_db]
        once qstat shows the jobs are done, run this ONCE -- no PBS
        needed, it's fast, it just reads the .jsonl files -- to combine
        everything into one database: <parent_dir>/smorfs_all_bins.sqlite
"""

import json
import re
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

# ---- config: same values as Functional_Framework_Analysis.py ----
PARENT_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL_ge6")
SCRIPT_PATH = Path(__file__).resolve()
LOG_DIR = PARENT_DIR / "logs_consolidate"
PYTHON_BIN = "/home/n12228516/.conda/envs/kruthi/bin/python3"

NAME_RE = re.compile(r"^COG\d+-([A-Za-z]+)$")


# ---------------------------------------------------------- parsing/walking --

def parse_gff(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) != 9:
                continue
            try:
                start = int(fields[3])
            except ValueError:
                print(f"Warning: bad coordinates in {path}; skipping line")
                continue
            attrs = dict(
                kv.strip().split("=", 1) for kv in fields[8].split(";") if "=" in kv
            )
            genes.append({
                "start": start,
                "strand": fields[6],
                "name": attrs.get("Name"),
                "target": attrs.get("target") == "1",
                "note": attrs.get("Note", ""),
            })
    return sorted(genes, key=lambda gene: gene["start"])


def get_cog_categories(name):
    match = NAME_RE.match(name or "")
    return list(match.group(1)) if match else []


def find_target(genes):
    for gene in genes:
        if gene["note"] == "TARGET_smORF":
            return gene
    return next((gene for gene in genes if gene["target"]), None)


def neighbour_genes_all(genes, target):
    ordered = genes[::-1] if target["strand"] == "-" else genes
    i = ordered.index(target)
    for side in (ordered[:i][::-1], ordered[i + 1:]):
        yield from enumerate(side, start=1)


def consolidate_bin(bin_dir):
    """Walk one bin folder once, yield one dict per smORF."""
    bin_dir = Path(bin_dir)
    smorf_dirs = sorted(p for p in bin_dir.iterdir() if p.is_dir())
    total = len(smorf_dirs)

    for number, smorf_dir in enumerate(smorf_dirs, 1):
        occurrence_dirs = [p for p in smorf_dir.iterdir() if p.is_dir()]
        occurrences = []

        for occ_dir in occurrence_dirs:
            gff_files = sorted(occ_dir.glob("*.gff"))
            if not gff_files:
                continue
            if len(gff_files) > 1:
                raise ValueError(
                    f"Expected exactly one .gff file in {occ_dir}, found "
                    f"{len(gff_files)}: {[p.name for p in gff_files]}"
                )

            genes = parse_gff(gff_files[0])
            target = find_target(genes)
            if target is None:
                continue

            neighbours = [
                [rank, get_cog_categories(gene["name"])]
                for rank, gene in neighbour_genes_all(genes, target)
            ]
            occurrences.append({
                "neighbours": neighbours,
                "target_categories": get_cog_categories(target["name"]),
            })

        yield {
            "smorf_id": smorf_dir.name,
            "n_occurrence_dirs": len(occurrence_dirs),
            "occurrences": occurrences,
        }

        if number % 500 == 0 or number == total:
            print(f"[{bin_dir.name}] {number:,}/{total:,} smORFs consolidated", flush=True)


# --------------------------------- mode 2: one bin -> its own .jsonl file --

def write_jsonl_for_bin(bin_dir):
    bin_dir = Path(bin_dir)
    out_file = bin_dir.parent / f"{bin_dir.name}.jsonl"
    t0 = time.time()
    n = 0
    with out_file.open("w") as fh:
        for record in consolidate_bin(bin_dir):
            fh.write(json.dumps(record) + "\n")
            n += 1
    print(f"\n[{bin_dir.name}] wrote {n:,} smORFs to {out_file.name} in {time.time() - t0:.1f}s", flush=True)


# ----------------------------- mode 1: one PBS job per bin, all parallel --

def submit_all():
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    bins = sorted(p for p in PARENT_DIR.glob("SmORF_neighbourhoods_*") if p.is_dir())
    if not bins:
        print(f"No SmORF_neighbourhoods_* folders found under {PARENT_DIR}")
        sys.exit(1)

    todo = [b for b in bins if not (PARENT_DIR / f"{b.name}.jsonl").exists()]
    skipped = len(bins) - len(todo)
    if skipped:
        print(f"{skipped} bin(s) already have a .jsonl -- skipping those.")
    if not todo:
        print("Every bin already has a .jsonl. Run with --merge to build the database.")
        return

    print(f"Submitting {len(todo)} job(s), one per bin, all running at the same time...")
    for bin_path in todo:
        bin_name = bin_path.name
        job_script = f"""#!/bin/bash
#PBS -l select=1:ncpus=2:mem=64gb
#PBS -l walltime=24:00:00
#PBS -q cpu_batch_exec

set -eo pipefail
export PYTHONNOUSERSITE=1

{PYTHON_BIN} -s -u "{SCRIPT_PATH}" "{bin_path}"
"""
        result = subprocess.run(
            ["qsub", "-N", f"jsonl_{bin_name}", "-o", str(LOG_DIR / f"{bin_name}.log"), "-j", "oe"],
            input=job_script, text=True, capture_output=True,
        )
        job_id = result.stdout.strip()
        print(f"  {bin_name} -> {job_id}")

    print("\nAll jobs submitted. Watch progress with:")
    print(f"  tail -f {LOG_DIR}/*.log")
    print("Check status with:")
    print("  qstat -u $USER")
    print(f"\nOnce qstat shows they're all done, run:")
    print(f"  python3 {SCRIPT_PATH.name} --merge")


# --------------------------- mode 3: combine every bin's data into 1 db --

def merge_to_sqlite(out_db=None):
    out_db = Path(out_db) if out_db else PARENT_DIR / "smorfs_all_bins.sqlite"

    bins = sorted(p for p in PARENT_DIR.glob("SmORF_neighbourhoods_*") if p.is_dir())
    if not bins:
        print(f"No SmORF_neighbourhoods_* folders found under {PARENT_DIR}")
        sys.exit(1)

    conn = sqlite3.connect(out_db)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS smorfs (
            bin_name TEXT NOT NULL,
            smorf_id TEXT NOT NULL,
            n_occurrence_dirs INTEGER NOT NULL,
            occurrences_json TEXT NOT NULL,
            PRIMARY KEY (bin_name, smorf_id)
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_bin ON smorfs(bin_name)")

    for bin_dir in bins:
        bin_name = bin_dir.name
        existing = conn.execute(
            "SELECT COUNT(*) FROM smorfs WHERE bin_name = ?", (bin_name,)
        ).fetchone()[0]
        if existing:
            print(f"[{bin_name}] already in {out_db.name} ({existing:,} smORFs) -- skipping")
            continue

        jsonl_path = PARENT_DIR / f"{bin_name}.jsonl"
        t0 = time.time()
        rows = []

        if jsonl_path.exists():
            print(f"[{bin_name}] reading {jsonl_path.name} ...", flush=True)
            with jsonl_path.open() as fh:
                for line in fh:
                    rec = json.loads(line)
                    rows.append((bin_name, rec["smorf_id"], rec["n_occurrence_dirs"],
                                  json.dumps(rec["occurrences"])))
        else:
            print(f"[{bin_name}] no .jsonl yet -- walking directly (one-time cost) ...", flush=True)
            for rec in consolidate_bin(bin_dir):
                rows.append((bin_name, rec["smorf_id"], rec["n_occurrence_dirs"],
                              json.dumps(rec["occurrences"])))

        conn.executemany(
            "INSERT INTO smorfs (bin_name, smorf_id, n_occurrence_dirs, occurrences_json) "
            "VALUES (?, ?, ?, ?)",
            rows,
        )
        conn.commit()
        print(f"[{bin_name}] {len(rows):,} smORFs added in {time.time() - t0:.1f}s")

    total = conn.execute("SELECT COUNT(*) FROM smorfs").fetchone()[0]
    conn.close()
    size_mb = out_db.stat().st_size / 1_000_000
    print(f"\nDone: {out_db} holds {total:,} smORFs across {len(bins)} bins ({size_mb:.1f} MB)")


def main():
    if len(sys.argv) == 1:
        submit_all()
        return

    if sys.argv[1] == "--merge":
        out_db = sys.argv[2] if len(sys.argv) > 2 else None
        merge_to_sqlite(out_db)
        return

    bin_dir = Path(sys.argv[1])
    if not bin_dir.is_dir():
        print("Usage:")
        print("  python3 build_all_bins_db.py                 # submit one PBS job per bin")
        print("  python3 build_all_bins_db.py <bin_dir>        # (run by each job) consolidate one bin")
        print("  python3 build_all_bins_db.py --merge [out_db] # combine everything into one sqlite")
        sys.exit(1)
    write_jsonl_for_bin(bin_dir)


if __name__ == "__main__":
    main()
