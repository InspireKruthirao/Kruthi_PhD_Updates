#!/usr/bin/env python3
"""
One script for the smORF neighbourhoods. Builds (or reuses) urban_soil.db,
then writes the neighbourhood of every occurrence of every smORF into ONE
separate database, smorf_neighbourhoods.sqlite (no millions of files):

    occurrences   one row per smORF occurrence: sample, contig, smORF
                  start/end, whether it overlaps a gene or is intergenic
    neighbours    one row per gene around it (WINDOW per side): rank, side,
                  ORF, Original_ID, coordinates, strand, distance to the
                  target gene in bp, COG ID and category, has eggNOG
    smorfs        one row per smORF in the format sweep_db.py reads
    files         one row per occurrence: the text of its neighbors.tsv,
                  EggNOG_annotations.tsv and <contig>.gff (zlib-compressed;
                  read them with the export command, or see read_files)

Genes are matched to eggNOG by sample + Original_ID (the eggNOG query ID),
as in Create_UrbanSoil_GFFs.py; the full eggNOG lines stay in urban_soil.db.

    python3 smorf_neighbourhoods_sqlite.py build       build urban_soil.db (skipped if complete)
    python3 smorf_neighbourhoods_sqlite.py check       eggNOG match rate per sample
    python3 smorf_neighbourhoods_sqlite.py chunk I N   piece I of N (0...N-1) -> parts/chunk_I_of_N.sqlite
    python3 smorf_neighbourhoods_sqlite.py merge N     the N pieces -> smorf_neighbourhoods.sqlite
    python3 smorf_neighbourhoods_sqlite.py export SMORF_ID [out_dir]   write the 3 files of each occurrence

On the HPC, run_smorf_neighbourhoods.sh does it all: a build job, then the
N pieces as one PBS array job (all at the same time), then the merge job.
"""

import gzip
import json
import lzma
import sqlite3
import sys
import time
import zlib
from bisect import bisect_left
from pathlib import Path
from urllib.parse import quote

SMORF_FILE = Path("/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil/UrbanSoil_SMORF_resource/UrbanSoil_100AA_SMORFs_origins.tsv.gz")
GENE_CATALOG = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Gene_Catalog_Urban_Soil/UrbanSoil.ORF.orig.tsv.xz")
EGGNOG_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/egg_nog/results")
DB_PATH = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL/urban_soil.db")
OUT_DIR = Path("/work/microbiome/users/kruthi/intermediate_results/urban_soil/Neighbourhood_Analysis_SQL")
OUT_DB = OUT_DIR / "smorf_neighbourhoods.sqlite"

MIN_OCCURRENCES = 6  # smORFs found in fewer samples are left out (as in _ge6)
BIN_WIDTH = 5
WINDOW = 5           # genes kept on each side of the smORF
BATCH = 100_000
PROGRESS_EVERY = 1000


# ================================================================ database ==

def connect(read_only=False):
    """read_only connections can be opened by many bin jobs at once"""
    if read_only:
        return sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True, timeout=120)
    connection = sqlite3.connect(str(DB_PATH), timeout=120)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-2097152")
    return connection


def complete():
    """Is there a database whose build finished?"""
    if not DB_PATH.is_file():
        return False
    try:
        connection = connect(read_only=True)
        row = connection.execute("SELECT completed FROM build_status LIMIT 1").fetchone()
        connection.close()
        return bool(row and row[0])
    except sqlite3.Error:
        return False


def insert_stream(connection, sql, rows, label):
    batch, total = [], 0
    for row in rows:
        batch.append(row)
        if len(batch) == BATCH:
            connection.executemany(sql, batch)
            connection.commit()
            total += len(batch)
            batch.clear()
            if total % 1_000_000 == 0:
                print(f"  {label}: {total:,} rows", flush=True)
    if batch:
        connection.executemany(sql, batch)
        connection.commit()
        total += len(batch)
    print(f"  {label}: {total:,} rows total", flush=True)


def smorf_rows():
    """(smorf_id, sample, base_contig, start, end) from the smORF origins file"""
    with gzip.open(SMORF_FILE, "rt") as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                start, end = map(int, parts[3].split("-", 1))
            except ValueError:
                continue
            yield parts[0], parts[1], parts[2].split("_polypolish", 1)[0], start, end


def gene_rows():
    """(orf, sample, base_contig, original_id, start, end, strand, partial)
    from the gene catalogue; original_id is the eggNOG query ID"""
    with lzma.open(GENE_CATALOG, "rt") as handle:
        header = next(handle).rstrip("\n").split("\t")
        col = {name: header.index(name) for name in
               ("ORF", "Sample", "Original_ID", "Start", "End", "Strand", "Partial")}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            try:
                original = parts[col["Original_ID"]].strip()
                yield (parts[col["ORF"]].strip(), parts[col["Sample"]].strip(),
                       original.split("_polypolish", 1)[0], original,
                       int(parts[col["Start"]]), int(parts[col["End"]]),
                       parts[col["Strand"]], parts[col["Partial"]])
            except (ValueError, IndexError):
                continue


def eggnog_rows(files, header_box):
    """(sample, query_id, annotation_line); the sample is the folder name"""
    for number, path in enumerate(files, 1):
        sample = path.parent.name
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#query"):
                    if not header_box:
                        header_box.append(line.rstrip("\n"))
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                yield sample, line.split("\t", 1)[0].strip(), line.rstrip("\n")
        print(f"  eggNOG file {number}/{len(files)}: {sample}", flush=True)


def build():
    """Build the database, unless a complete one already exists."""
    if complete():
        print(f"Complete database already exists: {DB_PATH}", flush=True)
        return
    if DB_PATH.exists():
        sys.exit(f"Incomplete database exists: {DB_PATH}; move it aside first")
    files = sorted(EGGNOG_DIR.glob("*/*.emapper.annotations.gz"))
    for path in (SMORF_FILE, GENE_CATALOG):
        if not path.is_file():
            sys.exit(f"Missing input: {path}")
    if not files:
        sys.exit(f"No eggNOG files under {EGGNOG_DIR}")

    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    connection = connect()
    connection.executescript("""
        CREATE TABLE smorf_occurrences (smorf_id TEXT NOT NULL, sample TEXT NOT NULL,
            base_contig TEXT NOT NULL, smorf_start INTEGER NOT NULL, smorf_end INTEGER NOT NULL);
        CREATE TABLE genes (orf TEXT NOT NULL, sample TEXT NOT NULL, base_contig TEXT NOT NULL,
            annotation_query TEXT NOT NULL, gene_start INTEGER NOT NULL, gene_end INTEGER NOT NULL,
            strand TEXT, partial TEXT);
        CREATE TABLE eggnog (sample TEXT NOT NULL, query_id TEXT NOT NULL,
            annotation_line TEXT NOT NULL, PRIMARY KEY (sample, query_id)) WITHOUT ROWID;
        CREATE TABLE eggnog_metadata (annotation_header TEXT NOT NULL);
        CREATE TABLE build_status (completed INTEGER NOT NULL);
    """)
    print("Loading smORF occurrences...", flush=True)
    insert_stream(connection, "INSERT INTO smorf_occurrences VALUES (?, ?, ?, ?, ?)", smorf_rows(), "smORFs")
    connection.execute("CREATE TABLE smorf_summary AS SELECT smorf_id, COUNT(*) AS occurrence_count "
                       "FROM smorf_occurrences GROUP BY smorf_id")
    connection.execute("CREATE UNIQUE INDEX idx_summary ON smorf_summary(smorf_id)")
    connection.execute("CREATE INDEX idx_occurrence ON smorf_occurrences(smorf_id)")
    print("Loading genes...", flush=True)
    insert_stream(connection, "INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?, ?)", gene_rows(), "genes")
    connection.execute("CREATE INDEX idx_gene_location ON genes(sample, base_contig, gene_start, gene_end)")
    print("Loading eggNOG...", flush=True)
    header = []
    insert_stream(connection, "INSERT OR REPLACE INTO eggnog VALUES (?, ?, ?)", eggnog_rows(files, header), "eggNOG")
    connection.execute("INSERT INTO eggnog_metadata VALUES (?)", (header[0] if header else "#query",))
    connection.execute("INSERT INTO build_status VALUES (1)")
    connection.commit()
    connection.close()
    print("Database ready. Run 'check' to see the eggNOG match rate.", flush=True)


def check():
    """eggNOG match rate: for every sample, how many genes have an eggNOG
    annotation (same sample, query ID = Original_ID). A low rate means the
    sample names or the IDs do not agree between the two sources."""
    connection = connect(read_only=True)
    total = matched = 0
    print("sample\tgenes\twith_eggNOG\tpercent", flush=True)
    for (sample,) in connection.execute("SELECT DISTINCT sample FROM eggnog").fetchall():
        n, m = connection.execute(
            "SELECT COUNT(*), SUM(EXISTS (SELECT 1 FROM eggnog e WHERE e.sample = g.sample "
            "AND e.query_id = g.annotation_query)) FROM genes g WHERE g.sample = ?", (sample,)).fetchone()
        m = m or 0
        total, matched = total + n, matched + m
        print(f"{sample}\t{n}\t{m}\t{100 * m / n if n else 0:.1f}", flush=True)
    print(f"ALL\t{total}\t{matched}\t{100 * matched / total if total else 0:.1f}", flush=True)


# ========================================================= neighbourhoods ==

def strand(value):
    value = str(value).strip()
    if value in {"1", "+", "1.0"}:
        return "+"
    if value in {"-1", "-", "-1.0"}:
        return "-"
    return "."


def cog_name(annotation_line):
    """COG1234-E from an eggNOG annotation line (first COG in eggNOG_OGs,
    with the COG_category), as in Create_UrbanSoil_GFFs.py; Unknown if none"""
    if not annotation_line:
        return "Unknown"
    fields = annotation_line.split("\t")
    if len(fields) < 7:
        return "Unknown"
    category = fields[6].strip()
    for item in fields[4].split(","):
        cog = item.split("@", 1)[0].strip()
        if cog.startswith("COG"):
            return f"{cog}-{category}" if category and category != "-" else cog
    return "Unknown"


def neighbourhood(genes, smorf_start, smorf_end, smorf_id):
    """The anchor gene (most overlap with the smORF, or for an intergenic
    smORF the next gene) and WINDOW genes on each side, as in the original
    _sql_submit.pbs. genes is sorted by start. Returns (genes, anchor index
    within them, note)."""
    overlaps = [max(0, min(g[2], smorf_end) - max(g[1], smorf_start)) for g in genes]
    best = max(overlaps)
    if best > 0:
        anchor = overlaps.index(best)
        note = f"TARGET ({smorf_id}) - overlaps by {best} bp"
    else:
        anchor = min(bisect_left([g[1] for g in genes], smorf_start), len(genes) - 1)
        note = f"TARGET ({smorf_id}) - intergenic"
    left = max(0, anchor - WINDOW)
    return genes[left:anchor + WINDOW + 1], anchor - left, note


def esc(value):
    return quote(str(value), safe="._:-")

def gap(a_start, a_end, b_start, b_end):
    """Base pairs between two genes; negative if they overlap."""
    return max(a_start, b_start) - min(a_end, b_end) - 1


def cog_parts(name):
    """('COG1234', ['E', 'G']) from 'COG1234-EG'; ('', []) if not a COG"""
    head, _, letters = name.partition("-")
    if not head.startswith("COG") or not head[3:].isdigit():
        return "", []
    return head, list(letters) if letters.isalpha() else []


def occurrence_rows(occ_id, genes, anchor, annotations):
    """One row per neighbour gene for the neighbours table, plus the
    occurrence as stored in the smorfs table (the format sweep_db.py reads):
    neighbours as [rank, [COG letters]] outwards from the target, upstream
    side first in the target's direction, and for each the same side,
    same_strand, dist_bp, gap_bp, cog_id as build_all_bins_db.py v2."""
    named = []
    for orf, start, end, s, partial, original in genes:
        line = annotations.get(original)
        named.append((orf, start, end, strand(s), partial, original, cog_name(line), line is not None))
    target = named[anchor]
    order = list(range(len(named)))[::-1] if target[3] == "-" else list(range(len(named)))
    i = order.index(anchor)
    side_of, rank_of, prev_of = {}, {}, {}
    for side, idxs in (("up", order[:i][::-1]), ("down", order[i + 1:])):
        prev = anchor
        for rank, j in enumerate(idxs, 1):
            side_of[j], rank_of[j], prev_of[j] = side, rank, prev
            prev = j

    rows, neighbours, info = [], [], []
    for j, (orf, start, end, s, partial, original, name, has_egg) in enumerate(named):
        cog_id, cats = cog_parts(name)
        is_target = j == anchor
        dist = 0 if is_target else gap(start, end, target[1], target[2])
        p = named[prev_of[j]] if not is_target else None
        gap_bp = None if is_target else gap(start, end, p[1], p[2])
        rows.append((occ_id, j, rank_of.get(j, 0), side_of.get(j, "target"), orf, original, start, end, s,
                     partial, int(is_target), dist, gap_bp, cog_id, "".join(cats), int(has_egg)))
    for side in ("up", "down"):
        for j in sorted((j for j in side_of if side_of[j] == side), key=lambda j: rank_of[j]):
            cog_id, cats = cog_parts(named[j][6])
            cog_id = cog_id if cats else ""  # as build_all_bins_db.py: a COG without category counts as none
            neighbours.append([rank_of[j], cats])
            p = named[prev_of[j]]
            info.append([side, int(named[j][3] == target[3]), gap(named[j][1], named[j][2], target[1], target[2]),
                         gap(named[j][1], named[j][2], p[1], p[2]), cog_id])
    occurrence = {"neighbours": neighbours, "neighbour_info": info,
                  "target_categories": cog_parts(target[6])[1]}
    return rows, occurrence


def occurrence_files(sample, contig, smorf_id, smorf_start, smorf_end,
                     genes, anchor, note, annotations, header):
    """The text of neighbors.tsv, EggNOG_annotations.tsv (one line per gene
    of neighbors.tsv, same order, "-" if no eggNOG) and <contig>.gff"""
    neighbors = ["ORF\tSample\tContig\tOriginal_ID\tStart\tEnd\tStrand\tPartial\tNote\n"]
    eggnog = ["ORF\t" + header.lstrip("#") + "\n"]
    gff = ["##gff-version 3\n"]
    for i, (orf, start, end, s, partial, original) in enumerate(genes):
        neighbors.append(f"{orf}\t{sample}\t{contig}\t{original}\t{start}\t{end}\t{strand(s)}\t"
                         f"{partial}\t{note if i == anchor else ''}\n")
        eggnog.append(f"{orf}\t{annotations.get(original) or '-'}\n")
        name = cog_name(annotations.get(original))
        attrs = [f"ID={esc(orf)}", f"Name={esc(name)}", f"target={int(i == anchor)}",
                 f"Original_ID={esc(original)}"]
        if i == anchor:
            attrs += ["Note=TARGET_smORF", "colour=1", f"smorf_id={esc(smorf_id)}",
                      f"smorf_start={smorf_start}", f"smorf_end={smorf_end}",
                      f"anchor={'overlap' if 'overlaps' in note else 'intergenic'}"]
        gff.append("\t".join([contig, "Prodigal_v2.6.3", "CDS", str(start), str(end),
                              "1" if name != "Unknown" else "0.0", strand(s), "0",
                              ";".join(attrs)]) + "\n")
    return "".join(neighbors), "".join(eggnog), "".join(gff)


def pack(text):
    return zlib.compress(text.encode(), 6)


def read_files(db_path, smorf_id):
    """{occ_id: (gff_name, neighbors_tsv, eggnog_annotations_tsv, gff)} for
    one smORF, as text. Example:
        files = read_files("smorf_neighbourhoods.sqlite", "GMSC10.SMORF.000123")"""
    db = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    rows = db.execute("SELECT f.occ_id, f.gff_name, f.neighbors_tsv, f.eggnog_annotations_tsv, f.gff "
                      "FROM files f JOIN occurrences o ON o.occ_id = f.occ_id WHERE o.smorf_id = ? "
                      "ORDER BY f.occ_id", (smorf_id,)).fetchall()
    return {r[0]: (r[1], *(zlib.decompress(b).decode() for b in r[2:])) for r in rows}


SCHEMA = """
    CREATE TABLE IF NOT EXISTS occurrences (
        occ_id TEXT PRIMARY KEY, smorf_id TEXT NOT NULL, bin_name TEXT NOT NULL,
        sample TEXT NOT NULL, contig TEXT NOT NULL, smorf_start INTEGER NOT NULL,
        smorf_end INTEGER NOT NULL, anchor TEXT NOT NULL, n_genes INTEGER NOT NULL);
    CREATE TABLE IF NOT EXISTS neighbours (
        occ_id TEXT NOT NULL, position INTEGER NOT NULL, rank INTEGER NOT NULL, side TEXT NOT NULL,
        orf TEXT NOT NULL, original_id TEXT NOT NULL, gene_start INTEGER NOT NULL,
        gene_end INTEGER NOT NULL, strand TEXT NOT NULL, partial TEXT, is_target INTEGER NOT NULL,
        dist_bp INTEGER NOT NULL, gap_bp INTEGER, cog_id TEXT NOT NULL, cog_category TEXT NOT NULL,
        has_eggnog INTEGER NOT NULL, PRIMARY KEY (occ_id, position));
    CREATE TABLE IF NOT EXISTS smorfs (
        bin_name TEXT NOT NULL, smorf_id TEXT NOT NULL, n_occurrence_dirs INTEGER NOT NULL,
        occurrences_json TEXT NOT NULL, PRIMARY KEY (bin_name, smorf_id));
    CREATE TABLE IF NOT EXISTS files (
        occ_id TEXT PRIMARY KEY, gff_name TEXT NOT NULL, neighbors_tsv BLOB NOT NULL,
        eggnog_annotations_tsv BLOB NOT NULL, gff BLOB NOT NULL);
"""


def bin_name(count):
    """The bin of a smORF found in count samples, e.g. SmORF_neighbourhoods_6_10"""
    low = 1 + ((count - 1) // BIN_WIDTH) * BIN_WIDTH
    return f"SmORF_neighbourhoods_{max(low, MIN_OCCURRENCES)}_{low + BIN_WIDTH - 1}"


def run_chunk(i, n):
    """Piece i of n: every n-th smORF (found in >= MIN_OCCURRENCES samples),
    starting at the i-th, so the n pieces take about the same time.
    Writes OUT_DIR/parts/chunk_<i>_of_<n>.sqlite."""
    connection = connect(read_only=True)
    smorfs = connection.execute(
        "SELECT smorf_id, occurrence_count FROM smorf_summary WHERE occurrence_count >= ? "
        "ORDER BY smorf_id", (MIN_OCCURRENCES,)).fetchall()[i::n]
    run_part(connection, [(s, bin_name(c)) for s, c in smorfs], f"chunk_{i}_of_{n}")


def run_part(connection, smorfs, label):
    """Write the neighbourhoods of the given (smorf_id, bin_name) into
    OUT_DIR/parts/<label>.sqlite"""
    header = connection.execute("SELECT annotation_header FROM eggnog_metadata LIMIT 1").fetchone()[0]
    parts = OUT_DIR / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    part = parts / f"{label}.sqlite"
    tmp = part.with_suffix(".sqlite.part")
    tmp.unlink(missing_ok=True)
    out = sqlite3.connect(str(tmp))
    out.execute("PRAGMA journal_mode=OFF")
    out.execute("PRAGMA synchronous=OFF")
    out.executescript(SCHEMA)

    print(f"[{label}] {len(smorfs):,} smORFs", flush=True)
    t0 = time.time()
    totals = {"occurrences": 0, "written": 0, "no_genes": 0, "genes": 0, "with_eggnog": 0}

    for number, (smorf_id, bin_name) in enumerate(smorfs, 1):
        occurrences = connection.execute(
            "SELECT sample, base_contig, smorf_start, smorf_end FROM smorf_occurrences "
            "WHERE smorf_id = ? ORDER BY sample, base_contig, smorf_start", (smorf_id,)).fetchall()
        occ_rows, gene_rows_, occ_json, file_rows = [], [], [], []
        for sample, contig, smorf_start, smorf_end in occurrences:
            genes = connection.execute(
                "SELECT orf, gene_start, gene_end, strand, partial, annotation_query FROM genes "
                "WHERE sample = ? AND base_contig = ? ORDER BY gene_start", (sample, contig)).fetchall()
            if not genes:
                totals["no_genes"] += 1
                continue
            genes, anchor, note = neighbourhood(genes, smorf_start, smorf_end, smorf_id)
            ids = [g[5] for g in genes]
            annotations = dict(connection.execute(
                f"SELECT query_id, annotation_line FROM eggnog WHERE sample = ? AND query_id IN "
                f"({','.join('?' * len(ids))})", (sample, *ids)))
            occ_id = f"{smorf_id}|{sample}|{contig}|{smorf_start}"
            rows, occurrence = occurrence_rows(occ_id, genes, anchor, annotations)
            occ_rows.append((occ_id, smorf_id, bin_name, sample, contig, smorf_start, smorf_end,
                             "overlap" if "overlaps" in note else "intergenic", len(genes)))
            gene_rows_.extend(rows)
            occ_json.append(occurrence)
            texts = occurrence_files(sample, contig, smorf_id, smorf_start, smorf_end,
                                     genes, anchor, note, annotations, header)
            file_rows.append((occ_id, f"{contig}.gff", *map(pack, texts)))
            totals["genes"] += len(ids)
            totals["with_eggnog"] += sum(i in annotations for i in ids)
        totals["occurrences"] += len(occurrences)
        totals["written"] += len(occ_rows)
        out.executemany("INSERT OR REPLACE INTO occurrences VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", occ_rows)
        out.executemany("INSERT OR REPLACE INTO neighbours VALUES "
                        "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", gene_rows_)
        out.executemany("INSERT OR REPLACE INTO files VALUES (?, ?, ?, ?, ?)", file_rows)
        out.execute("INSERT OR REPLACE INTO smorfs VALUES (?, ?, ?, ?)",
                    (bin_name, smorf_id, len(occ_rows), json.dumps(occ_json)))
        if number % PROGRESS_EVERY == 0 or number == len(smorfs):
            out.commit()
            print(f"[{label}] {number:,}/{len(smorfs):,} smORFs ({time.time() - t0:.0f}s)", flush=True)

    out.commit()
    out.close()
    tmp.replace(part)  # only a finished bin gets its final name
    pct = 100 * totals["with_eggnog"] / totals["genes"] if totals["genes"] else 0
    print(f"[{label}] DONE: {totals['written']:,} of {totals['occurrences']:,} occurrences written, "
          f"{totals['no_genes']:,} without genes; {totals['with_eggnog']:,} of {totals['genes']:,} "
          f"genes have eggNOG ({pct:.1f}%)", flush=True)


def merge(n):
    """Combine the n pieces into OUT_DB."""
    paths = [OUT_DIR / "parts" / f"chunk_{i}_of_{n}.sqlite" for i in range(n)]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        sys.exit(f"{len(missing)} piece(s) not finished yet: {missing[:10]}")
    OUT_DB.unlink(missing_ok=True)
    out = sqlite3.connect(str(OUT_DB))
    out.execute("PRAGMA journal_mode=OFF")
    out.executescript(SCHEMA)
    for number, path in enumerate(paths, 1):
        out.execute("ATTACH DATABASE ? AS part", (str(path),))
        for table in ("occurrences", "neighbours", "smorfs", "files"):
            out.execute(f"INSERT INTO {table} SELECT * FROM part.{table}")
        out.commit()
        out.execute("DETACH DATABASE part")
        print(f"  merged {number}/{n}", flush=True)
    out.execute("CREATE INDEX idx_occ_smorf ON occurrences(smorf_id)")
    out.commit()
    n_smorfs, n_occ, n_genes, n_egg = out.execute(
        "SELECT (SELECT COUNT(*) FROM smorfs), (SELECT COUNT(*) FROM occurrences), "
        "COUNT(*), SUM(has_eggnog) FROM neighbours").fetchone()
    out.close()
    print(f"\nDone: {OUT_DB} ({OUT_DB.stat().st_size / 1e6:.0f} MB): {n_smorfs:,} smORFs, "
          f"{n_occ:,} occurrences, {n_genes:,} genes, {n_egg:,} with eggNOG "
          f"({100 * n_egg / n_genes if n_genes else 0:.1f}%)")


def export(smorf_id, out_dir):
    """Write the stored neighbors.tsv, EggNOG_annotations.tsv and
    <contig>.gff of every occurrence of smorf_id into
    out_dir/<smorf_id>/<sample>_<contig>_<smORF start>/"""
    db = sqlite3.connect(f"file:{OUT_DB}?mode=ro", uri=True)
    occs = dict(((o, (s, c, st)) for o, s, c, st in db.execute(
        "SELECT occ_id, sample, contig, smorf_start FROM occurrences WHERE smorf_id = ?", (smorf_id,))))
    if not occs:
        sys.exit(f"{smorf_id} is not in {OUT_DB}")
    for occ_id, (gff_name, neighbors, eggnog, gff) in read_files(OUT_DB, smorf_id).items():
        sample, contig, start = occs[occ_id]
        folder = Path(out_dir) / smorf_id / f"{sample}_{contig}_{start}"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "neighbors.tsv").write_text(neighbors)
        (folder / "EggNOG_annotations.tsv").write_text(eggnog)
        (folder / gff_name).write_text(gff)
    print(f"Wrote the files of {len(occs)} occurrence(s) under {Path(out_dir) / smorf_id}")


def main():
    command = sys.argv[1] if len(sys.argv) > 1 else ""
    if command == "build":
        build()
    elif command == "check":
        check()
    elif command == "chunk" and len(sys.argv) == 4:
        run_chunk(int(sys.argv[2]), int(sys.argv[3]))
    elif command == "merge" and len(sys.argv) == 3:
        merge(int(sys.argv[2]))
    elif command == "export" and len(sys.argv) in (3, 4):
        export(sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else "export")
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
