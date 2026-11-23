#!/usr/bin/env python3
"""Compare two Urban Soil gene catalogues (old Sep 2026 run vs Gene_Catalog_v2).

IDs are renumbered between runs, so the comparison is done on what the IDs
stand for: DNA sequences (64-bit hashes), the original ORF IDs, and the way
ORFs are grouped into clusters.

  compare_catalogs.py parse KIND SIDE PATH WORKDIR   # one file -> parquet of ints/hashes
  compare_catalogs.py report WORKDIR OUTDIR          # compare old vs new, write report.md

KIND: fasta100 | fasta95 | clusters | matches100 | matches95 | emapper
"""
import array
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
import polars as pl

EMAPPER_TEXT_COLS = ["seed_ortholog", "eggNOG_OGs", "max_annot_lvl", "COG_category",
                     "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
                     "KEGG_Pathway", "CAZy", "PFAMs"]


# ---------------------------------------------------------------- helpers
def open_stream(path):
    cmd = ["xz", "-T2", "-dc", path] if path.endswith(".xz") else ["pigz", "-p", "2", "-dc", path]
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 22)


def close_stream(p):
    p.stdout.close()
    if p.wait() != 0:
        sys.exit(f"decompression failed for {p.args}")


def h64(b):
    return int.from_bytes(hashlib.blake2b(b, digest_size=8).digest(), "little")


def idnum(b):
    """b'US.ORF.002_801_681' -> 2801681 ; b'US.100NT.010558031' -> 10558031"""
    return int(b.rsplit(b".", 1)[1].replace(b"_", b""))


def fmt(prefix, n):
    if n is None or n < 0:
        return "NA"
    return f"US.{prefix}.{n // 1_000_000:03d}_{n // 1000 % 1000:03d}_{n % 1000:03d}"


def save(cols, kind, side, workdir, stats):
    df = pl.DataFrame({k: np.frombuffer(v, dtype=d) for k, (v, d) in cols.items()})
    df.write_parquet(os.path.join(workdir, f"{kind}.{side}.parquet"))
    stats["rows"] = df.height
    with open(os.path.join(workdir, f"{kind}.{side}.json"), "w") as f:
        json.dump(stats, f)


# ---------------------------------------------------------------- parsers
def parse_fasta(path, with_orf):
    ids, orfs, hs, lens = array.array("q"), array.array("q"), array.array("Q"), array.array("I")
    p = open_stream(path)
    cur, orf, chunks = None, -1, []

    def flush():
        seq = b"".join(chunks)
        ids.append(cur); orfs.append(orf); hs.append(h64(seq)); lens.append(len(seq))

    for line in p.stdout:
        if line[:1] == b">":
            if cur is not None:
                flush()
            parts = line[1:].split()
            cur = idnum(parts[0])
            orf = idnum(parts[1]) if with_orf and len(parts) > 1 else -1
            chunks = []
        else:
            chunks.append(line.rstrip(b"\r\n"))
    if cur is not None:
        flush()
    close_stream(p)
    cols = {"id": (ids, np.int64), "seqh": (hs, np.uint64), "len": (lens, np.uint32)}
    if with_orf:
        cols["rep_orf"] = (orfs, np.int64)
    return cols, {}


def parse_tsv_ids(path, ncols_expected):
    """clusters: ORF f1 100NT f2 95NT ; matches100: ORF f 100NT ; matches95: 100NT f 95NT"""
    n_id = (ncols_expected + 1) // 2
    idcols = [array.array("q") for _ in range(n_id)]
    flagcols = [array.array("B") for _ in range(n_id - 1)]
    bad = 0
    p = open_stream(path)
    for line in p.stdout:
        c = line.rstrip(b"\r\n").split(b"\t")
        if len(c) != ncols_expected:
            bad += 1
            continue
        try:
            vals = [idnum(c[2 * i]) for i in range(n_id)]
        except ValueError:      # header line
            bad += 1
            continue
        for i in range(n_id):
            idcols[i].append(vals[i])
        for i in range(n_id - 1):
            flagcols[i].append(c[2 * i + 1][0] if c[2 * i + 1] else 0)
    close_stream(p)
    return idcols, flagcols, {"skipped_lines": bad}


def parse_emapper(path):
    p = open_stream(path)
    idx, ids, ev, sc = None, array.array("q"), array.array("d"), array.array("d")
    hcols = {c: array.array("Q") for c in EMAPPER_TEXT_COLS}
    skipped = 0
    for line in p.stdout:
        c = line.rstrip(b"\r\n").split(b"\t")
        if line[:1] == b"#" or c[0] == b"US.95NT" or (len(c) > 1 and c[1].startswith(b"#")):
            if idx is None and b"seed_ortholog" in c:
                names = [x.decode().lstrip("#") for x in c]
                idx = {n: names.index(n) for n in EMAPPER_TEXT_COLS + ["evalue", "score"]}
            continue
        if idx is None:
            sys.exit("emapper header not found before data")
        try:
            ids.append(idnum(c[0]))
            ev.append(float(c[idx["evalue"]])); sc.append(float(c[idx["score"]]))
        except (ValueError, IndexError):
            skipped += 1
            continue
        for n in EMAPPER_TEXT_COLS:
            hcols[n].append(h64(c[idx[n]]))
    close_stream(p)
    cols = {"id": (ids, np.int64), "evalue": (ev, np.float64), "score": (sc, np.float64)}
    cols.update({n: (a, np.uint64) for n, a in hcols.items()})
    return cols, {"skipped_lines": skipped}


def parse(kind, side, path, workdir):
    t0 = time.time()
    if kind in ("fasta100", "fasta95"):
        cols, st = parse_fasta(path, with_orf=(kind == "fasta100"))
    elif kind == "emapper":
        cols, st = parse_emapper(path)
    else:
        spec = {"clusters": (5, ["orf", "nt100", "nt95"]),
                "matches100": (3, ["orf", "nt100"]),
                "matches95": (3, ["nt100", "nt95"])}[kind]
        idcols, flagcols, st = parse_tsv_ids(path, spec[0])
        cols = {n: (a, np.int64) for n, a in zip(spec[1], idcols)}
        cols.update({f"flag{i + 1}": (a, np.uint8) for i, a in enumerate(flagcols)})
    st.update(path=path, seconds=round(time.time() - t0))
    save(cols, kind, side, workdir, st)
    print(f"{kind} {side}: done in {st['seconds']} s", flush=True)


# ---------------------------------------------------------------- report
class Report:
    def __init__(self, outdir):
        self.outdir, self.lines = outdir, []

    def h(self, t):
        self.lines.append(f"\n## {t}\n")

    def p(self, t=""):
        self.lines.append(t)

    def table(self, rows, header):
        self.p("| " + " | ".join(header) + " |")
        self.p("|" + "---|" * len(header))
        for r in rows:
            self.p("| " + " | ".join(f"{x:,}" if isinstance(x, int) else str(x) for x in r) + " |")
        self.p()

    def examples(self, name, df):
        path = os.path.join(self.outdir, f"examples_{name}.tsv")
        df.head(50).write_csv(path, separator="\t")
        self.p(f"Examples: `{path}`\n")

    def write(self):
        with open(os.path.join(self.outdir, "report.md"), "w") as f:
            f.write("\n".join(self.lines) + "\n")


def pct(a, b):
    return f"{100 * a / b:.4f}%" if b else "NA"


def load(workdir, kind, side):
    return pl.read_parquet(os.path.join(workdir, f"{kind}.{side}.parquet"))


def set_overlap(R, label, old, new):
    o, n = old.unique(), new.unique()
    both = o.to_frame("k").join(n.to_frame("k"), on="k", how="inner").height
    R.table([[label, o.len(), n.len(), both, o.len() - both, n.len() - both]],
            ["what", "unique in old", "unique in new", "in both", "only old", "only new"])
    return both


def partition(R, label, df, a, b):
    """df: one row per ORF with old group a and new group b."""
    pair = df.group_by([a, b]).agg(pl.len().alias("n"))
    pair = (pair.join(df.group_by(a).agg(pl.len().alias("na")), on=a)
                .join(df.group_by(b).agg(pl.len().alias("nb")), on=b))
    exact = pair.filter((pl.col("n") == pl.col("na")) & (pl.col("n") == pl.col("nb")))
    n_old, n_new = df[a].n_unique(), df[b].n_unique()
    orfs_same = int(exact["n"].sum())
    R.table([[label, df.height, n_old, n_new, exact.height, orfs_same, pct(orfs_same, df.height)]],
            ["level", "ORFs compared", "clusters old", "clusters new",
             "identical clusters", "ORFs in identical clusters", "% ORFs"])
    return pair.filter(~((pl.col("n") == pl.col("na")) & (pl.col("n") == pl.col("nb"))))


def report(workdir, outdir):
    R = Report(outdir)
    R.p("# Urban Soil gene catalogue: old (Sep 2026) vs new (Gene_Catalog_v2)")
    R.p(f"\nGenerated {time.strftime('%Y-%m-%d %H:%M')}. Sequences compared by 64-bit hash "
        "(chance of any false match is negligible, ~1 in 1000 over the whole catalogue).")
    stats = {}
    for fn in sorted(os.listdir(workdir)):
        if fn.endswith(".json"):
            stats[fn[:-5]] = json.load(open(os.path.join(workdir, fn)))

    # ---- 1. ORF files (checksums from run.sh)
    R.h("1. ORF files (input to clustering)")
    rows = []
    for f in ["UrbanSoil.ORF.fna.xz", "UrbanSoil.ORF.orig.tsv.xz"]:
        r = {}
        for side in ("old", "new"):
            md5, lines = open(os.path.join(workdir, f"md5.{f}.{side}")).read().split()
            r[side] = (md5, int(lines))
        rows.append([f, r["old"][1], r["new"][1],
                     "IDENTICAL" if r["old"][0] == r["new"][0] else "DIFFERENT"])
    R.table(rows, ["file", "lines old", "lines new", "decompressed content"])

    # ---- 2. 100NT sequences
    R.h("2. 100NT catalogue (US.100NT.fna.xz)")
    f100 = {s: load(workdir, "fasta100", s) for s in ("old", "new")}
    R.table([[s, d.height, d["seqh"].n_unique(), int(d["len"].sum())] for s, d in f100.items()],
            ["side", "sequences", "unique sequences", "total bp"])
    set_overlap(R, "DNA sequences", f100["old"]["seqh"], f100["new"]["seqh"])
    set_overlap(R, "representative ORF IDs", f100["old"]["rep_orf"], f100["new"]["rep_orf"])
    j = f100["old"].select("rep_orf", "seqh").join(f100["new"].select("rep_orf", "seqh"),
                                                    on="rep_orf", how="inner", suffix="_new")
    R.p(f"Same representative ORF in both: {j.height:,}; of these, same sequence: "
        f"{(j['seqh'] == j['seqh_new']).sum():,} ({pct(int((j['seqh'] == j['seqh_new']).sum()), j.height)})\n")
    only_old = f100["old"].join(f100["new"], on="seqh", how="anti")
    R.examples("100NT_seq_only_in_old", only_old.select(
        pl.col("id").map_elements(lambda x: fmt("100NT", x), return_dtype=pl.Utf8).alias("old_100NT"),
        pl.col("rep_orf").map_elements(lambda x: fmt("ORF", x), return_dtype=pl.Utf8).alias("rep_orf"), "len").head(50))

    # ---- 3. 95NT sequences
    R.h("3. 95NT catalogue (US.95NT.fna.xz)")
    f95 = {s: load(workdir, "fasta95", s) for s in ("old", "new")}
    R.table([[s, d.height, d["seqh"].n_unique(), int(d["len"].sum())] for s, d in f95.items()],
            ["side", "sequences", "unique sequences", "total bp"])
    set_overlap(R, "DNA sequences", f95["old"]["seqh"], f95["new"]["seqh"])

    # ---- 4. clusters: ORF -> 100NT -> 95NT
    R.h("4. Cluster membership (US.clusters.tsv.xz)")
    cl = {}
    for s in ("old", "new"):
        d = load(workdir, "clusters", s)
        d = (d.join(f100[s].select(pl.col("id").alias("nt100"), pl.col("rep_orf").alias("rep100")), on="nt100", how="left")
              .join(f95[s].select(pl.col("id").alias("nt95"), pl.col("seqh").alias("rep95_seqh")), on="nt95", how="left"))
        cl[s] = d
        R.p(f"- {s}: {d.height:,} rows ({stats[f'clusters.{s}']['skipped_lines']} lines skipped), "
            f"{d['orf'].n_unique():,} unique ORFs; 100NT IDs not found in 100NT fasta: {d['rep100'].null_count():,}; "
            f"95NT IDs not found in 95NT fasta: {d['rep95_seqh'].null_count():,}")
    R.p()
    j = cl["old"].join(cl["new"], on="orf", how="full", suffix="_new", coalesce=True)
    both = j.filter(pl.col("nt100").is_not_null() & pl.col("nt100_new").is_not_null())
    R.p(f"ORFs in both: {both.height:,}; only old: {j['nt100_new'].null_count():,}; "
        f"only new: {j['nt100'].null_count():,}\n")
    R.p("**Same grouping?** A cluster counts as identical when exactly the same ORFs are in it in both runs.\n")
    diff100 = partition(R, "100NT", both, "nt100", "nt100_new")
    diff95 = partition(R, "95NT", both, "nt95", "nt95_new")
    s100 = int((both["rep100"] == both["rep100_new"]).sum())
    s95 = int((both["rep95_seqh"] == both["rep95_seqh_new"]).sum())
    R.table([["100NT representative (same ORF ID)", s100, pct(s100, both.height)],
             ["95NT representative (same sequence)", s95, pct(s95, both.height)]],
            ["per ORF", "ORFs with same", "%"])
    fl = []
    for c in ("flag1", "flag2"):
        for s in ("old", "new"):
            vc = cl[s][c].value_counts().sort(c)
            fl.append([c, s, ", ".join(f"'{chr(k) if k else ''}': {v:,}" for k, v in vc.iter_rows())])
    R.table(fl, ["column", "side", "values"])
    R.examples("100NT_differing_clusters", diff100.sort("n", descending=True).with_columns(
        pl.col("nt100").map_elements(lambda x: fmt("100NT", x), return_dtype=pl.Utf8),
        pl.col("nt100_new").map_elements(lambda x: fmt("100NT", x), return_dtype=pl.Utf8)).rename(
        {"n": "ORFs_shared", "na": "ORFs_in_old_cluster", "nb": "ORFs_in_new_cluster"}))
    R.examples("95NT_differing_clusters", diff95.sort("n", descending=True).with_columns(
        pl.col("nt95").map_elements(lambda x: fmt("95NT", x), return_dtype=pl.Utf8),
        pl.col("nt95_new").map_elements(lambda x: fmt("95NT", x), return_dtype=pl.Utf8)).rename(
        {"n": "ORFs_shared", "na": "ORFs_in_old_cluster", "nb": "ORFs_in_new_cluster"}))
    del j, both, cl, diff100, diff95

    # ---- 5. matches files
    R.h("5. Matches files")
    m = {s: load(workdir, "matches100", s) for s in ("old", "new")}
    R.p("US.100NT.matches.xz (ORF -> 100NT):\n")
    set_overlap(R, "ORF IDs listed", m["old"]["orf"], m["new"]["orf"])
    for s in m:
        m[s] = m[s].join(f100[s].select(pl.col("id").alias("nt100"), pl.col("rep_orf")), on="nt100", how="left")
    j = m["old"].join(m["new"], on="orf", how="inner", suffix="_new")
    same = int((j["rep_orf"] == j["rep_orf_new"]).sum())
    R.p(f"ORFs in both pointing to the same representative ORF: {same:,} of {j.height:,} ({pct(same, j.height)})\n")

    m = {s: load(workdir, "matches95", s) for s in ("old", "new")}
    R.p("US.95NT.matches.tsv.xz (100NT -> 95NT):\n")
    for s in m:
        m[s] = (m[s].join(f100[s].select(pl.col("id").alias("nt100"), pl.col("rep_orf")), on="nt100", how="left")
                    .join(f95[s].select(pl.col("id").alias("nt95"), pl.col("seqh")), on="nt95", how="left"))
        R.p(f"- {s}: {m[s].height:,} rows; 100NT IDs not found in 100NT fasta: {m[s]['rep_orf'].null_count():,}; "
            f"95NT IDs not found in 95NT fasta: {m[s]['seqh'].null_count():,}")
    R.p()
    j = m["old"].drop_nulls("rep_orf").join(m["new"].drop_nulls("rep_orf"), on="rep_orf", how="inner", suffix="_new")
    same = int((j["seqh"] == j["seqh_new"]).sum())
    R.p(f"Rows matched by representative ORF: {j.height:,}; same 95NT sequence: {same:,} ({pct(same, j.height)})\n")
    del m, j

    # ---- 6. eggNOG annotations
    R.h("6. eggNOG annotations (US.95NT.emapper.annotations.gz)")
    em = {}
    for s in ("old", "new"):
        d = load(workdir, "emapper", s).join(
            f95[s].select(pl.col("id"), pl.col("seqh")), on="id", how="left")
        R.p(f"- {s}: {d.height:,} annotated genes ({stats[f'emapper.{s}']['skipped_lines']} lines skipped); "
            f"IDs not found in 95NT fasta: {d['seqh'].null_count():,}")
        em[s] = d.drop_nulls("seqh").unique("seqh", keep="first")
    R.p()
    set_overlap(R, "annotated sequences", em["old"]["seqh"], em["new"]["seqh"])
    j = em["old"].join(em["new"], on="seqh", how="inner", suffix="_new")
    rows = []
    for c in EMAPPER_TEXT_COLS + ["evalue", "score"]:
        same = int((j[c] == j[c + "_new"]).sum())
        rows.append([c, same, j.height - same, pct(same, j.height)])
    d = (j["score_new"] - j["score"]).abs()
    rows.append(["score within 1 bit", int((d <= 1).sum()), int((d > 1).sum()), pct(int((d <= 1).sum()), j.height)])
    R.table(rows, ["column", "same", "different", "% same"])
    diff = j.filter(pl.col("seed_ortholog") != pl.col("seed_ortholog_new")).select(
        pl.col("id").map_elements(lambda x: fmt("95NT", x), return_dtype=pl.Utf8).alias("old_95NT"),
        pl.col("id_new").map_elements(lambda x: fmt("95NT", x), return_dtype=pl.Utf8).alias("new_95NT"),
        "score", "score_new")
    R.examples("emapper_different_seed_ortholog", diff)

    R.h("Run statistics")
    R.table([[k, v.get("rows", ""), v.get("seconds", ""), v.get("path", "")] for k, v in stats.items()],
            ["parse step", "rows", "seconds", "file"])
    R.write()
    print("\n".join(R.lines))


if __name__ == "__main__":
    if sys.argv[1] == "parse":
        parse(*sys.argv[2:6])
    elif sys.argv[1] == "report":
        report(*sys.argv[2:4])
    else:
        sys.exit(__doc__)
