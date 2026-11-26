# Urban Soil Gene Catalogue

This directory is the reconstructed, cleaned-up **framework for the
work already done** on building the Urban Soil non-redundant gene
catalogue: starting from per-sample Prodigal ORF predictions and
ending with a 100NT catalogue, a 95NT catalogue, an ORF-to-cluster
mapping table, and eggNOG functional annotations for the 95NT
representatives. These are exactly the inputs that the downstream
smORF / gene-neighbourhood analyses (see `../Neighbourhood_Analysis*`)
consume.

It was assembled from the actual scripts/PBS jobs run on the HPC
(most of them recovered verbatim from shell history; a couple of
`vim`-authored scripts were reconstructed from their invocations and
logs -- see "Provenance & caveats" below).

## Pipeline (run in numeric order)

```
Prodigal (external, per sample)
   <sample>/<sample>_orfs.fna, <sample>_proteins.faa
        │
        ▼
01  Build_urban_soil_orf_catalog.py(.sh)
        concatenate all samples, assign stable US.ORF.### IDs
        →  UrbanSoil.ORF.fna.xz
        →  UrbanSoil.ORF.orig.tsv.xz   (ORF → sample/original ID/coords)
        │
        ▼
02  run_mmseqs_100NT_splitmem.sh
        mmseqs easy-linclust --min-seq-id 1.0 -c 1.0   (exact dedup)
        →  US.mmseqs.100NT.splitmem_{rep_seq.fasta,cluster.tsv}
        │
        ▼
03  Run_US_100NT_catalog.py(.pbs)
        renumber representatives → US.100NT.###_###_###
        →  US.100NT.fna.xz
        →  US.100NT.matches.xz        (ORF → 100NT, cols: orf, =/C, rep100)
        │
        ▼
04  run_mmseqs_95nt_64bit.sh
        mmseqs easy-cluster --min-seq-id 0.95 -c 0.9   (cluster the
        100NT REPRESENTATIVES, not the raw ORFs)
        →  mmseqs_95nt_64bit/US.mmseqs.95NT_{rep_seq.fasta,cluster.tsv}
        │
        ▼
05  postprocess_mmseqs95.py(.pbs)
        renumber representatives → US.95NT.#########
        →  US.95NT.fna.xz
        →  US.95NT.matches.tsv.xz     (100NT → 95NT, cols: rep100, =/R, rep95)
        │
        ▼
06  make_US_cluster_table.py(.pbs)
        join 03's and 05's mapping files
        →  US.clusters.tsv.xz  (ORF, =/C, rep100, =/R, rep95)  <- the
           "cluster-mapping" deliverable
        │
        ▼
07  Run_Eggnog_Mapper.sh            (58 samples, one eggNOG-mapper run each)
        emapper.py -m diamond --itype proteins  (per-sample protein FASTA)
        →  egg_nog/results/<sample>/<sample>.emapper.annotations
        │
        ▼
08  US95_emapper.py(.pbs)
        for each 95NT cluster, pick its exact canonical ORF (the one
        that is both its own 100NT rep and its own 95NT rep) and copy
        that ORF's per-sample eggNOG row onto the 95NT ID
        →  US.95NT.emapper.annotations.gz   <- the "eggNOG annotation"
           deliverable
```

### Final deliverables (in `Gene_Catalog_Urban_Soil/`)

| File | Produced by | Description |
|---|---|---|
| `UrbanSoil.ORF.fna.xz` / `.orig.tsv.xz` | 01 | all ORFs + original sample/coords |
| `US.100NT.fna.xz` / `.matches.xz` | 03 | exact-dedup ("100NT") catalogue |
| `US.95NT.fna.xz` / `.matches.tsv.xz` | 05 | 95%-identity ("95NT") catalogue |
| `US.clusters.tsv.xz` | 06 | ORF → 100NT → 95NT cluster-mapping table |
| `US.95NT.emapper.annotations.gz` | 08 | eggNOG functional annotation of 95NT reps |

## Key design decision: hierarchical clustering

95NT clustering is run on the **100NT representative sequences**,
not on the raw ~265M ORFs directly (this is what makes the ~256.8M →
~199.2M reduction computationally tractable). `US.clusters.tsv.xz`
(Stage 06) is therefore a proper two-level join: every original ORF
maps to exactly one 100NT representative, and every 100NT
representative maps to exactly one 95NT representative.

## ID formats

Catalogue IDs use two slightly different zero-padded conventions
depending on which script assigned them:

- `US.100NT.###_###_###` (underscore-grouped, from `Run_US_100NT_catalog.py`'s
  own local `pad9()`)
- `US.95NT.#########` (plain 9-digit, from the shared `lib.pad9()`)

`06_make_US_cluster_table.py`'s `numeric_id()` helper strips
underscores before parsing, so it transparently accepts either style
-- this is intentional, not a bug, but worth knowing if you ever
touch ID formatting.

## Shared dependency: `lib.py`

`05_postprocess_mmseqs95.py` and `06_make_US_cluster_table.py` import
`pad9()`/`xz_out()` from a small shared helper module that lives in
the sibling **Shanghai Dogs** project:

```
/work/microbiome/shanghai_dogs/resource_generation/lib.py
```

A vendored copy is kept in this directory (`lib.py`) purely for
reference/fallback -- swap the `sys.path.insert(...)` + `from lib
import ...` lines for a local import if you ever need to run these
scripts somewhere that path isn't mounted.

## Provenance & caveats

Everything above other than Prodigal itself was reconstructed by
walking the interactive shell history of the sessions that actually
built this catalogue on the HPC (Aug-Sep). Two scripts are
**reconstructions**, not verbatim recoveries, because they were
authored directly with `vim` rather than via a logged heredoc, so
their exact original text isn't in shell history:

- `04_run_mmseqs_95nt_64bit.sh` -- parameters reconstructed from the
  surrounding commands/logs and cross-checked against what
  `05_postprocess_mmseqs95.py` actually consumed. It also documents
  *why* a statically-linked AVX2 mmseqs2 binary was needed (the
  conda build crashed with illegal instructions on some compute
  nodes -- likely an AVX-512-only build vs. an older CPU).
- `07_Run_Eggnog_Mapper.sh` -- reconstructed as a PBS job array over
  58 samples, using the exact `emapper.py` CLI flags that were
  tested and confirmed working in history (`-m diamond --itype
  proteins --data_dir ... --override`).

Everything else (`01`, `02`, `03`, `05`'s Python, `06`, `08`, and the
matching `.pbs` wrappers) was recovered byte-for-byte from the
heredocs (`cat > file <<'EOF' ... EOF`) that originally created them,
and each `.py` file has been verified to at least `py_compile`
cleanly.

**Not included / out of scope:** the Prodigal ORF prediction step
itself (per the task description, this pipeline starts *from*
Prodigal-predicted ORFs) and the downstream eggNOG database setup
(v2.1.12, `diamond` search mode, database copied to
`egg_nog/databases`) are treated as pre-existing prerequisites and
were not reconstructed here.

Before re-running any of this on the real data: update the hardcoded
`/work/microbiome/users/kruthi/...` paths if your layout differs,
and double-check the `EXPECTED_*` sanity-check constants in each
script (`EXPECTED_REPRESENTATIVES`, `N_100NT`, `EXPECTED_95NT`,
`EXPECTED_ORFS`, etc.) against your actual data -- they were sized
for this specific 58-sample Urban Soil run (~265.4M ORFs → ~256.8M
100NT reps → ~199.2M 95NT reps) and will need updating for a
different sample set.
