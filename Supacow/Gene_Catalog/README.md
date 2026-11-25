# Supacow (cow rumen) gene catalogue

The Urban Soil gene catalogue steps (`../../Urban_Soil/Gene_Catalog/`) as one
Snakemake workflow, run on Hongyu's 2,444 prokaryotic MAGs
(`/work/microbiome/share/for_kruthi_from_hongyu/supacow_smORFs/MAGs`).

Outputs: `/work/microbiome/users/kruthi/intermediate_results/supacow/Gene_Catalog/`

## Run

```bash
bash run_pipeline.sh --dry-run   # what would run
bash run_pipeline.sh             # one queue job, 32 CPUs, 192 GB, 4 days
```

Resubmitting skips the finished steps.

## Steps

| Step | Rule | Same as Urban Soil | Output |
|---|---|---|---|
| 00 | `prodigal` | new: Urban Soil started from Prodigal files made elsewhere. Single-genome mode (`-p single`), one run per MAG | `prodigal/<MAG>/<MAG>_orfs.fna`, `_proteins.faa`, `.gff` |
| 01 | `orf_catalog` | 01 | `SC.ORF.fna.xz`, `SC.ORF.faa.xz`, `SC.ORF.orig.tsv.xz` (ORF, MAG, Prodigal ID, start, end, strand, partial) |
| 02 | `mmseqs_100nt` | 02: `easy-linclust --min-seq-id 1.0 -c 1.0 --cov-mode 0` | (temporary) |
| 03 | `catalog_100nt` | 03 | `SC.100NT.fna.xz` (`>SC.100NT.id ORF`), `SC.100NT.matches.xz` (ORF, =/C, 100NT) |
| 04 | `mmseqs_95nt` | 04: `easy-cluster --min-seq-id 0.95 -c 0.9 --cov-mode 0` on the 100NT representatives | (temporary) |
| 05 | `catalog_95nt` | 05 | `SC.95NT.fna.xz` (`>SC.95NT.id 100NT`), `SC.95NT.matches.tsv.xz` (100NT, =/R, 95NT) |
| 06 | `cluster_table` | 06 | `SC.clusters.tsv.xz` (ORF, =/C, 100NT, =/R, 95NT) |
| 07a | `rep_proteins` | (part of 08) | `SC.95NT.faa.xz`: protein of each 95NT representative |
| 07 | `eggnog` | 07 + 08 in one: eggNOG-mapper 2.1.12 (diamond) run **once on the 95NT representatives**, not per sample, so no transfer step is needed | `SC.95NT.emapper.annotations.gz` |
| 08 | `summary` | new | `SC.catalogue_summary.tsv`: genes per step, eggNOG and COG coverage |

## Differences from Urban Soil

- **No hard-coded expected counts.** Each step checks that every ORF (or 100NT
  sequence) of the step before appears exactly once.
- **One ID style**: `SC.ORF.000_000_000`, `SC.100NT.000_000_000`,
  `SC.95NT.000_000_000` (Urban Soil mixed grouped and plain 9-digit IDs).
- **eggNOG 2.1.12**, as Urban Soil, so the two catalogues are comparable.
  Hongyu's own eggNOG (emapper 3.0 beta, database 7.0.0) was run on her own
  gene calls, whose files were not shared, so it cannot be linked to these
  genes. The database is copied to the node's `$TMPDIR` before use.
- MMseqs2 is the static AVX2 build from Urban Soil (the conda build crashed with
  illegal instructions on some nodes).

## Not included

The ~109K smORFs from ~1,300 protist genomes (Isotricha, Epidinium, ...) in
`smORF_for_kruthi_annotation.tsv`: those genomes were not shared.
