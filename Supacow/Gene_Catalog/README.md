# Supacow (cow rumen) gene catalogue and smORFs

The Urban Soil gene catalogue steps (`../../Urban_Soil/Gene_Catalog/`) as one
Snakemake workflow, run on Hongyu's 2,444 prokaryotic MAGs
(`/work/microbiome/share/for_kruthi_from_hongyu/supacow_smORFs/MAGs`).

Outputs: `/work/microbiome/users/kruthi/intermediate_results/supacow/Gene_Catalog/`

## Run

```bash
bash run_pipeline.sh --dry-run   # what would run
bash run_pipeline.sh             # one queue job, 32 CPUs, 400 GB, 4 days
```

Resubmitting skips the finished steps (and finished eggNOG / GMSC-mapper chunks).

## Steps

| Step | Rule | Same as Urban Soil | Output |
|---|---|---|---|
| 00 | `prodigal` | new: Urban Soil started from Prodigal files made elsewhere. Single-genome mode (`-p single`), one run per MAG | `prodigal/<MAG>/<MAG>_orfs.fna`, `_proteins.faa`, `.gff` |
| 01 | `orf_catalog` | 01 | `COW.ORF.fna.xz`, `COW.ORF.faa.xz`, `COW.ORF.orig.tsv.xz` (ORF, MAG, Prodigal ID, start, end, strand, partial) |
| 02 | `mmseqs_100nt` | 02: `easy-linclust --min-seq-id 1.0 -c 1.0 --cov-mode 0` | (temporary) |
| 03 | `catalog_100nt` | 03 | `COW.100NT.fna.xz` (`>COW.100NT.id ORF`), `COW.100NT.matches.xz` (ORF, =/C, 100NT) |
| 04 | `mmseqs_95nt` | 04: `easy-cluster --min-seq-id 0.95 -c 0.9 --cov-mode 0` on the 100NT representatives | (temporary) |
| 05 | `catalog_95nt` | 05 | `COW.95NT.fna.xz` (`>COW.95NT.id 100NT`), `COW.95NT.matches.tsv.xz` (100NT, =/R, 95NT) |
| 06 | `cluster_table` | 06 | `COW.clusters.tsv.xz` (ORF, =/C, 100NT, =/R, 95NT) |
| 07a | `rep_proteins` | (part of 08) | `COW.95NT.faa.xz`: protein of each 95NT representative |
| 07 | `eggnog_split`, `eggnog_chunk`, `eggnog_merge` | 07 + 08 in one: eggNOG-mapper 2.1.12 (diamond) run **once on the 95NT representatives**, not per sample, so no transfer step is needed. In 20 chunks of 8 CPUs, so finished chunks survive a time-out | `COW.95NT.emapper.annotations.gz` |
| 09 | `gmsc_input`, `gmsc_mapper`, `smorf_resource` | `Urban_Soil/SmORF_Catalog/Run_GMSC_MAPPER.sh` + `Create_resource_files.py`: GMSC-mapper 0.2.0 (`--tool diamond`) on the MAGs, in 16 groups of MAGs (contigs renamed `<MAG>\|<contig>`; GMSC-mapper finds genes contig by contig, so grouping does not change the result). Unique mapped smORFs, most common first, get pad9 IDs | `smorfs/COW_100AA_SMORFs_sequences.faa.gz`, `_origins.tsv.gz` (SmORF ID, Sample ID = MAG, Contig, Coordinates, Strand), `_habitat_taxonomy.tsv`, `_merge_log.txt` |
| 10 | `smorf_cdhit`, `smorf_90aa` | `Cluster_SmORFs.sh` + `Process_90AA_clusters.py`: CD-HIT `-c 0.9 -n 5` | `smorfs/COW_SMORF_Clusters.tsv.gz` (100AA, 90AA, similarity), `smorfs/COW_90AA_SMORFs.faa.gz` |
| 08 | `summary` | new | `COW.catalogue_summary.tsv`: genes per step, eggNOG and COG coverage, smORF counts |

## Differences from Urban Soil

- **No hard-coded expected counts.** Each step checks that every ORF (or 100NT
  sequence) of the step before appears exactly once.
- **One ID style (pad9, grouped)**: genes `COW.ORF.000_000_000`, `COW.100NT.000_000_000`,
  `COW.95NT.000_000_000`; smORFs `COW_SM.100AA.000_000_000`, `COW_SM.90AA.000_000_000`
  (as `US_SM.100AA.000_000_000` in Urban Soil).
- **Databases** (eggNOG 2.1.12, GMSC) are copied to the node's local disk once per job
  and shared by the parallel chunks.
- **eggNOG 2.1.12**, as Urban Soil, so the two catalogues are comparable.
  Hongyu's own eggNOG (emapper 3.0 beta, database 7.0.0) was run on her own
  gene calls, whose files were not shared, so it cannot be linked to these
  genes. The database is copied to the node's `$TMPDIR` before use.
- MMseqs2 is the static AVX2 build from Urban Soil (the conda build crashed with
  illegal instructions on some nodes).

## smORFs: ours vs Hongyu's

Hongyu's `smORF_for_kruthi.faa` merges five tools (Bakta, GMSC-mapper, davin,
proSmORF, smORFinder; 1.09M of 1.87M found by GMSC-mapper). This pipeline runs
GMSC-mapper itself, as for Urban Soil, so the cow and Urban Soil smORFs come from
the same tool and settings.

## Not included

The ~109K smORFs from ~1,300 protist genomes (Isotricha, Epidinium, ...) in
`smORF_for_kruthi_annotation.tsv`: those genomes were not shared.
