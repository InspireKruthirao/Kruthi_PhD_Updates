# Gene catalogue pipeline (Snakemake + PBS)

Assemblies → Prodigal → ORF catalogue → 100NT (MMseqs2 linclust) → 95NT (MMseqs2 cluster) → eggNOG.
One command runs every step as its own PBS job, in order, and resumes where it stopped after a failure.
eggNOG runs per sample at the same time as the MMseqs2 steps.

## Outputs (in `outdir`)
| File | Content |
|---|---|
| `UrbanSoil.ORF.fna.xz` | all ORFs, IDs `US.ORF.000_000_000` |
| `UrbanSoil.ORF.orig.tsv.xz` | ORF → sample, Prodigal ID, start, end, strand, partial |
| `US.100NT.fna.xz` / `US.100NT.matches.xz` | 100% identity catalogue / ORF → 100NT (`=` rep, `C` collapsed) |
| `US.95NT.fna.xz` / `US.95NT.matches.tsv.xz` | 95% catalogue / 100NT → 95NT (`=` rep, `R` reduced) |
| `US.clusters.tsv.xz` | ORF, 100NT, 95NT in one table |
| `US.95NT.emapper.annotations.gz` | eggNOG annotation of each 95NT representative |

Per-sample Prodigal and eggNOG files are kept in `prodigal/` and `eggnog/`. Large intermediate files go to `tmp/` and are deleted automatically.

## Run
Run on the login node inside `tmux` (Snakemake itself is light; every step runs as its own PBS job):
```bash
cd Gene_Catalog_Pipeline
~/.pixi/bin/snakemake --profile pbs --configfile config/urban_soil.yaml -n   # dry run: list the steps
tmux new -s catalog
~/.pixi/bin/snakemake --profile pbs --configfile config/urban_soil.yaml      # real run
# detach: Ctrl+B then D   |   reattach (same login node): tmux attach -t catalog
```
Progress: `qstat -u $USER`. Job logs go to `logs/`. If a step fails, fix the cause and run the same command again: finished steps are skipped.

## New dataset
Copy `config/urban_soil.yaml`, change `name`, `prefix`, `assembly_list`, `assembly_suffix`, `outdir`, and run with `--configfile` pointing to the new file. CPUs, memory and walltime per step are set in the same file.

## HPC storage rule
eggNOG jobs copy the eggNOG database to `${TMPDIR}` and read it only from there, and at most 8 run at once (`eggnog_slots` in `pbs/config.yaml`). Never point eggNOG at the database on `/work` directly.

## Checks
Every step stops with `CHECK FAILED: ...` if something does not add up, for example: ORF counts, every ORF placed in exactly one cluster, each representative present once, matching eggNOG headers, or one annotation per 95NT.
