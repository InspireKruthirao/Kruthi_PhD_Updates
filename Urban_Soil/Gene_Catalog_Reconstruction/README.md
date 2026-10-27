# Gene catalogue reconstruction

Goal: rebuild the Urban Soil gene catalogue (`../Gene_Catalog/`) in much less
time, and reuse the pipeline for larger data (SPIRE, Aim 2). Each new step is
checked against the existing catalogue before replacing the old one.

## Plan

| Step | Now | Idea |
|---|---|---|
| Compress Prodigal outputs | xz, 8 h per sample | zstd (multi-threaded) |
| 100% clustering | MMseqs2, 100 CPUs, 700 GB, up to 48 h | identical sequences = group by sequence hash |
| 100NT / 95NT post-processing, cluster table | Python on .xz text, up to 48 h | polars + parquet (or SQLite) |
| eggNOG | every ORF, per sample | cluster representatives only, copied to members |
| 95% clustering | MMseqs2, 100 CPUs, 1.2 TB | keep; tune settings |

## To do

- Real run times of the current steps (`resources_used.walltime` in the job logs)
- Ask Luis whether SPIRE already provides gene predictions and annotations
- Rebuild one step at a time, compare with the existing outputs
