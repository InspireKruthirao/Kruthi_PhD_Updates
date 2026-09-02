# Urban Soil Gene Catalogue

## Purpose

This workflow builds the **Urban Soil non-redundant gene catalogue** from Prodigal-predicted ORFs and generates 100NT, 95NT, cluster-mapping, and eggNOG annotation files for downstream smORF and gene-neighbourhood analyses.

## Locations

**Scripts**

```text
/work/microbiome/users/kruthi/intermediate_results/Kruthi_PhD_Updates/Urban_Soil/Gene_Catalog/
```

**Intermediate files**

```text
/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/
```

**Final catalogue**

```text
/work/microbiome/users/kruthi/intermediate_results/urban_soil/Gene_Catalog_Urban_Soil/
```

---

## Workflow

```text
Assemblies
   ↓
Prodigal ORF prediction
   ↓
Urban Soil ORF catalogue
   ↓
100NT MMseqs2 clustering
   ↓
100NT catalogue
   ↓
95NT MMseqs2 clustering
   ↓
95NT catalogue
   ↓
ORF → 100NT → 95NT mapping
   ↓
eggNOG annotation
```

---

## 1. ORF Prediction

**Script**

```text
run_prodigal.pbs
```

Runs Prodigal on the Urban Soil assemblies and generates per-sample ORF predictions.

---

## 2. Compress Prodigal Outputs

**Script**

```text
Compress_prodigal_outputs.pbs
```

Compresses the per-sample Prodigal output files.

---

## 3. Build the Urban Soil ORF Catalogue

**Scripts**

```text
Build_urban_soil_orf_catalog.py
Build_urban_soil_orf_catalog.sh
```

Combines all predicted ORFs and assigns standardized IDs:

```text
US.ORF.000_000_000
```

**Outputs**

```text
UrbanSoil.ORF.fna.xz
UrbanSoil.ORF.orig.tsv.xz
```

**Total ORFs**

```text
265,392,705
```

---

## 4. 100NT Clustering

**Script**

```text
run_mmseqs_linclust_splitmem.sh
```

Clusters the complete ORF catalogue with MMseqs2 at:

```text
100% nucleotide identity
100% coverage
```

**Raw files**

```text
Prodigal/mmseqs_test/US.mmseqs.100NT.splitmem_rep_seq.fasta
Prodigal/mmseqs_test/US.mmseqs.100NT.splitmem_cluster.tsv
```

**Result**

```text
Original ORFs          265,392,705
100NT representatives  256,829,048
Collapsed ORFs           8,563,657
```

---

## 5. Build the 100NT Catalogue

**Scripts**

```text
Run_US_100NT_catalog.py
Run_US_100NT_catalog.pbs
```

Converts the raw MMseqs2 output into standardized 100NT IDs:

```text
US.100NT.000_000_000
```

**Outputs**

```text
US.100NT.fna.xz
US.100NT.matches.xz
```

`US.100NT.fna.xz` contains the representative sequences.

`US.100NT.matches.xz` maps each original ORF to its 100NT representative.

---

## 6. 95NT Clustering

**Script**

```text
run_mmseqs_95nt_64bit.sh
```

Clusters the 100NT representatives with MMseqs2 at:

```text
95% nucleotide identity
90% coverage
```

**Raw files**

```text
Prodigal/mmseqs_95nt_64bit/US.mmseqs.95NT_rep_seq.fasta
Prodigal/mmseqs_95nt_64bit/US.mmseqs.95NT_cluster.tsv
```

**Result**

```text
100NT representatives  256,829,048
95NT representatives   199,172,868
Collapsed at 95NT       57,656,180
```

---

## 7. Build the 95NT Catalogue

**Scripts**

```text
postprocess_mmseqs95.py
run_postprocess_mmseqs95.pbs
```

**Outputs**

```text
US.95NT.fna.xz
US.95NT.matches.tsv.xz
```

Maps the 100NT representatives to their corresponding 95NT representatives.

---

## 8. Build the Cluster Table

**Scripts**

```text
make_US_cluster_table.py
run_US_cluster_table.pbs
```

Creates the complete mapping:

```text
US.ORF → US.100NT → US.95NT
```

**Output**

```text
US.clusters.tsv.xz
```

**Rows**

```text
265,392,705
```

---

## 9. eggNOG Annotation

**Script**

```text
Run_Eggnog_Mapper.sh
```

Runs eggNOG-mapper to generate functional annotations for the Urban Soil ORFs.

---

## 10. Transfer Annotations to 95NT

**Scripts**

```text
US95_emapper.py
run_US95_emapper.pbs
```

Transfers ORF-level eggNOG annotations to the corresponding 95NT representatives.

**Output**

```text
US.95NT.emapper.annotations.gz
```

**Result**

```text
95NT representatives  199,172,868
Annotated               13,656,625
Unannotated             185,516,243
Samples                          58
```

## Summary

| Stage | Count |
|---|---:|
| Original ORFs | 265,392,705 |
| 100NT representatives | 256,829,048 |
| 95NT representatives | 199,172,868 |
| Annotated 95NT representatives | 13,656,625 |
