# Urban Soil smORF/ORF Pipeline

## Overview

This directory contains the scripts used to build a non-redundant ORF and smORF catalogue for the Urban Soil metagenome dataset.

The workflow is based on the approach previously used for the Shanghai Dogs (SHD) dataset, with modifications to handle the much larger Urban Soil dataset of approximately **265 million predicted ORFs**.

## Pipeline

### 1. ORF Prediction

**Script:** `run_prodigal_array.sh`

Runs Prodigal on the Urban Soil assemblies to predict open reading frames (ORFs).

### 2. ORF Catalogue Generation

**Scripts:**
`build_urban_soil_orf_catalog.sh`
`build_urban_soil_orf_catalog.py`

Combines the predicted ORFs into a single Urban Soil ORF catalogue and compresses the output.

**Output:**

`UrbanSoil.ORF.fna.xz`

### 3. 100% Nucleotide Deduplication

**Scripts:**
`full_dedup_100NT.sh`
`postprocess_cdhit_clusters.py`

Removes identical ORF sequences using `cd-hit-est` at 100% nucleotide identity and processes the resulting cluster information.

### 4. smORF Identification

**Script:** `Run_GMSC_MAPPER.sh`

Maps the Urban Soil ORFs against the GMSC resource to identify candidate small ORFs (smORFs).

### 5. smORF Resource Generation

**Script:** `Create_resource_files.py`

Creates the main resource files required for downstream smORF analysis, including:

* Protein sequences (`.faa.gz`)
* Sequence origins, including sample, contig and genomic coordinates
* Habitat information
* Taxonomic information

### 6. smORF Clustering

**Scripts:**
`Cluster_SmORFs.sh`
`script.py`

Clusters smORFs at 90% sequence identity using CD-HIT.

**Main outputs:**

* `UrbanSoil_Clusters.tsv.gz`
* `UrbanSoil_90AA_SMORFs.faa.gz`

### 7. Occurrence Statistics

**Script:** `Run_Stats_data_US.py`

Calculates smORF occurrence statistics and examines their distribution across the Urban Soil dataset.

### 8. Functional Annotation

**Script:** `run_eggnog.sh`

Runs eggNOG-mapper to obtain functional annotations for the relevant sequences.

### 9. Gene-Neighborhood Analysis

**Script:** `Gene_Neighborhood.py`

Uses genomic neighborhood information to investigate and predict potential functions of small proteins based on their surrounding genes and conserved gene context.

## Workflow Summary

The Urban Soil analysis consists of four main stages:

1. **ORF prediction and catalogue generation**
2. **Sequence deduplication and smORF identification**
3. **smORF clustering and functional annotation**
4. **Gene-neighborhood-based functional prediction**

The final goal is to use conserved genomic context to support functional prediction of poorly annotated or unknown small proteins in the Urban Soil metagenome dataset.
