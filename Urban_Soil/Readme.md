Urban Soil smORF/ORF Pipeline

This directory contains the scripts used to build a non-redundant ORF and smORF catalog for the Urban Soil metagenome dataset, mirroring the approach used for the 
Shanghai Dogs (SHD) dataset, with adaptations for Urban Soil's much larger data volume 265M ORFs.

Assemblies
    │
    ▼
run_prodigal_array.sh          Predict ORFs from assemblies (Prodigal)
    │
    ▼
build_urban_soil_orf_catalog.sh/.py
    Build/compress the raw ORF catalog (UrbanSoil.ORF.fna.xz)
    │
    ▼
full_dedup_100NT.sh            Full-scale 100% identity dedup (cd-hit-est)
  ├─ postprocess_cdhit_clusters.py
    │
    ▼
Run_GMSC_MAPPER.sh             Map ORFs against GMSC to identify smORFs
    │
    ▼
Create_resource_files.py       Build 100AA smORF resource files:
                                  - sequences (.faa.gz)
                                  - origins (sample/contig/coords)
                                  - habitat & taxonomy
    │
    ▼
Cluster_SmORFs.sh / script.py  90% identity clustering of smORFs (cd-hit)
                                  -> UrbanSoil_Clusters.tsv.gz
                                  -> UrbanSoil_90AA_SMORFs.faa.gz
    │
    ▼
Run_Stats_data_US.py           Occurrence statistics / distribution
    │
    ▼
run_eggnog.sh                  Functional annotation (eggNOG-mapper)
    │
    ▼
[Gene_Neighborhood.py]         Neighbourhood-based function prediction
