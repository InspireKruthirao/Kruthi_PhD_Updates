# 🧬 Small Proteins Functional Analysis via Gene Synteny

⭐ Star us on GitHub — it motivates us a lot!

#add image here ( Tryhttps://github.com/matiassingers/awesome-readme?tab=readme-ov-file)

This repository contains code, data, and documentation for my **PhD research** focused on identifying **functions of small proteins** using **gene synteny** and related computational analyses.

---

## 🔬 Project Overview

Small proteins play critical roles in biological systems, but their functions are often poorly understood. This project explores **gene neighborhood patterns**, to predict functional associations.

## 📥 Input Data (smORF Catalog from Zenodo)

This project uses the **Small ORFs (smORFs) catalog** from Zenodo:

🔗 **Source:**  
https://zenodo.org/records/16356977  

### 📊 Dataset Summary

- **403,491** non-redundant smORFs  
- **273,065** clusters at **90% amino acid identity**  
- smORFs derived from multiple samples (identical sequences may occur in multiple origins)

### 📂 Files Used in This Project

| File | Description |
|------|-------------|
| `SHD1_SM.100AA.faa.gz` | FASTA of smORFs clustered at **100% amino acid identity** |
| `SHD1_SM.100AA_origins.tsv.gz` | Metadata mapping smORFs to original samples and genomic origins |
| `SHD1_SM.90AA.faa.gz` | FASTA of smORFs clustered at **90% amino acid identity** |
| `SHD1_SM.clusters.tsv.gz` | Cluster table for **90% identity clusters** |

###  Notes

- The same sequence identifier may appear multiple times in the metadata file because identical smORFs may be present in multiple samples.
- Use **100% AA clusters** for exact sequence-level analyses.
- Use **90% AA clusters** for functional grouping and conserved synteny analyses.

WIP
- Map neighboring genes to identify conserved synteny.
- Cluster ORFs to reveal functionally related protein groups.
- Annotate clusters using eggNOG for functional insights.
- Visualize gene clusters and synteny patterns.

---

## 📁 Repository Structure


## Authors and acknowledgment
