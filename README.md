# Small Protein Functional Analysis Using Gene Neighbourhoods

This repository contains scripts, workflows, figures, and results investigating the functional prediction of small proteins (smORFs) using genomic neighbourhood information.


## Workflow at a glance

![Illustrated workflow of the smORF gene-neighbourhood prediction pipeline](assets/smorf_neighbourhood_workflow.png)

*From metagenomic material to short-protein discovery, genomic-neighbourhood context, functional annotation and evidence-weighted prediction. Purple marks the target smORF; teal and blue mark neighbouring annotated genes.*

## Research objective

The project examines whether conserved neighbouring genes can help predict the functions of poorly characterised small proteins (≤100 amino acids).

## Datasets

### Shanghai Dogs

The Shanghai Dogs long-read metagenomic dataset is used as the proof-of-concept dataset for developing and validating the gene-neighbourhood prediction workflow.

The `Shanghai_Dogs/` directory contains:

- Neighbourhood extraction scripts
- Prediction and validation scripts
- Statistical analyses
- Figures
- Results

### Urban Soil

The Urban Soil dataset is used to scale the workflow to a substantially larger metagenomic dataset.

The `Urban_Soil/` directory contains:

- Gene catalogue construction
- smORF catalogue construction
- Sequence clustering
- Functional annotation
- Gene-neighbourhood analysis

## General workflow

1. Predict protein-coding genes and smORFs.
2. Construct non-redundant gene and smORF catalogues.
3. Annotate neighbouring genes using eggNOG-mapper.
4. Extract genomic neighbourhoods surrounding each smORF.
5. Predict smORF functions using neighbourhood-based voting.
6. Evaluate prediction confidence and accuracy.
7. Generate statistical summaries and figures.
