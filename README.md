# RNA-Benchmark

This repository contains resources and information for the RNA Benchmark Project, which evaluates the performance of various computational models in predicting RNA structures, including single RNA, RNA/RNA and RNA/Protein structures. Below, you'll find details about the dataset, selection criteria, and included models.

## Dataset

The benchmark dataset was derived from the Protein Data Bank (PDB) with a cutoff date of September 30, 2021. This date aligns with the latest data used for training the AlphaFold 3 and HelixFold 3 models. The dataset includes RNA structures that meet specific criteria to ensure consistency and quality:

**Sequence Identity Filtering:**  

Structures were filtered using CD-HIT-EST to exclude those with more than 80% sequence identity. Sequence identity was calculated by aligning sequences to the longer sequence, with up to 80% overlap considered.  

**Nucleotide Content:**  

Structures containing ambiguous nucleotides, such as 'X' or 'N', were excluded.

## Models Included in the Benchmark

The benchmark evaluates the performance of the following **nine** computational models: AlphaFold 3, HelixFold 3, RosettaFoldNA, DRFold, trRosetta, RhoFold+, NuFold, Chai and Boltz. HelixFold 3, Boltz, and Chai are derived from AlphaFold 3 and considered AlphaFold-Based Models, while RhoFold+ and similar models leverage large language models (LLMs) for structure prediction (LLM-Based Models).

