# RNA-Benchmark

This repository contains resources and information for the RNA Benchmark Project, which evaluates the performance of various computational models in predicting RNA structures, including single RNA, RNA/RNA and RNA/Protein structures. Below, you'll find details about the dataset, selection criteria, and included models.

## Dataset

The benchmark dataset was derived from the Protein Data Bank (PDB) with a cutoff date of September 30, 2021. This date aligns with the latest data used for training the AlphaFold 3 and HelixFold 3 models. The dataset includes RNA structures that meet specific criteria to ensure consistency and quality:

**Sequence Identity Filtering:**  

Structures were filtered using CD-HIT-EST to exclude those with more than 80% sequence identity. Sequence identity was calculated by aligning sequences to the longer sequence, with up to 80% overlap considered.  

**Nucleotide Content:**  

Structures containing ambiguous nucleotides, such as 'X' or 'N', were excluded. Sequences that are shorter than 10 nt and contain same nucleotide type (80-100% of the same nucleotide in the sequence) are most likely linear and have no modular structure, are short and are therefore excluded from the benchmark datset. This was mostly the case with some RNA/Protein synthetic complexes.

**MSA Generation:**

Multiple sequence alignments (MSAs) for the analysis were generated using **rMSA** for RNA sequences and **MMseqs2** for protein sequences. The MSAs are provided in the default output formats of these tools: *.afa* for RNA and *.a3m* for proteins. Depending on the specific model requirements, MSAs were converted into alternative formats as needed. It is important to note that rMSA converts the 'U' nucleotide to 'T' in RNA sequence alignments.

MSAs can be found at:
'''
cd MSA/RNA
cd MSA/PROTEIN
'''



## Models Included in the Benchmark

The benchmark evaluates the performance of the following **nine** computational models: 

1) AlphaFold 3
2) HelixFold 3
3) RosettaFoldNA
4) DRFold
5) trRosetta
6) RhoFold+
7) NuFold
8) Chai
9) Boltz

HelixFold 3, Boltz, and Chai are derived from AlphaFold 3 and considered AlphaFold-Based Models, while RhoFold+ and similar models leverage large language models (LLMs) for structure prediction (LLM-Based Models).

