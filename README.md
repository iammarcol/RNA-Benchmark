# RNA-Benchmark

This repository contains resources and information for the RNA Benchmark Project, which evaluates the performance of various computational models in predicting RNA structures, including single RNA, RNA/RNA and RNA/Protein structures. Below, you'll find details about the dataset, selection criteria, and included models.

## Dataset

The benchmark dataset was derived from the Protein Data Bank (PDB) with a cutoff date of September 30, 2021. This date aligns with the latest data used for training the AlphaFold 3 and HelixFold 3 models. The dataset includes RNA structures that meet specific criteria to ensure consistency and quality:

**Sequence Identity Filtering:**  

Structures were filtered using CD-HIT-EST to exclude those with more than 80% sequence identity. Sequence identity was calculated by aligning sequences to the longer sequence, with up to 80% overlap considered.  

**Nucleotide Content:**  

Structures containing ambiguous nucleotides, such as 'X' or 'N', were excluded. Sequences that are shorter than 10 nt and contain same nucleotide type (80-100% of the same nucleotide in the sequence) are most likely linear and have no modular structure, are short and are therefore excluded from the benchmark datset. This was mostly the case with some RNA/Protein synthetic complexes.

❗This filtering resulted in a dataset comprising 46 single RNA complexes, 16 RNA/RNA complexes, and 40 RNA/Protein complexes. (See benchmark_metadata.csv for more information)

**MSA Generation:**

Multiple sequence alignments (MSAs) for the analysis were generated using **rMSA** for RNA sequences and **MMseqs2** for protein sequences. The MSAs are provided in the default output formats of these tools: *.afa* for RNA and *.a3m* for proteins. Depending on the specific model requirements, MSAs were converted into alternative formats as needed. It is important to note that rMSA converts the 'U' nucleotide to 'T' in RNA sequence alignments.

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

HelixFold 3, Boltz, Chai and DRFold are derived from AlphaFold 3 and considered AlphaFold-Based Models, while RhoFold+ and similar models leverage large language models (LLMs) for structure prediction (LLM-Based Models). Many models allow only for single RNA inference (DRFold, NuFold, RhoFold+, trRosetta), while other are able to predict both monomeric and oligomeric RNAs, alone and in a complex with proteins.

❗Chai was excluded from the RNA complex analysis, as its inference is currently limited to a maximum of 2048 tokens, and excluding RNA complexes with >2048 tokens would further reduce the size of our comparable PDB data. Boltz failed to predict copmlexes consisting of 5+ protein+RNA subunits due to "out of memory" issue, even after the update to a newer 0.4.0 version. Therefore, such complexes were not used in our analysis of benchmark data. That reduced our benchmark dataset for RNA/RNA and RNA/Protein complexes from 56 to 43 examples.

## Data

```
RNA-Benchmark/
├── FASTA/                        # fasta sequences of single RNA, RNA/RNA and RNA/Protein complexes 
│   ├── single        
│   ├── complex       
├── PDB/                          # solved structures in pdb format downloaded from PDB
│   ├── true_pdb_complexes.zip 
│   ├── true_pdb_single.zip 
├── PRED/
│   ├── pred_pdb_complexes.zip
│   ├── pred_pdb_single.zip              # predicted structures for both single RNA and complexes
├── MSA/
│   ├── msas.zip                  # MSAs used to generate structures for both single RNA and complexes
├── meta/
│   ├── benchmarking_metadata.zip # information about sequences used in the analysis
```



