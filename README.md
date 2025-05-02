# RNA-Benchmark

This repository contains resources and information for the RNA Benchmark Project, which evaluates the performance of various computational models in predicting RNA structures, including single RNA, RNA/RNA and RNA/Protein structures. Below, you'll find details about the dataset, selection criteria, and included models.

## Dataset

The benchmark dataset was derived from the Protein Data Bank (PDB) with a cutoff date of September 30, 2021. This date aligns with the latest data used for training the AlphaFold 3 and most models. The dataset includes RNA structures that meet specific criteria to ensure consistency and quality:

**Sequence Identity Filtering:**  

Structures were filtered using CD-HIT-EST to exclude those with more than 80% sequence identity. Sequence identity was calculated by aligning sequences to the longer sequence, with up to 80% overlap considered.  

**Nucleotide Content:**  

Structures containing ambiguous nucleotides, such as 'X' or 'N', were excluded. Sequences that are shorter than 10 nt and contain same nucleotide type (80-100% of the same nucleotide in the sequence) are most likely linear and have no modular structure, are short and are therefore excluded from the benchmark datset. This was mostly the case with some RNA/Protein synthetic complexes.

**MSA Generation:**

Multiple sequence alignments (MSAs) for the analysis were generated using **rMSA** for RNA sequences and **MMseqs2** for protein sequences. The MSAs are provided in the default output formats of these tools: *.afa* for RNA and *.a3m* for proteins. Depending on the specific model requirements, MSAs were converted into alternative formats as needed. It is important to note that rMSA converts the 'U' nucleotide to 'T' in RNA sequence alignments.

## Models Included in the Benchmark

The benchmark evaluates the performance of the following **eight** computational models: 

1) AlphaFold 3
2) HelixFold 3
3) RosettaFoldNA
4) trRosetta
5) RhoFold+
6) NuFold
7) Chai
8) Boltz

HelixFold 3, Boltz, Chai and DRFold are derived from AlphaFold 3 and considered AlphaFold-Based Models, while RhoFold+ and similar models leverage large language models (LLMs) for structure prediction (LLM-Based Models). Many models allow only for single RNA inference (NuFold, RhoFold+, trRosetta), while others are able to predict both monomeric and oligomeric RNAs, alone and in a complex with proteins.

❗Chai was excluded from the RNA complex analysis, as its inference is currently limited to a maximum of 2048 tokens, and excluding RNA complexes with >2048 tokens would further reduce the size of our comparable PDB data. Boltz failed to predict copmlexes consisting of 5+ protein+RNA subunits due to "out-of-memory" issues, even after the update to a newer 0.4.0 version. Therefore, such complexes were not used in our analysis of benchmark data. This, along with filtering step resulted in a dataset comprising 50 single RNA (momnomeric RNA) and 46 RNA complexes. File **meta.csv** contains only those IDs that were used in the analysis, and should be used as reference.

## Data

```
RNA-Benchmark/
├── FASTA/                        # fasta sequences of single RNA, RNA/RNA and RNA/Protein complexes 
│   ├── single        
│   ├── complex       
├── PDB/
|   ├── pdb.zip                   # top-ranked predictions    
│   ├── true_pdb_complexes.zip    # solved structures in pdb format downloaded from PDB
│   ├── true_pdb_single.zip 
├── MSA/
│   ├── msas.zip                  # MSAs used to generate structures for both single RNA and complexes
├── OUTPUTS/                      # all CSV files of benchmarked data and scores
├── scr/                          # code to generate all figures
├── meta.csv                      # information about sequences/structures used in the final analysis
```

❗ To generate all the plots from the benchmark analysis, run 