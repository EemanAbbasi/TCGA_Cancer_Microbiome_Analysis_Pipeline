# Mutualism and Competition Analysis Pipeline

This repository contains a GitHub Actions pipeline to analyze mutualistic and competitive interactions between microbial species based on OTU (Operational Taxonomic Unit) data. The methodology is inspired by the sections *Determining co-occurring communities* and *Determining metabolic interactions* from relevant literature. The pipeline filters OTU tables, constructs co-occurrence networks, detects communities, reconstructs metabolic models using `CarveMe` [1], and evaluates metabolic interactions using `SMETANA` [2].

## Overview

The pipeline performs the following steps:

- **Filter OTU Table**: Retains species present in >50% of samples to ensure robust representation.
- **Correlation Analysis**: Computes Spearman correlations (|corr| > 0.6, adjusted p < 0.01) to construct a co-occurrence network.
- **Community Detection**: Applies the fast greedy algorithm (Clauset et al., 2004) to identify sub-communities with >10 species.
- **Genome Download**: Retrieves RefSeq genomes from NCBI for species in detected communities.
- **Metabolic Model Reconstruction**: Constructs genome-scale metabolic models using `CarveMe` with default, M9, and LB media settings.
- **Metabolic Interactions**: Calculates competition and mutualism scores using `SMETANA`.
- **In-Silico Simulations**: Simulates interactions of differentially enriched species with group sizes of 4, 6, and 8 species (10 repeats each).

## Prerequisites

- **Python**: Version 3.9 or higher.
- **Dependencies**:
  - `numpy`, `scipy`, `pandas`, `networkx` (for data processing and community detection).
  - `carveme` (for metabolic model reconstruction).
  - `smetana` (for metabolic interaction analysis).
  - `requests` (for downloading genomes from NCBI).
- **Input Data**: An OTU table (e.g., `input_data/otu_table.csv`) with species abundance across samples.

## Setup and Usage

1. **Carve a metaboluc model**  
   Instead of manually providing genome data, use `CarveMe` to fetch sequences from NCBI RefSeq by specifying an accession code. For example:
   ```bash
   carve --refseq GCF_000005845.2 -o ecoli_k12_mg1655.xml
   ```
    This command downloads the *Escherichia coli* K-12 MG1655 genome and builds a metabolic model.

3. **Gapfill with the relevant growth media**  
Refine the reconstructed model by filling metabolic gaps using a specified medium (e.g., M9):
```bash
carvme gapfill ecoli_k12_mg1655.xml -m M9 -o ecoli_k12_mg1655_gapfilled.xml
```

3. **Analyze Metabolic Interactions**  
Use `SMETANA` to compute interaction scores for communities. Provide the gapfilled models and a TSV file listing community members (e.g., `communities.tsv`):
```bash
smetana_main.py *.xml 
```
This leverages code from `smetana_main.py` in the metabolic interaction folder.

## References

1. Machado, D., Andrejev, S., Tramontano, M., & Patil, K. R. (2018). Fast automated reconstruction of genome-scale metabolic models for microbial species and communities. *Nucleic Acids Research*, 46(15), 7542–7553. https://doi.org/10.1093/nar/gky537 (`CarveMe`)  
2. Zelezniak, A., Andrejev, S., Ponomarova, O., Mende, D. R., Bork, P., & Patil, K. R. (2015). Metabolic dependencies drive species co-occurrence in diverse microbial communities. *Proceedings of the National Academy of Sciences*, 112(20), 6449–6454. https://doi.org/10.1073/pnas.1421834112 (`SMETANA`)  
