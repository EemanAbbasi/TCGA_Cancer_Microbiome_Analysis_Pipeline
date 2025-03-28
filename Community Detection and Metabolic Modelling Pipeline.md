Mutualism and Competition Analysis Pipeline
This repository contains a GitHub Actions pipeline to analyze mutualistic and competitive interactions between microbial species based on OTU data, as outlined in the sections Determining co-occurring communities and Determining metabolic interactions. The pipeline filters OTU tables, constructs co-occurrence networks, detects communities, reconstructs metabolic models using CarvMe, and evaluates metabolic interactions using SMETANA.

Overview
The pipeline performs the following steps:

Filter OTU Table: Retains species present in >50% of samples.
Correlation Analysis: Computes Spearman correlations (|corr| > 0.6, adjusted p < 0.01) to build a co-occurrence network.
Community Detection: Uses the fast greedy algorithm (Clauset et al., 2004) to identify sub-communities (>10 species).
Genome Download: Fetches RefSeq genomes from NCBI for species in detected communities.
Metabolic Model Reconstruction: Builds genome-scale metabolic models with CarvMe using default, M9, and LB media settings.
Metabolic Interactions: Calculates competition and mutualism scores with SMETANA (Zelezniak et al., 2015).
In-Silico Simulations: Simulates interactions of differentially enriched species with varying group sizes (4, 6, 8 species, 10 repeats).
