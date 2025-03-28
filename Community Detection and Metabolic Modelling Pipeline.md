Filter OTU Table: Keep species present in >50% of samples from otu_table.csv.
Correlation Analysis: Compute Spearman correlation matrix; filter pairs with |corr| > 0.6 and adjusted p-value < 0.01 to build a co-occurrence network.
Community Detection: Apply fast greedy algorithm (Clauset et al., 2004) to identify sub-communities (>10 species).
Download Genomes: Fetch reference genomes for species in sub-communities from NCBI.
Metabolic Model Reconstruction: Use CarvMe to build genome-scale metabolic models from genomes, testing default, M9, and LB media options.
Metabolic Interactions: Run SMETANA (Zelezniak et al., 2015) to calculate competition and mutualism scores for sub-communities.
In-Silico Simulations: Simulate interactions of differentially enriched species with 4, 6, or 8 randomly selected species (10 repeats each) and compare to baseline (random species excluding the target).
