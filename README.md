# GhostPopulationSimulation
Population Divergence with 3 populations, with 1 population acting as a ghost but still able to engage in gene flow

msprime3pop.py is an msprime tskit simulator that simulates lineages backwards in time. it also contains a tskit to vcf converter that creates vcf files that are able to be analyzed with the PPP.py

Models 1-15 are examples of varying levels of gene flow between populations
The figures for models 1-15 is labeled in "gene flow figures" bilateral and unilateral
Deep_Divergence_Simulated_Data is newly created data to help supplement "Accounting for Gene Flow from Unsampled ‘Ghost’ Populations while Estimating Evolutionay History under the Isolation with Migration Model"
Arun Sethuraman, Melissa Lynch
doi: https://doi.org/10.1101/733600
In this new data, models A-E are instead 1-5 with an added Model6 as control with an ancient divergence time of 2,000,000 years and a more "recent" second divergence at 500,000 years 
Recent_Divergence_Simulated_Data follows the same pattern above, but with an ancient divergence time of 18500 years and a recent divergence of 1850 years. 