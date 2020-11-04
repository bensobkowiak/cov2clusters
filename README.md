# cov2clustering (Clustering of SARS-CoV-2 using whole genome sequences)

This repository contains scripts written in R (cov2clustering.R) for clustering and assigning cluster names of SARS-Cov-2 sequences. 

The input is a SNP or genetic distance matrix, which can be computed rapidly using the ape package in R, PAUP, or with https://github.com/tseemann/snp-dists. 

Temporal and spatial metadata can be included to refine clusters in the form of a .csv file.

