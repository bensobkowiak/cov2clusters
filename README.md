# cov2clusters - Stable clustering of SARS-CoV-2 sequences from phylogenetic trees

This repository contains two R scripts (cov2clusters.R and clusteringFunctions.R) for genomic clustering and assigning cluster names of SARS-Cov-2 sequences. There is also a folder containing demo data and results.

The input is an untimed phylogenetic tree in Newick format and collection dates, either in a .csv file or .json file from a timed phylogeny. A further variable can be added into the metadata file to restrict clustering to samples sharing this feature.

To run with with default parameters and dates in a .csv file with two columns (ID, date):

> cov2clusters(treeName = "tree.nwk", metafile = "dates.csv", json_dates = FALSE)

Or with a .json file:

> cov2clusters(treeName = "tree.nwk", json_file = dates.json)

Please cite the following preprint if you choose to use this approach in published work:

https://www.medrxiv.org/content/10.1101/2022.03.10.22272213v2

This code has been tested on R version 3.6.0 (2019-04-26).
