# cov2clusters - Stable clustering of SARS-CoV-2 sequences from phylogenetic trees

This repository contains scripts written in R (cov2clusters.R) for genomic clustering and assigning cluster names of SARS-Cov-2 sequences. 

The input is an untimed phylogenetic tree in Newick format and collection dates, either in a .csv file or .json file from a timed phylogeny. A further variable can be added into the metadata file to restrict clustering to samples sharing this feature.

To run with with default parameters and dates in a .csv file with two columns (sample name, date):

cov2clusters(treeName = "tree.nwk", metafile = "dates.csv", json_dates = FALSE)

Or with a .json file:

cov2clusters(treeName = "tree.nwk", json_file = dates.json)

