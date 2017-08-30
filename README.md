# maRge (Motif Analysis in R for Genomic Exploration)

<img src="https://media.giphy.com/media/l2JehBcFwVdTlFRu0/giphy.gif" width="250">

`[WARNING: PACKAGE UNDER CONSTRUCTION]`


## Description
maRge is a package developed for simplifying motif analysis of genomic regions. 
On one side, it provides a wrap to HOMER de novo motif analysis. On the other side, 
it provides easy tools to query a region for transcription factor binding motifs.

## Installation
```
## install.packages("devtools")
devtools::install_github("mirthelle/maRge")
library(maRge)
```
## Functionality
### Run HOMER from R
We have implemented several functions for running HOMER programs from (they are actually just calls to the system with the command line specified in the arguments of the function).

- `deNovoMotifHOMER`. Runs `findMotifsGenome.pl` to obtain lists of motifs enriched in your set of regions.
- `countSignHOMER`. This function is usefull when you want to use **de novo** motifs for your analysis. It returns the number of significant **de novo** motifs that were found in your data, so you can extract the rows from the dataframe in downstream analysis.
- `catSignMotifsHOMER`. Creates a single `.motif` file containing all the sequences that were found significant in the analysis.

## Query a single region for TF Motifs


## Generation of Transcription Factor Regulatory Networks
