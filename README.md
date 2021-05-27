# lncRna
is the package for lncRNA identification i a simple way.

This package allows to carry out the lncRNA identification pipeline in a few very simple steps. In general the pipieline consist of three stages:
- extraction of known lncRNA from reference GTF
- filtering potential (unknown) lncRNA based on the features such as sequence length, exon number
- lncRNA classification based on the coding potential and the simmilarity with known nucleotide and peptide sequences (homologues in databases).


## How to install
First you need to install and load the '[devtools!](https://github.com/r-lib/devtools)' package in R. 
```
# Install devtools from CRAN
install.packages("devtools")
# load the library
library(devtools)
```
Next you should download the 'lncRna' library from GitHub and load it in your environment:
```
install_github("prodakt/lncRna")
library("lncRna")
```

## rapid pipeline (first fast run)

1. First you need to install/load the '[rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)' library to import any GTF and GFF files in a acurate format.
```
install.packages("BiocManager")
library("rtracklayer")
```

2. Run `??lncRna` to list the available functions in the library.
3. xxx
