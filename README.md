# lncRna <img src="inst/img/lncRna_logo.png" align="right" height = 150/>

![GitHub](https://img.shields.io/github/license/prodakt/lncRna?style=plastic)
![GitHub top language](https://img.shields.io/github/languages/top/prodakt/lncRna?style=plastic)
![GitHub R package version](https://img.shields.io/github/r-package/v/prodakt/lncRna?color=gree&label=ver&style=plastic)
![GitHub repo file count](https://img.shields.io/github/directory-file-count/prodakt/lncRna?color=purple&style=plastic)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/prodakt/lncRna?color=white&label=activ)
![GitHub last commit](https://img.shields.io/github/last-commit/prodakt/lncRna?color=%23FF0000%20&label=last)

## Introduction

**lncRna** is an R package designed to simplify and streamline the process of long non-coding RNA (lncRNA) identification and functional analysis.

Accurate identification and functional annotation of lncRNAs are critical for understanding their regulatory roles in complex biological processes. This package provides a comprehensive, user-friendly interface for processing transcriptomic data, implementing lncRNA identification pipelines, and systematically evaluating the performance of multiple coding-potential prediction tools.

## Core functionality

* **Easy Pipeline Setup:** Filter transcripts from GTF files based on structural features and expression levels.

* **Coding Potential Assessment:** Integrate and evaluate results from multiple external prediction tools (e.g., CPC2, PLEK, CPAT).

* **Statistical Evaluation:** Compute key metrics (accuracy, precision, sensitivity) to select the most reliable lncRNA identification strategy (including tool combinations and "at least N" consensus).

* **Functional Annotation:** Identify potential *cis*- and *trans*-acting interactions with protein-coding genes and perform functional enrichment analysis.

* **Visualization:** Generate informative plots (radar plots, clock plots, interactive Sankey diagrams, Venn diagrams) to visualize performance and interaction results.

## Example Visualizations

The package allows for easy generation of various plots, such as Venn diagrams for coding potential results:

![venn_all](https://github.com/prodakt/lncRna/blob/main/inst/img/vennFull.png)
![venn_4](https://github.com/prodakt/lncRna/blob/main/inst/img/vennSel.png)

## Installation

You can install the `lncRna` package using `BiocManager` from Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("lncRna")
```
or directly from Github
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("prodakt/lncRna")
```

## Getting Started

To get started with `lncRna`, please read the package vignette which contains a complete, step-by-step tutorial on how to use the package functions:

```r
library(lncRna)
vignette("lncRna")
