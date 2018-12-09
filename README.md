
<!-- README.md is generated from README.Rmd. Please edit that file -->
proDD
=====

Differential Detection for Label-free (LFQ) Mass Spectometry Data

The tool fits a probabilistic ropout model to an intensity matrix from from label-free quantification (LFQ). Dropouts in LFQ data occur if the protein has a low intensity. Our model takes the non-random missingness into account, by constructing a Bayesian hierarchical model. After fitting the model you can sample from the posterior distribution of the means from each protein and condition. The posterior are a useful element to calculate all kind of statistics/metrics including the probability that the intensity of a protein in one condition is smaller than in the control (similar to the one-sided p-value).

Installation
============

Install the latest version directly from GitHub (make sure that `devtools` is installed)

``` r
devtools::install("const-ae/proDD")
```

Disclaimer
==========

I am still actively working on the project and although the algorithm is working fine at this point, the API
might still be subject to change.
