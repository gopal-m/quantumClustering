# Quantum clustering

This is an implementation of David Horn's quantum clustering algorithm in R.  It does not suffer from the curse of dimensionality and takes advantage of eignenfunctions for non-linear clustering.  It was adapted from http://horn.tau.ac.il/software/qc.m with help from the following paper: Algorithm for Data Clustering Pattern Recognition Problems Based on Quantum Mechanics (Physical Review Letters 88, 2002).

## Dependencies

* foreach

## Suggested

* doMC -- strongly recommended, for running the algorithm on multiple cores with greater speed
* devtools -- optional, for building the package
* testthat -- optional, for running tests to verify functionality

## Usage:

```
library(quantumClustering)
clusters <- qc(dataset, sigma, steps=21, min_d_factor=2,
               n_clusters_max=1000, verbose=FALSE)
```

The qc() function returns a vector of the designated clusters for each row in dataset.  See help(qc) for more information or tests/ for more examples.


## Installation:

```
# If you do not already have the "devtools" package installed
install.packages("devtools")

devtools::install_github("rdtaylor/quantumClustering")
```
