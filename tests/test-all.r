#!/usr/bin/Rscript

library(testthat, quietly=T)
library(quantumClustering, quietly=T)

# Comment out if you do not have doMC installed and/or do not want to
# run processes in parallel
library(doMC, quietly=T)
registerDoMC(cores=4)

# Get rid of the species in order to cluster only on numeric data
iris_clustering <- iris[, setdiff(colnames(iris), "Species")]

test_that("Achieves proper iris clusters for sigma_factor=1", {
    sigma <- calcSD(iris_clustering)
    clusters <- qc(iris_clustering, sigma=sigma, min_d_factor=2)
    expected_out <- c(rep(2, 50), rep(1, 100))
    expect_equal(clusters, expected_out)
})

test_that("Achieves proper iris clusters for sigma_factor=1.85", {
    sigma <- calcSD(iris_clustering) / 4
    clusters <- qc(iris_clustering, sigma=sigma, min_d_factor=2)
    expected_out <-
        c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
          2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,
          3,3,3,4,1,4,3,3,3,3,3,1,1,3,3,5,4,1,3,1,4,1,3,3,1,1,3,3,4,5,3,1,1,4,
          3,3,1,3,3,3,1,3,3,3,1,3,3,1)
    expect_equal(clusters, expected_out)
})

test_that("Cluster limitation removes first cluster", {
    sigma <- calcSD(iris_clustering) / 2
    # Should be the same as sigma_factor=1, or clusters of 50 and 100
    # Since the second is bigger, only it will be accepted.  Therefore,
    # cluster 1 becomes 0 (unclustered) and cluster 2 becomes 1.
    clusters <- qc(iris_clustering, sigma=sigma, min_d_factor=2,
                   n_clusters_max=1)

    expected_out <- c(rep(0, 50), rep(1, 100))
    expect_equal(clusters, expected_out)
})

print("Success!")
