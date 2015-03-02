# Adapted from http://horn.tau.ac.il/software/qc.m
#
# David Horn and Assaf Gottlieb.  Algorithm for Data Clustering Pattern
# Recognition Problems Based on Quantum Mechanics (Physical Review Letters
# 88, 2002).

#' Summarize standard deviation
#'
#' Summarize the standard deviation of a data set with multiple columns into a single value.  This function accomplishes this task by finding the length of a vector containing the standard deviation of each column.  This function may be useful in establishing a \code{sigma} for qc().
#'
#' @param dataset The cleaned input data to be clustered in either data.frame or matrix format.  This should contain only numeric data and must not have any factors, strings, NAs, etc.  It also should not contain any irrelevant columns such as observation ID or redundant data.
#'
#' @return The function returns a double containing the approximated standard deviation of the entire data set.
#'
#' @examples 
#' data(iris)
#' iris_data <- iris[, -5] # Drop the species factor
#' sigma <- 0.4 * calcSD(iris_data)
#'
#' @export
calcSD <- function(dataset) {
    sd <- c()
    p <- ncol(dataset)

    if(is.null(p)) {
        sd <- sd(dataset)
    } else {
        for(i in seq(1, p)) {
            sd <- c(sd, sd(dataset[, i]))
        }
        sd <- (sum(sd^2)) ^ (1/2)
    }

    sd
}

#' Cluster data
#'
#' Cluster data using a non-linear clustering algorithm which finds organic groups of data from eigenfunctions and does not suffer from the curse of dimensionality.  The algorithm is approximately O(n^2).
#'
#' @param dataset The cleaned input data to be clustered in either data.frame or matrix format.  This should contain only numeric data and must not have any factors, strings, NAs, etc.  It also should not contain any irrelevant columns such as observation ID or redundant data.
#' @param sigma A double which controls how closely related data in clusters should be.  The smaller the number, more clusters will be created with fewer observations in each.  If sigma is too small, observations either will not be clustered or will be in their own individual clusters.  If sigma is too large, most -- if not all -- observations will be in the first cluster.
#' @param steps An integer specifying the number of expectation-maximization steps to take.  If faster, less accurate results are required, this may be reduced from the default of 21.
#' @param min_d_factor A double which controls how close data points must be in order to be considered in the same cluster.  Specifically, this value is the number of sigmas of distance to be within said threshold.  This value should probably not be changed unless there is a strong reason to do so.
#' @param n_clusters_max An integer specifying the maximum number of clusters to return.  These clusters will always be the most common clusters with the most observations in them.
#' @param verbose A boolean value which toggles the algorithm's verbosity for details as to how far along it is.
#'
#' @return The function returns a vector of the numeric clusters assigned to each row in \code{dataset}.
#'
#' @examples 
#' # Set up parallel execution
#' library(doMC)
#' registerDoMC(cores=2) # Replace `2` with the number of cores in your machine
#'
#' data(iris)
#' iris_data <- iris[, -5] # Remove the classification of the iris data set
#' clusters <- qc(iris_data, 1.0)
#' secondary_cluster_rows <- which(clusters == 2)
#' print(iris[secondary_cluster_rows, ])
#'
#' @export
qc <- function(dataset, sigma, steps=21, min_d_factor=2,
               n_clusters_max=1000, verbose=FALSE) {

    ETA_DECAY <- 0.9

    main <- function(dataset, sigma, steps, min_d_factor,
                     n_clusters_max) {
        min_d <- getMinD(sigma, min_d_factor)
        eta <- getEta(sigma)
        q <- getQ(sigma)

        qc_log(paste0("Starting gradient decent with sigma ", sigma,
                      " and q ", q))
        D <- gradDesc(dataset, dataset, q, steps, eta, steps)
        qc_log("Gradient descent done")

        derived_clusters <- assignClusters(D, min_d, n_clusters_max)
        derived_clusters
    }

    gradDesc <- function(data_points, D, q, init_steps, eta, steps) {

        vecNorm <- function(v) {
            squared_sum <- sum(v^2)
            if(squared_sum == 0) {
                return(v)
            }
        
            # A new variable would be preferred, but for performance purposes,
            # overwrite in order to avoid reallocating a new vector
            v <- v / sqrt(squared_sum)
            v
        }

        start <- Sys.time()
        dV <- gradient(data_points, q, D)
        dV <- t(apply(dV, MARGIN=1, FUN=vecNorm))
        qc_log(paste0("Gradient descent step: ", (init_steps - steps) + 1,
                      " done")) 
        endtime <- Sys.time()
        qc_log(endtime - start)

        # Faster to overwrite and avoid allocating a new matrix
        D <- D - eta * dV
        eta <- eta * ETA_DECAY

        if (steps <= 1) {
            D
        }
        else {
            gradDesc(data_points, D, q, init_steps, eta, (steps - 1))
        }
    }

    gradient <- function(ri, q, r=ri) {
        n <- nrow(ri)
        p <- ncol(ri)

        V <- matrix(0, n, 1)
        dP2 <- matrix(0, n, 1)
        P <- matrix(0, n, 1)

        dV1 <- matrix(0, n, p)
        dV2 <- matrix(0, n, p)
        dV <- matrix(0, n, p)

        point_iter <- seq(1, n)
        results <- foreach(point=point_iter,
                           .combine=c, .inorder=TRUE) %dopar% {
            calcPoint(point, r, ri, n, q, p)
        }

        # Merge the results from the list returned above
        # (No clean multiple assignment in R working with foreach results)
        result_offset <- 4
        for(point in point_iter) {
            list_index <- point - 1
            curr_index <- list_index * result_offset

            P[point] <- results[[curr_index + 1]]
            dP2[point] <- results[[curr_index + 2]]
            dV1[point, ] <- results[[curr_index + 3]]
            dV2[point, ] <- results[[curr_index + 4]]
        }

        P_zero_ind = which(P == 0)
        if(length(P_zero_ind) > 0) {
            P_non_zero_min <- min(P[-P_zero_ind])
            P[P_zero_ind] <- P_non_zero_min
        }

        V <- -p/2 + q*dP2/P
        E <- -min(V)
        V <- V + E

        for(dim in seq(1, p)) {
            dV[, dim] <- -q*dV1[, dim] + (V - E + (p + 2) / 2) * dV2[, dim]
        }
        dV[P_zero_ind, ] <- 0

        dV
    }

    calcPoint <- function(point, r, ri, n, q, p) {

        repRow <- function(X, row_reps, col_reps){
            p <- length(X)
            matrix(X, row_reps, p, byrow=T)
        }

        empty_col <- matrix(0, n, 1)
        empty_row <- matrix(0, 1, p)
        empty_matrix <- matrix(0, n, p)

        single_laplace <- empty_col
        dV1 <- empty_row
        dV2 <- empty_row

        squared_matrix <- repRow(r[point, ], n, 1)
        squared_matrix <- squared_matrix - ri
        unsquared_mat <- squared_matrix
        squared_matrix <- squared_matrix^2

        D2 <- apply(squared_matrix, MARGIN=1, FUN=sum)
        single_point <- exp(-q * D2)

        for(dim in seq(1, p)) {
            single_laplace <- single_laplace + squared_matrix[, dim] *
                              single_point
        }

        for(dim in seq(1, p)) {
            curr_col <- unsquared_mat[, dim]
            dV1[, dim] <- dV1[, dim] + sum(curr_col * single_laplace)
            dV2[, dim] <- dV2[, dim] + sum(curr_col * single_point)
        }


        P <- sum(single_point)
        dP2 <- sum(single_laplace)

        results <- list(P=P, dP2=dP2, dV1=dV1, dV2=dV2)
        results
    }

    assignClusters <- function(data_points, min_d, n_clusters_max) {

        isAssigned <- function(clusters, index) {
            clusters[index] != 0
        }

        qc_log("Assigning clusters")

        n <- nrow(data_points)
        p <- ncol(data_points)
        clusters <- rep(0, n)
        label <- 1
        start_index = 1

        for (start_index in seq(1, n)) {

            if (min(clusters) > 0) {
                break
            }

            if (!isAssigned(clusters, start_index)) {

                clusters[start_index] <- label
                for(row in seq(start_index, n)) {
                    if (!isAssigned(clusters, row)) {
                        ref_dist <- distance(data_points[start_index, ],
                                             data_points[row, ])
                        if (ref_dist <= min_d) {
                            clusters[row] <- label
                        }
                    }
                }

                label <- label + 1
            }
        }

        clusters
    }

    getMinD <- function(sigma, min_d_factor=2) {
        min_d_factor * sigma
    }

    getEta <- function(sigma) {
        sigma / 2.5
    }

    getQ <- function(sigma) {
        1 / (2 * sigma^2)
    }

    distance <- function(x, y) {
        sqrt(sum((x - y) ^ 2))
    }

    qc_log <- function(out_text) {
        if(verbose) {
            print(out_text)
        }
    }

    orderClusters <- function(training, derived_clusters, n_clusters_max) {
    
        getPopularClusters <- function(derived_clusters, n_clusters_max) {
            popular_clusters_table <- tail(sort(table(derived_clusters)),
                                           n=n_clusters_max)
            popular_clusters <- rev(as.integer(names(popular_clusters_table)))
            popular_clusters
        }

        popular_clusters <- getPopularClusters(derived_clusters,
                                               n_clusters_max)

        clusters <- rep(0, nrow(training))
        i <- 1
        for(label in popular_clusters) {
            clusters[which(derived_clusters == label)] <- i
            i <- i + 1
        }

        clusters
    }

    dataset <- as.matrix(dataset)
    derived_clusters <- main(dataset, sigma, steps, min_d_factor,
                             n_clusters_max)
    ordered_clusters <- orderClusters(dataset, derived_clusters,
                                      n_clusters_max=n_clusters_max)
    ordered_clusters
}
