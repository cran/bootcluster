#' @name ensemble.cluster.multi
#' @title Multi-Method Ensemble Clustering with Graph-based Consensus
#'
#' @description Implements ensemble clustering by combining multiple clustering methods 
#' (k-means, hierarchical, and spectral clustering) using a graph-based consensus approach.
#'
#' @param x data.frame or matrix where rows are observations and columns are features
#' @param k_km number of clusters for k-means clustering
#' @param k_hc number of clusters for hierarchical clustering
#' @param k_sc number of clusters for spectral clustering
#' @param n_ref number of reference distributions for stability assessment (default: 3)
#' @param B number of bootstrap samples for stability estimation (default: 100)
#' @param hc.method hierarchical clustering method (default: "ward.D")
#' @param dist_method distance method for spectral clustering (default: "euclidean")
#'
#' @return A list containing:
#' \describe{
#'   \item{membership}{Final cluster assignments from ensemble consensus}
#'   \item{k_consensus}{Number of clusters found in consensus}
#'   \item{individual_results}{List of results from individual clustering methods}
#'   \item{stability_measures}{Stability measures for each method}
#'   \item{graph}{igraph object of the ensemble graph}
#' }
#'
#' @details
#' This function implements a multi-method ensemble clustering approach that:
#' 1. Applies multiple clustering methods (k-means, hierarchical, spectral)
#' 2. Assesses stability of each clustering through bootstrapping
#' 3. Constructs a weighted bipartite graph representing all clusterings
#' 4. Uses fast greedy community detection for final consensus
#'
#' @importFrom kernlab specc
#' @importFrom stats hclust cutree dist cor var runif prcomp kmeans
#' @importFrom flexclust dist2
#' @importFrom network set.vertex.attribute
#' @importFrom grid unit.c
#' @importFrom dplyr left_join
#'
#' @examples
#' \donttest{
#' data(iris)
#' df <- iris[,1:4]
#' result <- ensemble.cluster.multi(df, k_km=3, k_hc=3, k_sc=3)
#' plot(df[,1:2], col=result$membership, pch=16)
#' }
#'
#' @export
ensemble.cluster.multi <- function(x, k_km, k_hc, k_sc, n_ref=3, B=100, hc.method="ward.D", dist_method="euclidean") {
    # Input validation
    if (!is.matrix(x) && !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    if (k_km < 1 || k_hc < 1 || k_sc < 1)
        stop("Number of clusters must be positive")
    
    result <- list()
    
    # K-means clustering and stability
    results_km <- esmbl.stability(x, k_km, scheme="kmeans", B=B)
    Smin_ref_km_frame <- matrix(NA, nrow=1, ncol=n_ref)
    for (i_ref in 1:n_ref) {
        df_ref <- ref_dist(x)
        results_km_ref <- esmbl.stability(df_ref, k_km, B=50, scheme="kmeans")
        Smin_ref_km_frame[,i_ref] <- results_km_ref$Smin
    }
    Smin_ref_km_mean <- mean(Smin_ref_km_frame, na.rm=TRUE)
    Smin_delta_km <- max(results_km$Smin - Smin_ref_km_mean, 0.01)
    
    # Hierarchical clustering and stability
    results_hc <- esmbl.stability(x, k_hc, scheme="hc", B=B, hc.method=hc.method, dist_method = dist_method)
    Smin_ref_hc_frame <- matrix(NA, nrow=1, ncol=n_ref)
    for (i_ref in 1:n_ref) {
        df_ref <- ref_dist(x)
        results_hc_ref <- esmbl.stability(df_ref, k_hc, B=50, scheme="hc")
        Smin_ref_hc_frame[,i_ref] <- results_hc_ref$Smin
    }
    Smin_ref_hc_mean <- mean(Smin_ref_hc_frame, na.rm=TRUE)
    Smin_delta_hc <- max(results_hc$Smin - Smin_ref_hc_mean, 0.01)
    
    # Spectral clustering and stability
    results_sc <- esmbl.stability(x, k_sc, scheme="spectral", B=B)
    Smin_ref_sc_frame <- matrix(NA, nrow=1, ncol=n_ref)
    for (i_ref in 1:n_ref) {
        df_ref <- ref_dist(x)
        results_sc_ref <- esmbl.stability(df_ref, k_sc, B=50, scheme="spectral")
        Smin_ref_sc_frame[,i_ref] <- results_sc_ref$Smin
    }
    Smin_ref_sc_mean <- mean(Smin_ref_sc_frame, na.rm=TRUE)
    Smin_delta_sc <- max(results_sc$Smin - Smin_ref_sc_mean, 0.01)
    
    # Create incidence matrix for graph
    inc <- matrix(0, nrow=sum(k_km, k_hc, k_sc), ncol=nrow(x))
    
    # Add k-means edges
    for (i in 1:k_km) {
        clust_mem <- which(results_km$membership == i)
        if (length(clust_mem) > 0) {
            weights <- results_km$obs_wise[clust_mem] * Smin_delta_km
            inc[i, clust_mem] <- weights
        }
    }
    
    # Add hierarchical clustering edges
    for (i in 1:k_hc) {
        clust_mem <- which(results_hc$membership == i)
        if (length(clust_mem) > 0) {
            weights <- results_hc$obs_wise[clust_mem] * Smin_delta_hc
            inc[i + k_km, clust_mem] <- weights
        }
    }
    
    # Add spectral clustering edges
    for (i in 1:k_sc) {
        clust_mem <- which(results_sc$membership == i)
        if (length(clust_mem) > 0) {
            weights <- results_sc$obs_wise[clust_mem] * Smin_delta_sc
            inc[i + k_km + k_hc, clust_mem] <- weights
        }
    }
    
    # Remove any empty rows
    inc <- inc[rowSums(inc) > 0, , drop=FALSE]
    
    if (nrow(inc) == 0 || ncol(inc) == 0) {
        stop("No valid clustering structure found")
    }
    
    # Create graph and find communities
    graph <- igraph::graph_from_biadjacency_matrix(inc, weighted=TRUE)
    communities <- igraph::cluster_fast_greedy(graph, weights=igraph::E(graph)$weight)
    
    # Extract final clustering
    n_vertices <- nrow(inc)
    consensus_membership <- igraph::membership(communities)[(n_vertices + 1):(ncol(inc) + n_vertices)]
    
    # Prepare return values
    result$membership <- consensus_membership
    result$k_consensus <- length(unique(consensus_membership))
    result$individual_results <- list(
        kmeans = results_km,
        hierarchical = results_hc,
        spectral = results_sc
    )
    result$stability_measures <- list(
        kmeans_delta = Smin_delta_km,
        hierarchical_delta = Smin_delta_hc,
        spectral_delta = Smin_delta_sc
    )
    result$graph <- graph
    
    return(result)
}
