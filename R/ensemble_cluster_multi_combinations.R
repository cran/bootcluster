#' Multi-Method Ensemble Clustering with Multiple Stability Combinations
#'
#' @description Implements ensemble clustering using multiple methods for combining stability measures,
#' generating separate consensus results for each combination method.
#'
#' @param x data.frame or matrix where rows are observations and columns are features
#' @param k_km number of clusters for k-means clustering
#' @param k_hc number of clusters for hierarchical clustering
#' @param k_sc number of clusters for spectral clustering
#' @param n_ref number of reference distributions for stability assessment (default: 3)
#' @param B number of bootstrap samples for stability estimation (default: 100)
#' @param hc.method hierarchical clustering method (default: "ward.D")
#' @param dist_method distance method for spectral clustering (default: "euclidean")
#' @param alpha weight for weighted combination (default: 0.5)
#'
#' @return A list containing results for each combination method:
#' \describe{
#'   \item{product}{Results using product combination}
#'   \item{arithmetic}{Results using arithmetic mean}
#'   \item{geometric}{Results using geometric mean}
#'   \item{harmonic}{Results using harmonic mean}
#'   \item{weighted}{Results using weighted combination}
#' }
#' Each method's results contain:
#' \describe{
#'   \item{fastgreedy}{Results from fast greedy community detection}
#'   \item{metis}{Results from METIS (leading eigenvector) community detection}
#'   \item{hmetis}{Results from hMETIS (Louvain) community detection}
#'   \item{graph}{igraph object of the ensemble graph}
#'   \item{edge_weights}{Edge weights of the graph}
#'   \item{individual_results}{Results from individual clustering methods}
#'   \item{stability_measures}{Stability measures}
#'   \item{incidence_matrix}{Incidence matrix used for graph construction}
#' }
#' Each community detection method's results contain:
#' \describe{
#'   \item{membership}{Final cluster assignments}
#'   \item{k_consensus}{Number of clusters found}
#' }
#' The function also returns comparison statistics for each community detection method:
#' \describe{
#'   \item{comparison$fastgreedy}{Comparison stats for fast greedy results}
#'   \item{comparison$metis}{Comparison stats for METIS results}
#'   \item{comparison$hmetis}{Comparison stats for hMETIS results}
#' }
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
#' results <- ensemble_cluster_multi_combinations(df, k_km=3, k_hc=3, k_sc=3)
#' # Compare cluster assignments from different methods
#' table(product = results$product$membership, 
#'       arithmetic = results$arithmetic$membership)
#' }
#'
#' @export
ensemble_cluster_multi_combinations <- function(x, k_km, k_hc, k_sc, n_ref=3, B=100, 
                                             hc.method="ward.D", dist_method="euclidean",
                                             alpha=0.25) {
    # Input validation
    if (!is.matrix(x) && !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    if (k_km < 1 || k_hc < 1 || k_sc < 1)
        stop("Number of clusters must be positive")
    
    # Start timing for common operations
    common_start_time <- Sys.time()
    
    # Define combination methods
    combine_methods <- define_combination_methods(alpha)
    
    # Track runtime
    runtime <- list()
    
    # Calculate stability measures for each clustering method
    message("Start calculating k-means stability measures")
    km_stability <- calculate_stability_measures(x, k=k_km, scheme="kmeans", 
                                              B=B, n_ref=n_ref)
    
    message("Start calculating hclust stability measures")
    hc_stability <- calculate_stability_measures(x, k=k_hc, scheme="hc", 
                                              B=B, n_ref=n_ref,
                                              hc.method=hc.method, 
                                              dist_method=dist_method)
    
    message("Start calculating spectral stability measures")
    sc_stability <- calculate_stability_measures(x, k=k_sc, scheme="spectral", 
                                              B=B, n_ref=n_ref)
    
    # Extract results and stability measures
    results_km <- km_stability$results
    results_hc <- hc_stability$results
    results_sc <- sc_stability$results
    
    Smin_delta_km <- km_stability$Smin_delta
    Smin_delta_hc <- hc_stability$Smin_delta
    Smin_delta_sc <- sc_stability$Smin_delta
    
    # End timing for common operations
    common_end_time <- Sys.time()
    common_runtime <- as.numeric(difftime(common_end_time, common_start_time, units="secs"))
    message(paste("Common operations completed in", round(common_runtime, 2), "seconds"))
    
    # Generate results for each combination method
    results <- list()
    
    for (method_name in names(combine_methods)) {
        message(paste("Processing", method_name, "combination method"))
        combine_fn <- combine_methods[[method_name]]
        
        # Create incidence matrix with error checking
        tryCatch({
            message("Creating incidence matrix...")
            inc <- create_incidence_matrix(
                results_km, results_hc, results_sc,
                k_km, k_hc, k_sc,
                Smin_delta_km, Smin_delta_hc, Smin_delta_sc,
                combine_fn, nrow(x)
            )
            message("Incidence matrix created successfully")
            
            # Check for invalid values
            if (any(is.na(inc)) || any(is.infinite(inc))) {
                stop("Invalid values in incidence matrix: NAs or Infinities detected")
            }
            if (any(inc < 0)) {
                stop("Invalid values in incidence matrix: Negative values detected")
            }
            message("Incidence matrix validation passed")
            
        }, error = function(e) {
            stop("Error creating incidence matrix: ", e$message)
        })
        
        # Skip if no valid clustering structure
        if (nrow(inc) == 0 || ncol(inc) == 0) {
            warning(paste("No valid clustering structure found for method:", method_name))
            next
        }
        
        # Create graph and find communities with error checking
        tryCatch({
            message("Creating graph and finding communities...")
            graph_results <- create_graph_and_communities(inc, method = "all")
            message("Graph creation and community detection successful")
            
            # Validate graph results
            if (is.null(graph_results$graph)) {
                stop("Graph creation failed: NULL graph returned")
            }
            if (is.null(graph_results$fastgreedy) || is.null(graph_results$metis) || is.null(graph_results$hmetis)) {
                stop("One or more community detection methods failed")
            }
            message("Graph validation passed")
            
        }, error = function(e) {
            stop("Error in graph creation or community detection: ", e$message)
        })
        
        # Store results for each community detection method
        method_results <- list(
            fastgreedy = list(
                membership = graph_results$fastgreedy$membership,
                k_consensus = length(unique(graph_results$fastgreedy$membership)),
                graph = graph_results$fastgreedy$communities,
                runtime = common_runtime + graph_results$fastgreedy$runtime
            ),
            metis = list(
                membership = graph_results$metis$membership,
                k_consensus = length(unique(graph_results$metis$membership)),
                graph = graph_results$metis$communities,
                runtime = common_runtime + graph_results$metis$runtime
            ),
            hmetis = list(
                membership = graph_results$hmetis$membership,
                k_consensus = length(unique(graph_results$hmetis$membership)),
                graph = graph_results$hmetis$communities,
                runtime = common_runtime + graph_results$hmetis$runtime
            )
        )
        
        # Store common graph and edge weights
        method_results$graph <- graph_results$graph
        method_results$edge_weights <- igraph::E(graph_results$graph)$weight
        
        # Add common results
        method_results$individual_results <- list(
            kmeans = list(
                membership = results_km$membership,
                runtime = results_km$runtime
            ),
            hierarchical = list(
                membership = results_hc$membership,
                runtime = results_hc$runtime
            ),
            spectral = list(
                membership = results_sc$membership,
                runtime = results_sc$runtime
            )
        )
        
        method_results$stability_measures <- list(
            kmeans_delta = Smin_delta_km,
            hierarchical_delta = Smin_delta_hc,
            spectral_delta = Smin_delta_sc
        )
        
        method_results$incidence_matrix <- inc
        method_results$common_runtime <- common_runtime
        
        # Store all results for this combination method
        results[[method_name]] <- method_results
    }
    
    # Add comparison statistics
    results$comparison <- list()
    for (comm_method in c("fastgreedy", "metis", "hmetis")) {
        results$comparison[[comm_method]] <- calculate_comparison_stats(
            lapply(results, function(x) {
                list(
                    membership = x[[comm_method]]$membership,
                    k_consensus = x[[comm_method]]$k_consensus,
                    graph = x$graph,
                    edge_weights = x$edge_weights
                )
            }),
            names(combine_methods)
        )
    }
    
    # Add runtime information
    results$runtime <- list(
        kmeans = results_km$runtime,
        hierarchical = results_hc$runtime,
        spectral = results_sc$runtime,
        common = common_runtime
    )
    
    return(results)
}
