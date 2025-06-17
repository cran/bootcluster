#' Calculate Stability Measures for a Clustering Method
#'
#' @param x Input data matrix
#' @param k Number of clusters
#' @param scheme Clustering scheme ("kmeans", "hc", or "spectral")
#' @param B Number of bootstrap samples
#' @param n_ref Number of reference distributions
#' @param hc.method Hierarchical clustering method
#' @param dist_method Distance method
#' @return List containing stability measures and clustering results
#' @keywords internal
calculate_stability_measures <- function(x, k, scheme, B=100, n_ref=3,
                                      hc.method="ward.D", dist_method="euclidean") {
    # Input validation
    if (!is.matrix(x) && !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    if (k < 1)
        stop("Number of clusters must be positive")
    if (B < 1)
        stop("Number of bootstrap samples must be positive")
    if (n_ref < 1)
        stop("Number of reference distributions must be positive")
    
    # Calculate base stability with error checking
    tryCatch({
        results <- esmbl.stability(x, k=k, scheme=scheme, B=B,
                                 hc.method=hc.method, dist_method=dist_method)
        
        # Validate results
        if (is.null(results$Smin) || is.na(results$Smin) || is.infinite(results$Smin))
            stop("Invalid Smin value in base stability")
        if (any(is.na(results$obs_wise)) || any(is.infinite(results$obs_wise)))
            stop("Invalid observation-wise stability measures")
        
        #message("Base stability calculation successful")
    }, error = function(e) {
        stop("Error in base stability calculation: ", e$message)
    })
    
    # Calculate reference stability with error checking
    Smin_ref_frame <- matrix(NA, nrow=1, ncol=n_ref)
    for (i_ref in 1:n_ref) {
        tryCatch({
            df_ref <- ref_dist(x)
            results_ref <- esmbl.stability(df_ref, k=k, B=10, scheme=scheme)
            
            # Validate reference results
            if (is.null(results_ref$Smin) || is.na(results_ref$Smin) || is.infinite(results_ref$Smin))
                stop(paste("Invalid Smin value in reference", i_ref))
            
            Smin_ref_frame[,i_ref] <- results_ref$Smin
            #message("Reference stability ", i_ref, " calculation successful")
        }, error = function(e) {
            message("Error in reference ", i_ref, ": ", e$message)
            # Continue with other references
        })
    }
    
    # Calculate final stability delta with validation
    Smin_ref_mean <- mean(Smin_ref_frame, na.rm=TRUE)
    if (is.na(Smin_ref_mean) || is.infinite(Smin_ref_mean))
        stop("Invalid reference stability mean")
    
    Smin_delta <- max(results$Smin - Smin_ref_mean, 0.01)
    message("Final stability delta: ", Smin_delta)
    
    return(list(
        results = results,
        Smin_delta = Smin_delta
    ))
}

#' Define Stability Combination Methods
#'
#' @param alpha Weight for weighted combination
#' @return List of combination functions
#' @keywords internal
define_combination_methods <- function(alpha = 0.5) {
    list(
        product = function(s_delta, s_obs) s_delta * s_obs,
        arithmetic = function(s_delta, s_obs) (s_delta + s_obs) / 2,
        geometric = function(s_delta, s_obs) sqrt(s_delta * s_obs),
        harmonic = function(s_delta, s_obs) 2 / (1/s_delta + 1/s_obs),
        weighted = function(s_delta, s_obs) alpha * s_delta + (1 - alpha) * s_obs
    )
}

#' Create Incidence Matrix for Graph Construction
#'
#' @param results_km K-means results
#' @param results_hc Hierarchical clustering results
#' @param results_sc Spectral clustering results
#' @param k_km Number of k-means clusters
#' @param k_hc Number of hierarchical clusters
#' @param k_sc Number of spectral clusters
#' @param Smin_delta_km K-means stability delta
#' @param Smin_delta_hc Hierarchical clustering stability delta
#' @param Smin_delta_sc Spectral clustering stability delta
#' @param combine_fn Function to combine stability measures
#' @param n_samples Number of samples
#' @return Incidence matrix for graph construction
#' @keywords internal
create_incidence_matrix <- function(results_km, results_hc, results_sc,
                                  k_km, k_hc, k_sc,
                                  Smin_delta_km, Smin_delta_hc, Smin_delta_sc,
                                  combine_fn, n_samples) {
    # Input validation
    validate_stability_results <- function(results, method) {
        if (is.null(results$membership) || is.null(results$obs_wise)) {
            stop(paste(method, "results missing required components"))
        }
        if (any(is.na(results$obs_wise)) || any(is.infinite(results$obs_wise))) {
            stop(paste(method, "stability measures contain NA/Inf values"))
        }
        if (any(results$obs_wise < 0)) {
            stop(paste(method, "stability measures contain negative values"))
        }
    }
    
    validate_stability_results(results_km, "K-means")
    validate_stability_results(results_hc, "Hierarchical")
    validate_stability_results(results_sc, "Spectral")
    
    # Initialize matrix
    inc <- matrix(0, nrow=sum(k_km, k_hc, k_sc), ncol=n_samples)
    
    # Helper function to safely compute weights
    safe_combine <- function(obs_wise, s_delta, method) {
        tryCatch({
            weights <- combine_fn(obs_wise, s_delta)
            if (any(is.na(weights)) || any(is.infinite(weights))) {
                stop(paste("Invalid weights in", method))
            }
            weights
        }, error = function(e) {
            stop(paste("Error computing weights for", method, ":", e$message))
        })
    }
    
    # Add k-means edges
    for (i in 1:k_km) {
        clust_mem <- which(results_km$membership == i)
        if (length(clust_mem) > 0) {
            weights <- safe_combine(results_km$obs_wise[clust_mem], 
                                  Smin_delta_km, "k-means")
            inc[i, clust_mem] <- weights
        }
    }
    
    # Add hierarchical clustering edges
    for (i in 1:k_hc) {
        clust_mem <- which(results_hc$membership == i)
        if (length(clust_mem) > 0) {
            weights <- safe_combine(results_hc$obs_wise[clust_mem], 
                                  Smin_delta_hc, "hierarchical")
            inc[i + k_km, clust_mem] <- weights
        }
    }
    
    # Add spectral clustering edges
    for (i in 1:k_sc) {
        clust_mem <- which(results_sc$membership == i)
        if (length(clust_mem) > 0) {
            weights <- safe_combine(results_sc$obs_wise[clust_mem], 
                                  Smin_delta_sc, "spectral")
            inc[i + k_km + k_hc, clust_mem] <- weights
        }
    }
    
    # Remove empty rows
    inc <- inc[rowSums(inc) > 0, , drop=FALSE]
    return(inc)
}

#' Create Graph and Find Communities Using Different Methods
#'
#' @param inc Incidence matrix
#' @param method Community detection method: "fastgreedy", "metis", "hmetis", or "all" (default: "fastgreedy")
#' @return List containing graph and community detection results from specified method(s)
#' @keywords internal
create_graph_and_communities <- function(inc, method = "fastgreedy") {
    if (nrow(inc) == 0 || ncol(inc) == 0) {
        stop("Invalid incidence matrix: no valid clustering structure found")
    }
    
    # Create graph
    graph <- igraph::graph_from_biadjacency_matrix(inc, weighted=TRUE)
    n_vertices <- nrow(inc)
    
    # Helper function to extract membership
    extract_membership <- function(communities) {
        igraph::membership(communities)[(n_vertices + 1):(ncol(inc) + n_vertices)]
    }
    
    # Function to get fastgreedy results with timing
    get_fastgreedy <- function() {
        start_time <- Sys.time()
        communities <- igraph::cluster_fast_greedy(graph, weights=igraph::E(graph)$weight)
        end_time <- Sys.time()
        runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
        
        list(
            communities = communities,
            membership = extract_membership(communities),
            runtime = runtime
        )
    }
    
    # Function to get METIS results with timing
    get_metis <- function() {
        start_time <- Sys.time()
        # Use leading eigenvector method which can handle disconnected graphs
        communities <- igraph::cluster_leading_eigen(graph, weights=igraph::E(graph)$weight)
        
        # If the method fails, fall back to fastgreedy
        if (is.null(communities)) {
            communities <- igraph::cluster_fast_greedy(graph, weights=igraph::E(graph)$weight)
        }
        end_time <- Sys.time()
        runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
        
        list(
            communities = communities,
            membership = extract_membership(communities),
            runtime = runtime
        )
    }
    
    # Function to get hMETIS results with timing
    get_hmetis <- function() {
        start_time <- Sys.time()
        # hMETIS is specifically designed for hypergraphs
        # We'll use a multilevel approach similar to hMETIS
        communities <- igraph::cluster_louvain(graph, weights=igraph::E(graph)$weight)
        end_time <- Sys.time()
        runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
        
        list(
            communities = communities,
            membership = extract_membership(communities),
            runtime = runtime
        )
    }
    
    # Get results based on method
    results <- if (method == "all") {
        list(
            fastgreedy = get_fastgreedy(),
            metis = get_metis(),
            hmetis = get_hmetis(),
            graph = graph
        )
    } else {
        # Only run the requested method
        result <- switch(method,
            "fastgreedy" = get_fastgreedy(),
            "metis" = get_metis(),
            "hmetis" = get_hmetis(),
            stop("Invalid method. Choose 'fastgreedy', 'metis', 'hmetis', or 'all'")
        )
        result$graph <- graph
        result
    }
    
    return(results)
}

#' Calculate Comparison Statistics
#'
#' @param results List of results for each method
#' @param method_names Names of combination methods
#' @return List of comparison statistics
#' @keywords internal
calculate_comparison_stats <- function(results, method_names) {
    list(
        k_consensus = sapply(results[method_names], function(x) x$k_consensus),
        modularity = sapply(results[method_names], function(x) {
            comm <- igraph::cluster_fast_greedy(x$graph)
            membership <- as.numeric(igraph::membership(comm))
            igraph::modularity(x$graph, membership)
        }),
        edge_weight_summary = lapply(results[method_names], function(x) {
            list(
                min = min(x$edge_weights),
                max = max(x$edge_weights),
                mean = mean(x$edge_weights),
                sd = sd(x$edge_weights)
            )
        })
    )
}
