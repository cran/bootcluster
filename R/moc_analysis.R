#' Multi-Method Ensemble Clustering Analysis for Multiple-Objective Clustering (MOC) Datasets
#'
#' @description Performs ensemble clustering analysis on multiple datasets using
#' different clustering methods and compares their performance.
#'
#' @param datasets List of datasets to analyze
#' @param selected Indices of datasets to analyze
#' @param n_ref Number of reference distributions (default: 3)
#' @param B Number of bootstrap samples (default: 100)
#' @param plot Whether to generate plots (default: TRUE)
#' @param plot_file Output file for plots (default: NULL)
#'
#' @return A list containing:
#' \describe{
#'   \item{results}{Results for each dataset}
#'   \item{ari_table}{Adjusted Rand Index comparison table}
#'   \item{runtime_table}{Runtime comparison table}
#'   \item{plots}{List of generated plots if plot=TRUE}
#' }
#'
#' @import mclust
#' @importFrom stats prcomp
#' @importFrom grDevices rainbow recordPlot pdf dev.off dev.new
#' @importFrom graphics par plot title boxplot
#' @importFrom progress progress_bar
#'
#' @export
analyze_moc_datasets <- function(datasets, selected, n_ref=3, B=100, 
                               plot=TRUE, plot_file=NULL) {
    results <- list()
    plots <- list()
    
    # Setup plotting if requested
    if (plot) {
        if (!is.null(plot_file)) {
            pdf(file=plot_file, width=15, height=9)
        }
        # Each row represents a dataset, with 3 columns for different methods
        n_rows <- length(selected)
        par(mfrow=c(n_rows, 3))
    }
    
    # Setup progress bar
    message("Starting analysis of ", length(selected), " datasets...")
    pb <- progress::progress_bar$new(
        format = "  Processing :what [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
        total = length(selected),
        clear = FALSE,
        width = 100,
        show_after = 0
    )
    
    # Process each dataset
    for (i in selected) {
        current_dataset <- names(datasets)[i]
        message("\nAnalyzing dataset: ", current_dataset)
        pb$tick(tokens = list(what = current_dataset))
        df <- datasets[[i]]
        true_clusters = df[,1]
        k <- max(true_clusters)  # number of clusters from first column
        df <- df[,-1]     # remove cluster column
        scaled_df <- scale(df)
        
        # Create PCA plot data
        if (ncol(df) > 2) {
            pca <- prcomp(df, scale.=TRUE)
            pca_plot <- pca$x[,1:2]
        } else {
            pca_plot <- df
        }
        
        # Run ensemble clustering with all methods
        tryCatch({
            tryCatch({
                # If stability measures work, try the ensemble
                result <- ensemble_cluster_multi_combinations(
                    scaled_df, 
                    k_km=k, k_hc=k, k_sc=k,
                    n_ref=n_ref, B=B,
                    hc.method="centroid"
                )
            }, error = function(e) {
                message("Error in stability calculation: ", e$message)
                stop("Stability calculation failed")
            })
            
            # Store results
            results[[names(datasets)[i]]] <- list(
                ensemble_results = result,
                pca_plot = pca_plot,
                true_clusters = true_clusters
            )
            
            # Generate plots if requested
            if (plot) {
                # Plot fastgreedy results
                plot(pca_plot,
                     col = rainbow(k)[result$product$fastgreedy$membership],
                     pch = result$product$fastgreedy$membership,
                     main = paste(names(datasets)[i], "k =", k, "fg"))
                
                # Plot METIS results
                plot(pca_plot,
                     col = rainbow(k)[result$product$metis$membership],
                     pch = result$product$metis$membership,
                     main = paste(names(datasets)[i], "k =", k, "metis"))
                
                # Plot hMETIS results
                plot(pca_plot,
                     col = rainbow(k)[result$product$hmetis$membership],
                     pch = result$product$hmetis$membership,
                     main = paste(names(datasets)[i], "k =", k, "hmetis"))
                
                # Store plots
                plots[[names(datasets)[i]]] <- recordPlot()
            }
        }, error = function(e) {
            message("\nError processing dataset ", names(datasets)[i], ": ", e$message)
            # Store partial results
            results[[names(datasets)[i]]] <- list(
                error = e$message,
                pca_plot = pca_plot,
                true_clusters = true_clusters
            )
        }, warning = function(w) {
            message("\nWarning in dataset ", names(datasets)[i], ": ", w$message)
        })
        
    }
    
    if (plot && !is.null(plot_file)) {
        dev.off()
    }
    
    # Calculate ARI table for successful results
    ari_table <- do.call(rbind, lapply(names(results), function(name) {
        res <- results[[name]]
        if (!is.null(res$error)) {
            # Create empty rows for each combination method
            methods <- c("product", "arithmetic", "geometric", "harmonic", "weighted")
            do.call(rbind, lapply(methods, function(method) {
                data.frame(
                    Dataset = paste(name, method, sep="."),
                    KMeans = NA,
                    Hierarchical = NA,
                    Spectral = NA,
                    Fastgreedy = NA,
                    METIS = NA,
                    hMETIS = NA
                )
            }))
        } else {
            true_clusters <- res$true_clusters
            # Create a row for each combination method
            methods <- c("product", "arithmetic", "geometric", "harmonic", "weighted")
            do.call(rbind, lapply(methods, function(method) {
                data.frame(
                    Dataset = paste(name, method, sep="."),
                    KMeans = mclust::adjustedRandIndex(res$ensemble_results[[method]]$individual_results$kmeans$membership, 
                                                      true_clusters),
                    Hierarchical = mclust::adjustedRandIndex(res$ensemble_results[[method]]$individual_results$hierarchical$membership, 
                                                           true_clusters),
                    Spectral = mclust::adjustedRandIndex(res$ensemble_results[[method]]$individual_results$spectral$membership, 
                                                        true_clusters),
                    Fastgreedy = mclust::adjustedRandIndex(res$ensemble_results[[method]]$fastgreedy$membership, 
                                                          true_clusters),
                    METIS = mclust::adjustedRandIndex(res$ensemble_results[[method]]$metis$membership, 
                                                     true_clusters),
                    hMETIS = mclust::adjustedRandIndex(res$ensemble_results[[method]]$hmetis$membership, 
                                                      true_clusters)
                )
            }))
        }
    }))
    
    # Calculate runtime table for successful results
    runtime_table <- do.call(rbind, lapply(names(results), function(name) {
        res <- results[[name]]
        if (!is.null(res$error)) {
            # Create empty rows for each combination method
            methods <- c("product", "arithmetic", "geometric", "harmonic", "weighted")
            do.call(rbind, lapply(methods, function(method) {
                data.frame(
                    Dataset = paste(name, method, sep="."),
                    KMeans = NA,
                    Hierarchical = NA,
                    Spectral = NA,
                    Fastgreedy = NA,
                    METIS = NA,
                    hMETIS = NA
                )
            }))
        } else {
            # Create a row for each combination method
            methods <- c("product", "arithmetic", "geometric", "harmonic", "weighted")
            do.call(rbind, lapply(methods, function(method) {
                # Get individual method runtimes
                kmeans_runtime <- res$ensemble_results[[method]]$individual_results$kmeans$runtime
                hierarchical_runtime <- res$ensemble_results[[method]]$individual_results$hierarchical$runtime
                spectral_runtime <- res$ensemble_results[[method]]$individual_results$spectral$runtime
                
                # Get community detection method runtimes (includes common runtime)
                fastgreedy_runtime <- res$ensemble_results[[method]]$fastgreedy$runtime
                metis_runtime <- res$ensemble_results[[method]]$metis$runtime
                hmetis_runtime <- res$ensemble_results[[method]]$hmetis$runtime
                
                data.frame(
                    Dataset = paste(name, method, sep="."),
                    KMeans = kmeans_runtime,
                    Hierarchical = hierarchical_runtime,
                    Spectral = spectral_runtime,
                    Fastgreedy = fastgreedy_runtime,
                    METIS = metis_runtime,
                    hMETIS = hmetis_runtime
                )
            }))
        }
    }))
    
    return(list(
        results = results,
        ari_table = ari_table,
        runtime_table = runtime_table,
        plots = if(plot) plots else NULL
    ))
}

#' Plot MOC Results
#'
#' @param results Results from analyze_moc_datasets
#' @param dataset_names Name or vector of names of datasets to plot
#' @param methods Methods to plot (default: c("fastgreedy", "metis", "hmetis"))
#' @param plot_file Output file for plots (default: NULL)
#' @param max_plots_per_page Maximum number of plots per page (default: 12)
#' @param mar Margins for plots (default: c(2, 2, 2, 1))
#'
#' @export
plot_moc_results <- function(results, dataset_names, 
                           methods=c("kmeans", "hierarchical", "spectral", 
                                   "fastgreedy", "metis", "hmetis"),
                           plot_file=NULL, max_plots_per_page=12, 
                           mar=c(2, 2, 2, 1)) {
    # Convert single dataset name to vector if needed
    if (!is.vector(dataset_names) || is.character(dataset_names) && length(dataset_names) == 1) {
        dataset_names <- c(dataset_names)
    }
    
    # Filter out any dataset names that don't exist in results
    valid_datasets <- names(results$results)
    dataset_names <- dataset_names[dataset_names %in% valid_datasets]
    
    if (length(dataset_names) == 0) {
        stop("No valid dataset names provided")
    }
    
    # Open PDF device if plot_file is provided
    if (!is.null(plot_file)) {
        pdf(file=plot_file, width=15, height=10)
    }
    
    # Calculate total number of plots
    n_methods <- length(methods)
    total_plots <- length(dataset_names) * n_methods
    
    # Calculate number of pages needed
    n_pages <- ceiling(total_plots / max_plots_per_page)
    
    # Calculate plots per page
    plots_per_page <- min(max_plots_per_page, total_plots)
    
    # Calculate rows and columns for layout
    if (n_methods <= 3) {
        # If few methods, use them as columns
        n_cols <- n_methods
        n_rows <- min(ceiling(plots_per_page / n_cols), ceiling(length(dataset_names)))
    } else {
        # Otherwise, use a more balanced layout
        n_cols <- min(3, n_methods)
        n_rows <- min(ceiling(plots_per_page / n_cols), ceiling(length(dataset_names) * n_methods / n_cols))
    }
    
    # Set smaller margins to fit more plots
    old_par <- par(mfrow=c(n_rows, n_cols), mar=mar)
    on.exit(par(old_par), add=TRUE)
    
    # Counter for plots
    plot_count <- 0
    page_count <- 1
    
    # Plot each dataset
    for (dataset_name in dataset_names) {
        result <- results$results[[dataset_name]]
        if (is.null(result) || !is.null(result$error)) {
            message("Skipping dataset ", dataset_name, " due to errors")
            next
        }
        
        pca_plot <- result$pca_plot
        k <- length(unique(result$true_clusters))
        
        for (method in methods) {
            # Check if we need a new page
            if (plot_count >= plots_per_page && plot_count < total_plots) {
                if (!is.null(plot_file)) {
                    dev.new()
                } else {
                    # For interactive plotting, just create a new plot window
                    dev.new(width=15, height=10)
                    par(mfrow=c(n_rows, n_cols), mar=mar)
                }
                plot_count <- 0
                page_count <- page_count + 1
            }
            
            # Get membership data
            if (method %in% c("kmeans", "hierarchical", "spectral")) {
                membership <- result$ensemble_results$product$individual_results[[method]]$membership
            } else {
                membership <- result$ensemble_results$product[[method]]$membership
            }
            
            # Create the plot
            plot(pca_plot,
                 col = rainbow(k)[membership],
                 pch = membership,
                 main = paste(dataset_name, "k =", k, method))
            
            # Increment plot counter
            plot_count <- plot_count + 1
        }
    }
    
    # Close PDF device if it was opened
    if (!is.null(plot_file)) {
        dev.off()
    }
    
    # Return information about the plots
    invisible(list(
        datasets = dataset_names,
        methods = methods,
        pages = page_count
    ))
}

#' Compare MOC Results
#'
#' @param results Results from analyze_moc_datasets
#' @param metric Metric to compare ("ari", "runtime", or "modularity")
#' @param plot Whether to generate comparison plot (default: TRUE)
#'
#' @export
compare_moc_results <- function(results, metric="ari", plot=TRUE) {
    if (metric == "ari") {
        comparison <- results$ari_table
        
        # Split Dataset column into Dataset and Method
        split_names <- strsplit(as.character(comparison$Dataset), "\\.")
        comparison$Dataset <- sapply(split_names, "[", 1)
        comparison$Method <- sapply(split_names, "[", 2)
        
        # Reorder columns
        comparison <- comparison[, c("Dataset", "Method", "KMeans", "Hierarchical", 
                                   "Spectral", "Fastgreedy", "METIS", "hMETIS")]
        
        if (plot) {
            # Create separate boxplots for each combination method
            methods <- unique(comparison$Method)
            n_methods <- length(methods)
            old_par <- par(mfrow=c(2, ceiling(n_methods/2)), mar=c(8,4,4,2))
            on.exit(par(old_par))
            
            for (method in methods) {
                subset_data <- comparison[comparison$Method == method, 
                                       c("KMeans", "Hierarchical", "Spectral", 
                                         "Fastgreedy", "METIS", "hMETIS")]
                boxplot(subset_data, 
                       main=paste("Method:", method),
                       ylab="ARI",
                       las=2)  # Rotate x-axis labels
            }
        }
    } else if (metric == "runtime") {
        comparison <- results$runtime_table
        
        # Split Dataset column into Dataset and Method
        split_names <- strsplit(as.character(comparison$Dataset), "\\.")
        comparison$Dataset <- sapply(split_names, "[", 1)
        comparison$Method <- sapply(split_names, "[", 2)
        
        # Reorder columns
        comparison <- comparison[, c("Dataset", "Method", "KMeans", "Hierarchical", 
                                   "Spectral", "Fastgreedy", "METIS", "hMETIS")]
        
        if (plot) {
            # Create separate boxplots for each combination method
            methods <- unique(comparison$Method)
            n_methods <- length(methods)
            old_par <- par(mfrow=c(2, ceiling(n_methods/2)), mar=c(8,4,4,2))
            on.exit(par(old_par))
            
            for (method in methods) {
                subset_data <- comparison[comparison$Method == method, 
                                       c("KMeans", "Hierarchical", "Spectral", 
                                         "Fastgreedy", "METIS", "hMETIS")]
                boxplot(subset_data, 
                       main=paste("Method:", method),
                       ylab="Runtime (seconds)",
                       las=2)  # Rotate x-axis labels
            }
        }
    } else if (metric == "modularity") {
        comparison <- do.call(rbind, lapply(names(results$results), function(name) {
            res <- results$results[[name]]
            if (!is.null(res$error)) return(NULL)
            
            # Get modularity for each combination method
            methods <- c("product", "arithmetic", "geometric", "harmonic", "weighted")
            do.call(rbind, lapply(methods, function(method) {
                data.frame(
                    Dataset = paste(name, method, sep="."),
                    Fastgreedy = res$ensemble_results$comparison$fastgreedy$modularity,
                    METIS = res$ensemble_results$comparison$metis$modularity,
                    hMETIS = res$ensemble_results$comparison$hmetis$modularity
                )
            }))
        }))
        
        if (plot) {
            # Split Dataset column into Dataset and Method
            split_names <- strsplit(as.character(comparison$Dataset), "\\.")
            comparison$Dataset <- sapply(split_names, "[", 1)
            comparison$Method <- sapply(split_names, "[", 2)
            
            # Create separate boxplots for each combination method
            methods <- unique(comparison$Method)
            n_methods <- length(methods)
            old_par <- par(mfrow=c(2, ceiling(n_methods/2)), mar=c(8,4,4,2))
            on.exit(par(old_par))
            
            for (method in methods) {
                subset_data <- comparison[comparison$Method == method, 
                                       c("Fastgreedy", "METIS", "hMETIS")]
                boxplot(subset_data,
                       main=paste("Method:", method),
                       ylab="Modularity",
                       las=2)  # Rotate x-axis labels
            }
        }
    }
    
    return(comparison)
}
