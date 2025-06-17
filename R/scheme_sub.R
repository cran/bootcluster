#' Subsampling-based K-means Clustering
#'
#' @param df data.frame or matrix of the dataset
#' @param nk number of clusters
#' @param B number of bootstrap samples
#' @param cut_ratio ratio for subsampling (default: 0.5)
#' @param r number of k-means runs (default: 5)
#' @return List containing clustering results and stability measures
#' @importFrom stats kmeans
#' @keywords internal
scheme_sub_km <- function(df, nk, B=10, cut_ratio=0.5, r=5) {
    results <- list()
    n <- nrow(df)
    clst.mat <- matrix(NA, nrow=B+1, ncol=n)
    
    # Initial k-means clustering with timing and error handling
    start_time <- Sys.time()
    tryCatch({
        # Try with default settings first
        km <- stats::kmeans(df, centers=nk, nstart=r, iter.max=300)
    }, warning = function(w) {
        # If we get a warning about iterations, try with more iterations
        if (grepl("Quick-TRANSfer", w$message) || grepl("did not converge", w$message)) {
            message("Warning in k-means: ", w$message, ". Trying with increased iterations...")
            km <<- stats::kmeans(df, centers=nk, nstart=r+5, iter.max=500, algorithm="Lloyd")
        } else {
            # For other warnings, just log them
            message("Warning in k-means: ", w$message)
        }
    })
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
    
    clst <- km$cluster
    clst.mat[1,] <- clst
    
    # Store subsamples
    subsamples.store <- list()
    subsamples.store[[1]] <- 1:n
    
    # Perform subsampling and clustering
    for (b in 1:B) {
        n_sub <- ceiling(cut_ratio * n)
        subsamples <- sample(1:n, n_sub, replace=FALSE)
        subsamples.store[[b+1]] <- subsamples
        df.sub <- df[subsamples,]
        
        km.sub <- stats::kmeans(df.sub, centers=nk, nstart=r)
        clst.mat[(b+1), subsamples] <- km.sub$cluster
    }
    
    # Calculate stability measures
    stab.mat <- matrix(NA, nrow=B+1, ncol=n)
    clst.ref <- clst.mat[1,]
    for (i in 1:B+1) {
        stab.mat[i,] <- agreement_nk(clst.ref[subsamples.store[[i]]],
                                   clst.mat[i,subsamples.store[[i]]],
                                   nk)
    }
    
    results$membership <- clst.ref
    results$obs_wise <- colMeans(stab.mat, na.rm=TRUE)
    results$cluster.matrix <- clst.mat
    results$runtime <- runtime
    return(results)
}

#' Subsampling-based Hierarchical Clustering
#'
#' @param df data.frame or matrix of the dataset
#' @param nk number of clusters
#' @param B number of bootstrap samples
#' @param cut_ratio ratio for subsampling (default: 0.5)
#' @param hc.method hierarchical clustering method (default: "ward.D")
#' @param dist_method distance method (default: "euclidean")
#' @return List containing clustering results and stability measures
#' @importFrom stats hclust cutree dist
#' @keywords internal
scheme_sub_hc <- function(df, nk, B=10, cut_ratio=0.5, hc.method="ward.D", dist_method="euclidean") {
    results <- list()
    n <- nrow(df)
    clst.mat <- matrix(NA, nrow=B+1, ncol=n)
    
    # Initial hierarchical clustering with timing
    start_time <- Sys.time()
    hc <- stats::hclust(stats::dist(df, method=dist_method), method=hc.method)
    clst <- stats::cutree(hc, k=nk)
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
    
    clst.mat[1,] <- clst
    
    # Store subsamples
    subsamples.store <- list()
    subsamples.store[[1]] <- 1:n
    
    # Perform subsampling and clustering
    for (b in 1:B) {
        n_sub <- ceiling(cut_ratio * n)
        subsamples <- sample(1:n, n_sub, replace=FALSE)
        subsamples.store[[b+1]] <- subsamples
        df.sub <- df[subsamples,]
        
        hc.sub <- stats::hclust(stats::dist(df.sub, method=dist_method), method=hc.method)
        clst.mat[(b+1), subsamples] <- stats::cutree(hc.sub, k=nk)
    }
    
    # Calculate stability measures
    stab.mat <- matrix(NA, nrow=B+1, ncol=n)
    clst.ref <- clst.mat[1,]
    for (i in 1:B+1) {
        stab.mat[i,] <- agreement_nk(clst.ref[subsamples.store[[i]]],
                                   clst.mat[i,subsamples.store[[i]]],
                                   nk)
    }
    
    results$membership <- clst.ref
    results$obs_wise <- colMeans(stab.mat, na.rm=TRUE)
    results$cluster.matrix <- clst.mat
    results$runtime <- runtime
    return(results)
}

#' Subsampling-based Spectral Clustering
#'
#' @param df data.frame or matrix of the dataset
#' @param nk number of clusters
#' @param B number of bootstrap samples
#' @param cut_ratio ratio for subsampling (default: 0.5)
#' @return List containing clustering results and stability measures
#' @importFrom kernlab specc
#' @keywords internal
scheme_sub_spectral <- function(df, nk, B=10, cut_ratio=0.5) {
    results <- list()
    n <- nrow(df)
    clst.mat <- matrix(NA, nrow=B+1, ncol=n)
    
    # Helper function for spectral clustering with improved error handling
    do_spectral <- function(data) {
        # Scale and clean the data
        scaled_data <- scale(as.matrix(data))
        
        # Handle any remaining NAs or Infs
        if (any(is.na(scaled_data)) || any(is.infinite(scaled_data))) {
            # Replace NAs with column means
            for (j in 1:ncol(scaled_data)) {
                na_idx <- is.na(scaled_data[,j]) | is.infinite(scaled_data[,j])
                if (any(na_idx)) {
                    scaled_data[na_idx,j] <- mean(scaled_data[!na_idx,j])
                }
            }
        }
        
        # Use spectral clustering with warning handling
        result <- NULL
        
        withCallingHandlers({
            # Try with increased iterations first
            sc <- kernlab::specc(scaled_data, centers=nk, iterations=1000)
            result <- sc@.Data
        }, warning = function(w) {
            # Handle convergence warnings
            if (grepl("did not converge", w$message)) {
                # Remove the warning
                invokeRestart("muffleWarning")
                
                # Try with a different kernel and parameters
                tryCatch({
                    sc <<- kernlab::specc(scaled_data, centers=nk, kernel="rbfdot", 
                                        kpar=list(sigma=0.2), iterations=2000)
                    result <<- sc@.Data
                }, error = function(e) {
                    # If that fails too, try with yet another approach
                    sc <<- kernlab::specc(scaled_data, centers=nk, kernel="polydot", 
                                        kpar=list(degree=2, scale=1, offset=1))
                    result <<- sc@.Data
                })
            }
        })
        
        return(result)
    }
    
    # Initial spectral clustering with timing and error handling
    start_time <- Sys.time()
    
    # Set up warning handler to capture warnings
    withCallingHandlers({
        # Try with increased iterations first
        sc <- kernlab::specc(as.matrix(df), centers=nk, iterations=1000)
        clst <- sc@.Data
    }, warning = function(w) {
        # Handle convergence warnings
        if (grepl("did not converge", w$message)) {
            message("Warning in spectral clustering: ", w$message, ". Trying with different parameters...")
            # Remove the warning
            invokeRestart("muffleWarning")
            
            # Try with a different kernel and parameters
            tryCatch({
                sc <<- kernlab::specc(as.matrix(df), centers=nk, kernel="rbfdot", 
                                    kpar=list(sigma=0.2), iterations=2000)
                clst <<- sc@.Data
            }, error = function(e) {
                # If that fails too, try with yet another approach
                message("Second attempt failed: ", e$message, ". Trying one more approach...")
                sc <<- kernlab::specc(as.matrix(df), centers=nk, kernel="polydot", 
                                    kpar=list(degree=2, scale=1, offset=1))
                clst <<- sc@.Data
            })
        }
    })
    
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
    
    clst.mat[1,] <- clst
    
    # Store subsamples
    subsamples.store <- list()
    subsamples.store[[1]] <- 1:n
    
    # Perform subsampling and clustering with retries
    for (b in 1:B) {
        success <- FALSE
        max_attempts <- 10  # Maximum number of attempts per subsample
        attempt <- 1
        
        while (!success && attempt <= max_attempts) {
            tryCatch({
                n_sub <- ceiling(cut_ratio * n)
                subsamples <- sample(1:n, n_sub, replace=FALSE)
                df.sub <- df[subsamples,]
                
                clst.mat[(b+1), subsamples] <- do_spectral(df.sub)
                subsamples.store[[b+1]] <- subsamples
                success <- TRUE
            }, error = function(e) {
                if (attempt == max_attempts) {
                    stop("Failed to perform spectral clustering on subsample dataset after ", max_attempts, " attempts")
                }
                attempt <- attempt + 1
            })
        }
        
        if (!success) {
            stop("Failed to complete spectral clustering for subsample dataset ", b)
        }
    }
    
    # Calculate stability measures
    stab.mat <- matrix(NA, nrow=B+1, ncol=n)
    clst.ref <- clst.mat[1,]
    for (i in 1:B+1) {
        stab.mat[i,] <- agreement_nk(clst.ref[subsamples.store[[i]]],
                                   clst.mat[i,subsamples.store[[i]]],
                                   nk)
    }
    
    results$membership <- clst.ref
    results$obs_wise <- colMeans(stab.mat, na.rm=TRUE)
    results$cluster.matrix <- clst.mat
    results$runtime <- runtime
    return(results)
}
