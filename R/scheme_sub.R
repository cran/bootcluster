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
scheme_sub_km <- function(df, nk, B=20, cut_ratio=0.5, r=5) {
    results <- list()
    n <- nrow(df)
    clst.mat <- matrix(NA, nrow=B+1, ncol=n)
    
    # Initial k-means clustering
    km <- stats::kmeans(df, centers=nk, nstart=r)
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
scheme_sub_hc <- function(df, nk, B=20, cut_ratio=0.5, hc.method="ward.D", dist_method="euclidean") {
    results <- list()
    n <- nrow(df)
    clst.mat <- matrix(NA, nrow=B+1, ncol=n)
    
    # Initial hierarchical clustering
    hc <- stats::hclust(stats::dist(df, method=dist_method), method=hc.method)
    clst <- stats::cutree(hc, k=nk)
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
scheme_sub_spectral <- function(df, nk, B=20, cut_ratio=0.5) {
    results <- list()
    n <- nrow(df)
    clst.mat <- matrix(NA, nrow=B+1, ncol=n)
    
    # Initial spectral clustering
    sc <- kernlab::specc(as.matrix(df), centers=nk)
    clst <- sc@.Data
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
        
        sc.sub <- kernlab::specc(as.matrix(df.sub), centers=nk)
        clst.mat[(b+1), subsamples] <- sc.sub@.Data
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
    return(results)
}
