#' Generate Reference Distribution
#'
#' @title Generate reference distribution for stability assessment
#'
#' @description Generates a reference distribution by sampling from uniform distributions
#' with ranges determined by the original data.
#'
#' @param df data.frame or matrix of the original dataset
#' @return A scaled matrix containing the reference distribution
#'
#' @importFrom stats runif
#'
#' @examples
#' \donttest{
#' data(iris)
#' df <- iris[,1:4]
#' ref <- ref_dist(df)
#' }
#'
#' @export
ref_dist <- function(df) {
    # Generate reference data from uniform distributions
    ref_dist_range <- apply(df, 2, range)
    ref_dist <- matrix(NA, nrow=nrow(df), ncol=ncol(df))
    
    for (i in 1:ncol(ref_dist)) {
        ref_dist[,i] <- runif(nrow(ref_dist),
                             min=ref_dist_range[1,i],
                             max=ref_dist_range[2,i])
    }
    
    # Use base::scale instead of importing scale
    df_ref <- base::scale(ref_dist)
    return(df_ref)
}

#' Generate Reference Distribution using PCA
#'
#' @title Generate PCA-based reference distribution
#'
#' @description Generates a reference distribution in PCA space by sampling from uniform
#' distributions with ranges determined by the PCA-transformed data.
#'
#' @param df data.frame or matrix of the original dataset
#' @return A scaled matrix containing the reference distribution in PCA space
#'
#' @importFrom stats prcomp runif
#'
#' @examples
#' \donttest{
#' data(iris)
#' df <- iris[,1:4]
#' ref <- ref_dist_pca(df)
#' }
#'
#' @export
ref_dist_pca <- function(df) {
    # Transform to PCA space
    pca <- prcomp(df, scale.=FALSE)
    pca_df <- pca$x
    ref_dist_range <- apply(pca_df, 2, range)
    ref_dist <- matrix(NA, nrow=nrow(pca_df), ncol=ncol(pca_df))
    
    for (i in 1:ncol(ref_dist)) {
        ref_dist[,i] <- runif(nrow(ref_dist),
                             min=ref_dist_range[1,i],
                             max=ref_dist_range[2,i])
    }
    
    # Use base::scale instead of importing scale
    df_ref <- base::scale(ref_dist)
    return(df_ref)
}

#' Generate Binary Reference Distribution
#'
#' @title Generate reference distribution for binary data
#'
#' @description Generates a reference distribution by randomly permuting each column
#' of the original binary dataset.
#'
#' @param df data.frame or matrix of the original binary dataset
#' @return A matrix containing the permuted binary reference distribution
#'
#' @examples
#' \donttest{
#' binary_data <- matrix(sample(0:1, 100, replace=TRUE), ncol=5)
#' ref <- ref_dist_bin(binary_data)
#' }
#'
#' @export
ref_dist_bin <- function(df) {
    df_ref <- data.frame(matrix(data=0, nrow=nrow(df), ncol=ncol(df)))
    n_rows <- nrow(df)
    n_cols <- ncol(df)
    
    for (col in 1:n_cols) {
        resample_indexes <- sample(1:n_rows, n_rows, replace=FALSE)
        df_ref[,col] <- df[resample_indexes,col]
    }
    
    return(df_ref)
}
