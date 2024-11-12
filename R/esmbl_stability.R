#' @name esmbl.stability
#'
#' @title  Estimate the stability of a clustering based on non-parametric bootstrap 
#' out-of-bag scheme, with option for subsampling scheme
#'
#' @details This function estimates the stability through out-of-bag observations 
#' It estimate the stability at the 
#' (1) observation level, (2) cluster level, and (3) overall. 
#'
#' @param x \code{data.frame} of the data set where rows are observations and columns are features  
#' @param k number of clusters for which to estimate the stability   
#' @param scheme clustering method to use ("kmeans", "hc", or "spectral")
#' @param B number of bootstrap re-samples
#' @param hc.method hierarchical clustering method (default: "ward.D")
#' @param cut_ratio ratio for subsampling (default: 0.5)
#' @param dist_method distance method for spectral clustering (default: "euclidean")
#' 
#' @return
#' \describe{
#' \item{membership}{vector of membership for each observation from the reference clustering}
#' \item{obs_wise}{vector of estimated observation-wise stability}
#' \item{clust_wise}{vector of estimated cluster-wise stability}
#' \item{overall}{numeric estimated overall stability}
#' \item{Smin}{numeric estimated Smin through out-of-bag scheme}
#' }
#'
#' @author Tianmou Liu
#' 
#' @import cluster mclust fpc plyr 
#' @importFrom flexclust dist2
#' @importFrom kernlab specc
#' @importFrom stats hclust cutree dist kmeans
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' data(iris)
#' df <- iris[,1:4]
#' result <- esmbl.stability(df, k=3, scheme="kmeans")
#' }
#' 
#' @export
esmbl.stability <- function(x, k, scheme="kmeans", B=100, hc.method="ward.D", cut_ratio=0.5, dist_method="euclidean") {
    result <- list()
    
    if (scheme == "kmeans") {
        output <- scheme_sub_km(df=x, nk=k, B=B, cut_ratio=cut_ratio)
    } else if (scheme == "hc") {
        output <- scheme_sub_hc(df=x, nk=k, B=B, cut_ratio=cut_ratio, hc.method=hc.method, dist_method=dist_method)
    } else if (scheme == "spectral") {
        output <- scheme_sub_spectral(df=x, nk=k, B=B, cut_ratio=cut_ratio)
    } else {
        stop("Invalid scheme. Use 'kmeans', 'hc', or 'spectral'")
    }
    
    cluster_membership <- output$membership
    observation_stability <- output$obs_wise
    
    # Calculate cluster-level stability
    cluster_level_stability <- matrix(NA, nrow=1, ncol=k)
    colnames(cluster_level_stability) <- paste("cluster", seq(1:k), sep="")
    
    for (i in 1:k) {
        cluster_i_observations <- which(cluster_membership == i)
        if (length(cluster_i_observations) > 0) {
            cluster_level_stability[i] <- mean(observation_stability[cluster_i_observations])
        } else {
            cluster_level_stability[i] <- 0
        }
    }
    
    overall_stability <- mean(observation_stability)
    Smin <- min(cluster_level_stability)
    
    result$membership <- cluster_membership
    result$obs_wise <- observation_stability
    result$clust_wise <- cluster_level_stability
    result$overall <- overall_stability
    result$Smin <- Smin
    
    return(result)
}
