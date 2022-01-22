#' @name ob.stability
#'
#' @title  Estimate the stability of a clustering based on non-parametric bootstrap 
#' out-of-bag scheme, with option for subsampling scheme
#'
#' @details This function estimates the stability through out-of-bag observations 
#' It estimate the stability at the 
#' (1) observation level, (2) cluster level, and (3) overall. 
#'
#' @param x \code{data.frame} of the data set where the rows as observations and columns as dimensions of features  
#' @param k number of clusters for which to estimate the stability   
#' @param B number of bootstrap re-samples
#' @param r integer parameter in the kmeansCBI() funtion
#' @param subsample logical parameter to use the subsampling scheme option in the resampling process (instead of bootstrap)
#' @param cut_ratio numeric parameter between 0 and 1 for subsampling scheme training set ratio
#' 
#' 
#' @return
#' \describe{
#' \item{\code{membership}}{\code{vector} of membership for each observation from the reference clustering}
#' \item{\code{obs_wise}}{\code{vector} of estimated observation-wise stability}
#' \item{\code{clust_wise}}{\code{vector} of estimated cluster-wise stability}
#' \item{\code{overall}}{\code{numeric} estimated overall stability}
#' \item{\code{Smin}}{\code{numeric} estimated Smin through out-of-bag scheme}
#' }
#'
#'
#'
#' @author Tianmou Liu
#' 
#' @import cluster mclust fpc plyr
#' @importFrom flexclust dist2
#'
#' @references Bootstrapping estimates of stability for clusters, observations and model selection.
#' Han Yu, Brian Chapman, Arianna DiFlorio, Ellen Eischen, David Gotz, Matthews Jacob and Rachael Hageman Blair.
#' 
#' 
#' @examples
#' \donttest{
#' set.seed(123)
#' data(iris)
#' df <- data.frame(iris[,1:4])
#' # You can choose to scale df before clustering by 
#' # df <- scale(df)
#' ob.stability(df, k = 2, B=500, r=5)
#' }
#' 
#' @export

ob.stability <- function(x, k, B=500, r=5,
                         subsample = FALSE, cut_ratio = 0.5){
  
  result <- list()
  
  #Take the out-of-bag stability output and calculate the cluster-wise and overall stability
  if(subsample == FALSE){
    # If we don't use data split, then use the previous function scheme_ob.R
    output <- scheme_ob(df = x, k = k, B = B, r = r)
  }
  else{
    # If we use data split, then use the adpated version scheme_subsample.R
    output <- scheme_subsample(df = x, k = k, B = B, r = r, cut_ratio)
    
  }
  cluster_membership    <- output$membership
  observation_stability <- output$obs_wise

  # work on the cluster-level stability calculation
  cluster_level_stability <- matrix(NA,nrow = 1, ncol = k)
  colnames(cluster_level_stability) <- paste("cluster ",seq(1:k),sep = "")
  
  for (i in 1:k) {
    cluster_i_observations     <- which(cluster_membership==i)
    cluster_level_stability[i] <- mean(observation_stability[cluster_i_observations])
  }
  
  overall_stability <- mean(observation_stability)
  
  Smin = min(cluster_level_stability)
  
  result$membership <- cluster_membership
  result$obs_wise   <- observation_stability
  result$clust_wise <- cluster_level_stability
  result$overall    <- overall_stability
  result$Smin       <- Smin
  
  return(result)
  
}