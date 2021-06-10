
#' Estimate clustering stability of k-means
#'
#' Estimate of k-means bootstrapping stability
#'
#' @details This function estimates the clustering stability through bootstrapping approach. 
#' Two schemes are provided. Scheme 1 uses the clustering of the original data as the reference 
#' for stability calculations. Scheme 2 searches acrossthe clustering samples that gives the 
#' most stable clustering.
#'
#' @param x a \code{data.frame} of the data set
#' @param k a \code{integer} number of clusters
#' @param B number of bootstrap re-samplings
#' @param r number of runs of k-means
#' @param scheme_2 \code{logical} \code{TRUE} if scheme 2 is used, \code{FASLE} if scheme 1 is used
#'
#' @return
#' \describe{
#' \item{\code{membership}}{a \code{vector} of membership for each observation from the reference clustering}
#' \item{\code{obs_wise}}{\code{vector} of estimated observation-wise stability}
#' \item{\code{overall}}{\code{numeric} estimated overall stability}
#' }
#'
#' @author Han Yu
#'
#' @references Bootstrapping estimates of stability for clusters, observations and model selection.
#' Han Yu, Brian Chapman, Arianna DiFlorio, Ellen Eischen, David Gotz, Matthews Jacob and Rachael Hageman Blair.
#'
#' @import cluster mclust fpc plyr
#' @importFrom flexclust dist2
#'
#' @examples
#'  \donttest{
#' set.seed(1)
#' data(wine)
#' x0 <- wine[,2:14]
#' x <- scale(x0)
#' stability(x, k = 3, B=20, r=5, scheme_2 = TRUE)
#' }
#' @export


stability <- function(x, k, B=20, r=5, scheme_2 = TRUE) {

  result <- list()

  output <- scheme2(df = x, nk = k, B = B, r = r)

  clst.mat <- output$cluster.matrix
  agree.mat <- output$agree.matrix

  ref <- 1

  if (scheme_2) {
    ref <- output$ref.cluster
  }

  result$membership <- clst.mat[ref, , drop = FALSE]
  result$obs_wise <- output$obs_wise
  result$overall <- mean(result$obs_wise)

  return(result)

}
