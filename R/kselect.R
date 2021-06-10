#' Estimate number of clusters
#'
#' Estimate number of clusters by bootstrapping stability
#'
#' @details This function estimates the number of clusters through a bootstrapping
#' approach, and a measure Smin, which is based on an observation-wise similarity
#' among clusterings. The number of clusters k is selected as the largest number of
#' clusters, for which the Smin is greater than a threshold. The threshold is often
#' selected between 0.8 ~ 0.9. Two schemes are provided. Scheme 1 uses the clustering
#' of the original data as the reference for stability calculations. Scheme 2 searches 
#' acrossthe clustering samples that gives the most stable clustering.
#'
#' @param x a \code{data.frame} of the data set
#' @param range a \code{vector} of \code{integer} values, of the possible numbers of clusters k
#' @param B number of bootstrap re-samplings
#' @param r number of runs of k-means
#' @param threshold the threshold for determining k
#' @param scheme_2 \code{logical} \code{TRUE} if scheme 2 is used, \code{FASLE} if scheme 1 is used
#'
#' @return
#' \describe{
#' \item{\code{profile}}{a \code{vector} of Smin measures for determining k}
#' \item{\code{k}}{\code{integer} estimated number of clusters}
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
#' \donttest{
#' set.seed(1)
#' data(wine)
#' x0 <- wine[,2:14]
#' x <- scale(x0)
#' k.select(x, range = 2:10, B=20, r=5, scheme_2 = TRUE)
#' }
#' @export

k.select <- function(x, range=2:7, B=20, r=5, threshold=0.8, scheme_2 = TRUE){
  df <- x
  result <- list()
  crit <- c() # matrix(NA, nrow=length(range), ncol=2) ##
  j <- 1
  for (k in range)
  {
    min.agr <- c()

    s2 <- scheme2(df=df, nk=range[j], B=B, r=r)
    c.mat <- s2$cluster.matrix

    ref <- 1
    if (scheme_2) {
      ref <- s2$ref.cluster
    }

    ref.scheme1 <- c.mat[1,] #
    c.mat.scheme1 <- c.mat[-1,] #
    min.agr.scheme1 <- c()

    ref.cl <- c.mat[ref,]
    c.mat.1 <- c.mat[-ref,]
    for (i in 1:B)
    {
      min.agr[i] <- min.agreement(ref.cl, agreement(ref.cl, c.mat.1[i,]))
      # min.agr.scheme1[i] <- #
       #  min.agreement(ref.scheme1, agreement(ref.scheme1, c.mat.scheme1[i,])) #
    }
    # crit[j,1] <- mean(min.agr.scheme1)
    crit[j] <- mean(min.agr)
    j <- j+1
  }
  names(crit) <- range

  result$profile <- crit
  result$k <- range[max(which(crit > threshold))]
  # colnames(crit) <- c("s1", "s2")
  return(result)
}

# s2 <- k.select(iris[,1:4], range=2:7, B=20, r=5, scheme_2 = TRUE)
