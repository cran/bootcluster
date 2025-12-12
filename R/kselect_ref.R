#' Estimate number of clusters
#'
#' Estimate number of clusters by bootstrapping stability
#'
#' @details This function uses the out-of-bag scheme to estimate the number of clusters
#' in a dataset. The function calculate the Smin of the dataset and at the same time, generate
#' a reference dataset with the same range as the original dataset in each dimension and calculate
#' the Smin_ref. The differences between Smin and Smin_ref at each k,Smin_diff(k), is taken into consideration as well as the 
#' standard deviation of the differences. We choose the k to be the argmax of ( Smin_diff(k) - ( Smin_diff(k+1) + (Smin_diff(k+1)) ) ).
#' If Smin_diff(k) less than 0.1 for all k in k_range, we say k = 1
#' 
#'
#' @param df \code{data.frame} of the input dataset
#' @param k_range \code{integer} valued \code{vector} of the numbers of clusters k to be tested upon
#' @param n_ref number of reference distribution to be generated
#' @param B number of bootstrap re-samples
#' @param B_ref number of bootstrap resamples for the reference distributions
#' @param r number of runs of k-means
#'
#' @return
#' \describe{
#' \item{\code{profile}}{\code{vector} of ( Smin_diff(k) - ( Smin_diff(k+1) + se(Smin_diff(k+1)) ) ) measures for researchers's inspection}
#' \item{\code{k}}{estimated number of clusters}
#' }
#'
#' @author Tianmou Liu
#'
#' @references Bootstrapping estimates of stability for clusters, observations and model selection.
#' Han Yu, Brian Chapman, Arianna DiFlorio, Ellen Eischen, David Gotz, Matthews Jacob and Rachael Hageman Blair.
#'
#' @import cluster mclust fpc
#' @importFrom flexclust dist2
#' @importFrom stats runif sd
#' @examples
#' \donttest{
#' set.seed(1)
#' data(iris)
#' df <- data.frame(iris[,1:4])
#' df <- scale(df)
#' k.select_ref(df, k_range = 2:7, n_ref = 5, B=500, B_ref = 500, r=5)
#' }
#' @export

k.select_ref <- function(df, k_range = 2:7, n_ref = 5, B=100, B_ref = 50, r=5){
  
  result <- list()
  
  
  std <- function(x) sd(x)/sqrt(length(x))
  
  # This is the generated reference data
  ref_dist_range <- apply(df,2,range)
  ref_dist <- matrix(NA, nrow=nrow(df), ncol=ncol(df))
  for (i in 1:ncol(ref_dist)) {
    ref_dist[,i] <- runif(nrow(ref_dist),
                          min = ref_dist_range[1,i],
                          max = ref_dist_range[2,i])
  }

  df_ref <- ref_dist
  
  df_ref <- scale(df_ref)

  # Building the Smin matrixes
  # The reference distribution should check until k = k_range+1, so the matrix should have
  # one more column
  
  # k_length is for ease of indicating matrix positions
  k_range_Smin = min(k_range):(max(k_range)+1)
  
  k_length = length(k_range_Smin)
  
  # Declaring the matrixes
  Smin_matrix <- matrix(NA, nrow = 1, ncol = k_length)
  
  Smin_ref_matrix <- matrix(NA, nrow = n_ref, ncol = k_length)
  colnames(Smin_matrix) <- paste("k:",k_range_Smin,sep = "")
  colnames(Smin_ref_matrix) <- colnames(Smin_matrix)
  
  for (i in 1:k_length){
    # k_smin is the k that we check in the following calculations
    k_Smin = k_range_Smin[i]
    output_k = ob.stability(df,k_Smin,B = B)
    Smin_matrix[1,i] <- output_k$Smin
    for (ref in 1:n_ref) {
      
      print(paste(k_Smin,ref))
      output_k_ref = ob.stability(df_ref,k_Smin,B=B)
      Smin_ref_matrix[ref,i] = output_k_ref$Smin
    }
  }
  

  # Now lets try the plus se method
  
  # First check if the smin minus reference is all negative, if so then nk = 1
  
  Smin_diff = Smin_matrix - colMeans(Smin_ref_matrix)
  if(sum(Smin_diff > 0) == 0){
    nkoobse = 1
  } else{
    
    # Let's prepare the Smin matrix as same nrows as the reference simulations
    # to calculate difference mean and se
    Smin_matrix_collection = data.frame()
    for (i in 1:n_ref) {
      Smin_matrix_collection <- rbind(Smin_matrix_collection, Smin_matrix)
      
    }
    
    colnames(Smin_matrix_collection) = colnames(Smin_matrix)
    
    # Make the difference matrix
    smin_minus_ref = Smin_matrix_collection - Smin_ref_matrix
    
    # Take away the last column because there is no nk+2 to compare
    Smin_nk_for_compare = colMeans(smin_minus_ref)[-(k_length)]
    
    # Calculate the mean + se data, and take away the k = 2 column because we compare with plus 1
    Smin_nkplus1_plus_se = (colMeans(smin_minus_ref) + apply(smin_minus_ref,2,std))[-1]
    
    # This is finding where the Smin(k) - ( Smin(k+1) + se(k+1) ) is maxed and positive
    hunt_for_elbow = Smin_nk_for_compare - Smin_nkplus1_plus_se
    
    # because the huntfor elbow start from k = 2
    nkoobse = which.max(hunt_for_elbow)+1
    
    if(sum(hunt_for_elbow> .1) == 0){
      nkoobse =1
    }
    
    
  }
  
  result$profile <- hunt_for_elbow
  result$k <- nkoobse
  
  return(result)
  
}
