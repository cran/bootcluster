myoverallscheme <- function(data.input, thresh, B=20,
                            large.size = large.size,
                            cor.method = cor.method){
  results <- list()
  # define a matrix of cluster memberships
  boost.result <- boost.community(
    thresh = thresh,
    Boot = B,
    large.size = large.size,
    data.input = data.input,
    cor.method = cor.method
  )
  
  
  clst.mat<-do.call(rbind,boost.result$c.mat)
  
  
  B1 <- B+1
  agree.mat <- matrix(NA, nrow=B1, ncol=B1)
  diag(agree.mat) <- 1
  
  for (i in 1:(B1-1))
  {
    for (j in (i+1):B1)
    {
      agree.mat[i,j] <- mean(agreement(clst.mat[i,], clst.mat[j,]))
      agree.mat[j,i] <- agree.mat[i,j]
    }
  }
  
  mean.agr<- rowMeans(agree.mat)
  ref <- which.max(mean.agr)
  
  # stab.mat <- stab.mat[-ref, , drop = FALSE]
  
  
  results$cluster.matrix <- clst.mat
  results$agree.matrix <- agree.mat
  results$ref.cluster <- ref
  results$mscore<-boost.result$m.mat
  results$data<-boost.result$data_keep
  results$graph<-boost.result$graph_keep
  results$adjacency<-boost.result$adjacency_keep
  return(results)
}

