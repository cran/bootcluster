

scheme2 <- function(df, nk, B=20, r=5){
  results <- list()
  n <- nrow(df)
  # define a matrix of cluster memberships
  clst.mat <- matrix(NA, nrow=B+1, ncol=n)

  km <- kmeansCBI(data=df, k=nk, scaling=FALSE, runs=r)
  center <- km$result$centers
  clst <- km$result$cluster

  clst.mat[1,] <- clst

  stab.matrix <- matrix(NA, nrow=n, ncol=B)

  for (b in 1:B){
    resample <- sample(1:n, replace = TRUE)

    df.star <- df[resample,]
    km.star <- kmeansCBI(data=df.star, k=nk, scaling=FALSE, runs=r)
    center.star <- km.star$result$centers

    class.star <- mapping.Euclidean(center.star, df, 1:nk)
    clst.mat[(b+1),] <- class.star
  }

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
  clst.ref <- clst.mat[ref, , drop = TRUE]

  stab.mat <- matrix(NA, nrow = B1, ncol = n)

  for (i in 1:B1) {
    stab.mat[i, ] <- agreement(clst.ref, clst.mat[i,])
  }

  # stab.mat <- stab.mat[-ref, , drop = FALSE]

  results$membership <- clst.ref
  results$obs_wise <- colMeans(stab.mat)
  results$cluster.matrix <- clst.mat
  results$agree.matrix <- agree.mat
  results$ref.cluster <- ref
  return(results)
}







