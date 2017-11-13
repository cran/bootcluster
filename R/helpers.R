
mapping.Euclidean <- function(center.ori, center.map, label.ori)
{
  dist.mat <- dist2(center.map, center.ori)
  mapping.pre <- apply(dist.mat, 1, which.min)
  mapping <- label.ori[mapping.pre]
  return(mapping)
}

jaccard <- function(set1, set2)
{
  # require("sets")
  jaccard <- length(intersect(set1, set2))/length(union(set1, set2))
  return(jaccard)
}

## obtain agreement between two clustering results
agreement <- function(clst1, clst2)
{
  nk <- length(unique(clst1))
  n <- length(clst1)

  cluster.sets.1 <- list()
  cluster.sets.2 <- list()
  for (i in 1:nk)
  {
    cluster.sets.1[[i]] <- which(clst1==i)
    cluster.sets.2[[i]] <- which(clst2==i)
  }

  jaccard.matrix <- matrix(NA, nrow=nk, ncol=nk)
  for (i in 1:nk){
    for (j in 1:nk){
      jaccard.matrix[i,j] <- jaccard(cluster.sets.1[[i]], cluster.sets.2[[j]])
    }
  }

  stab.vec <- c()

  for (i in 1:n){
    memb <- clst1[i]
    memb.star <- clst2[i]
    stab.vec[i] <- jaccard.matrix[memb, memb.star]
  }
  return(stab.vec)
}

min.agreement <- function(clst, agrmt){
  clst.list <- unique(clst)
  clst.sta <- c()
  n <- length(clst.list)
  for(i in 1:n) {
    clst.sta[i] <- mean(agrmt[clst==clst.list[i]])
  }
  min.agrmt <- min(clst.sta)
  return(min.agrmt)
}

