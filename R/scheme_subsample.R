scheme_subsample <- function(df, k, B=500, r=5, cut_ratio = 0.5){
  results <- list()
  n <- nrow(df)
  
  # define a matrix of cluster memberships
  clst.ref <- matrix(NA, nrow=n, ncol=1)
  
  clst.mat.OOB.m.OOB <- matrix(rep(NA, B), nrow = n, ncol = B, byrow = TRUE)
  # new object to store OOB inforation
  clst.mat.OOB.m.IN<- matrix(rep(NA, B), nrow = n, ncol = B, byrow = TRUE) 
  # new object to store OOB -> IN inforation
  
  
  km <- kmeansCBI(data=df, k=k, scaling=FALSE, runs=r)
  centers <- km$result$centers
  clst <- km$result$cluster
  
  clst.ref <- clst
  
  #stab.matrix <- matrix(NA, nrow=n, ncol=B)
  
  resample_store <- c(1:n) # original data indices
  
  
  inbag_points.list <- list()
  oob_points.list <- list()
  
  centers_in.list <- list()
  centers_out.list <- list()
  
  oob_stab_store <- matrix(rep(NA,n*B), nrow = n, ncol = B, byrow = TRUE)
  
  oob_clust_stability_store <- matrix(NA, nrow = k, ncol = B)
  
  for (b in 1:B){
    
    # Resample the dataset with the cut_ratio
    some_size <- ceiling(n*cut_ratio)
    some_size_resample <- sample(1:n, some_size, replace = FALSE)
    resample <- some_size_resample
    
    inbag_points <- resample
    oob_points <- setdiff(1:n, inbag_points) # resample the data with replacement
    
    df.star <- df[inbag_points, ] # the bootstrap indicis
    df.out <- df[oob_points, ] # the out of bag indicis
    
    # Store the inbag and oob points into a list
    inbag_points.list[[b]] <- inbag_points
    oob_points.list[[b]] <- oob_points
    
    # calculate the bootstrap clusterings    
    {
      C_in <- kmeansCBI(data=df.star, k=k, scaling=FALSE, runs=r) 
      #clustering the in bag
      
      C_out <- kmeansCBI(data=df.out, k=k, scaling=FALSE, runs=r) 
      #clustering the out of bag
      
      # Calculating centers
      centers_in <- C_in$result$centers 
      # grab the centers for the bootstrapped clustering
      
      centers_out <- C_out$result$centers 
      # grab the centers for the bootstrapped clustering
      
      # storing the centers into the i_j_k_l.list
      centers_in.list[[b]] <- centers_in
      centers_out.list[[b]] <- centers_out
    }
    
    # note this is pushing through the OOB data to the OOB Centers
    # cluster membership of the OOB data mapped to OOB centers
    mapper_out <- mapping.Euclidean(centers_out, df.out, 1:k) 
    
    # storing cluster membership for each OOB sample mapped to the OOB centers
    clst.mat.OOB.m.OOB[oob_points, b] <- mapper_out 
    
    # note this is pushing through the OOB data to the Bootstrap Centers
    # cluster membership of the OOB data mapped to Bootstrapped centers
    mapper_out_in <- mapping.Euclidean(centers_in, df.out, 1:k) 
    # storing cluster membership for each OOB sample mapped to the Bootstrapped centers
    clst.mat.OOB.m.IN[oob_points, b] <- mapper_out_in 
    
    # store the inbag points
    resample_store <- cbind(resample_store, resample)
    
    ##################
    # OOB -bag stability
    # grab the OOB bag indices for the cluster sets
    cluster_set1_oob <- clst.mat.OOB.m.OOB[,b] #OOB to OOB mappings
    cluster_set2_oob <- clst.mat.OOB.m.IN[,b] #OOB to IN mappings
    
    oob_stab <- agreement_nk(cluster_set1_oob[oob_points], cluster_set2_oob[oob_points], k)
    oob_stab_store[oob_points, b] <- oob_stab
    
  }
  
  
  
  
  results$membership        <- clst.ref
  results$obs_wise          <- apply(oob_stab_store, 1, 'mean', na.rm = TRUE)
  results$oob_stab_matrix   <- oob_stab_store
  results$oob_points_list   <- oob_points.list
  
  
  return(results)
  
}

