
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
  
  n1 <- length(clst1)
  n1_occur <- data.frame(table(clst1))
  nk1<- as.numeric(n1_occur[n1_occur$Freq >1,]$clst1)
  
  #nk2 <- length(unique(clst2))
  n2 <- length(clst2)
  n2_occur <- data.frame(table(clst2))
  nk2<- as.numeric(n2_occur[n2_occur$Freq >1,]$clst2)
  if(n1!=n2)warning('sample size is not equal')
  n<-n1
  if(length(nk1)==0|length(nk2)==0){
    stab.vec<- rep(0,n)
  }else{
    cluster.sets.1 <- list()
    cluster.sets.2 <- list()
    for (i in nk1)
    {
      cluster.sets.1[[i]] <- which(clst1==i)
      
    }
    
    for (i in nk2)
    {
      cluster.sets.2[[i]] <- which(clst2==i)
      
    }
    
    
    jaccard.matrix <- matrix(0, nrow=length(unique(clst1)), ncol=length(unique(clst2)))
    for (i in 1:length(nk1)){
      for (j in 1:length(nk2)){
        jaccard.matrix[i,j] <- jaccard(cluster.sets.1[[i]], cluster.sets.2[[j]])
      }
    }
    
    stab.vec <- c()
    for (i in 1:n){
      memb <- which(levels(n1_occur$clst1)==clst1[i])
      memb.star <- which(levels(n2_occur$clst2)==clst2[i])
      
      stab.vec[i] <- jaccard.matrix[memb, memb.star]
    }
    
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

community_coexist<-function(net.igraph,dats_met,kept.nodes){
  cmt <- fastgreedy.community(net.igraph) 
  modularity_score<-mean(cmt$modularity)
  if(length(kept.nodes)==dim(dats_met)[2]) cluster.output <- data.frame(node = kept.nodes, 
                                                                        cluster = cmt$membership)
  if(length(kept.nodes)!=dim(dats_met)[2]){
    cluster.temp <- data.frame(node = kept.nodes, cluster = cmt$membership)
    cluster.orig<-data.frame(node=colnames(dats_met))
    cluster.output <- suppressWarnings(left_join(cluster.orig,cluster.temp,by='node')) 
    index<-which(is.na(cluster.output$cluster))
    max.cluster<- max(cluster.output$cluster,na.rm = TRUE)
    
    for( i in 1:length(index)){
      cluster.output[index[i],]$cluster<- max.cluster+i
    }
    
    
  }
  
  result<-list(
    cluster.output=cluster.output,
    modularity_score=modularity_score)
  return(result)
}

community_membership<-function(cor_mat,thresh,large.size,dats_met){
  A <- cor_mat ## starting point
  A[abs(cor_mat)>thresh] <- 1
  A[abs(cor_mat)<=thresh] <- 0 # kill everything under a thresh
  A[is.na(cor_mat)] <- 0
  diag(A) <- 0
  # make it into a graph
  igraph.0 <- graph_from_adjacency_matrix(A, mode = c("undirected"))
  
  node.0 <- V(igraph.0)$name
  cls <- clusters(igraph.0)
  large.clusters <- which(cls$csize >= large.size)
  orig.membership<-ifelse( cls$membership %in% large.clusters,large.clusters,-1)
  # only retain nodes in clusters size >= 5
  kept.ind <- which(cls$membership %in% large.clusters)
  kept.nodes <- node.0[kept.ind]
  # sub.path <- metab_annotate_2017$SUPER_PATHWAY[kept.ind]
  orig.memb <- data.frame(node=node.0,cluster=orig.membership)
  net.igraph <- induced.subgraph(igraph.0, kept.nodes)
  data_keep<- dats_met[,colnames(dats_met)%in%kept.nodes,drop=FALSE]
  
  # community detection
  if(length(kept.nodes)!=0){
    r1<-community_coexist(net.igraph=net.igraph,dats_met=data_keep,kept.nodes=kept.nodes)
    
    cluster.output<-r1$cluster.output
    modularity_score<-r1$modularity_score
    
    
  }else{
    
    modularity_score=0
    
    cluster.output<-0
    
  }
  result<-list(cluster.output=cluster.output,
               modularity_score=modularity_score,
               data_keep=data_keep,
               igraphstore=net.igraph,
               adjacencymatrix= A)
  gc()
  
  return(result)
  
}

community_membership_boot<-function(cor_mat,thresh,dats_met){
  A <- cor_mat ## starting point
  A[abs(cor_mat)>thresh] <- 1
  A[abs(cor_mat)<=thresh] <- 0 # kill everything under a thresh
  A[is.na(cor_mat)] <- 0
  diag(A) <- 0
  # make it into a graph
  igraph.0 <- graph_from_adjacency_matrix(A, mode = c("undirected"))
  kept.nodes<-V(igraph.0)$name
  
  #@@@ Generate a ramdom graph for substraction
  coexist_seq_random<-NULL
  # community detection
  if(length(kept.nodes)!=0){
    r1<-community_coexist(net.igraph=igraph.0,dats_met=dats_met,kept.nodes=kept.nodes)
    cluster.output<-r1$cluster.output
    modularity_score<-r1$modularity_score
    
  }else{
    
    modularity_score=0
    subtract.mat<-0
    cluster.output<-0
    
  }
  
  result<-list(cluster.output=cluster.output,
               modularity_score=modularity_score,
               igraphstore=igraph.0,
               adjacencymatrix= A)
  gc()
  
  return(result)
  
}

boost.community<-function(thresh,Boot,data.input,large.size=5,
                          cor.method='spearman'){
  
  dats_met <- data.input
  colnames(dats_met)<-paste('M',c(1:dim(dats_met)[2]),sep = "")
  aux <- apply(dats_met, 2, var)
  if(length(which(aux == 0))!=0){
    dats_met_zeroremove <- dats_met[,-which(aux == 0)]}else{
      dats_met_zeroremove<-dats_met
    }
  
  cor.mat <- cor(dats_met_zeroremove, method = cor.method)
  ori_result<-community_membership(cor_mat = cor.mat,
                                   thresh = thresh,
                                   large.size = large.size,
                                   dats_met=dats_met_zeroremove)
  
  cluster.output_ori<-ori_result$cluster.output$cluster
  modularity_score_ori<-ori_result$modularity_score
  data_keep<-list()
  data_keep[[1]]<-ori_result$data_keep
  data_orig<- data_keep[[1]]
  graph_keep<-list()
  graph_keep[[1]]<-ori_result$igraphstore
  adjacency_keep<-list()
  adjacency_keep[[1]]<-ori_result$adjacencymatrix
  #################### Bootstrap and  Cluster ######################
  n = dim(dats_met)[1]
  sample.seq <- seq(1, n, 1)
  boot.index = matrix(sample(sample.seq, size = Boot * n, replace = TRUE),
                      Boot, n)
  cluster.track<-NULL
  cluster.track[[1]]<-cluster.output_ori
  modularity_score_boot<-NULL
  for(i in 1:Boot){
    dats_met_sample<-data_orig[boot.index[i,],,drop=FALSE]
    aux.sample <- apply(dats_met_sample, 2, var)
    zero.cov<-which(aux.sample == 0)
    if(length(zero.cov)==0)A_sample<-cor(dats_met_sample, method = cor.method)
    if(length(zero.cov)!=0)A_sample<-cor(dats_met_sample[,-zero.cov], method = cor.method)
    
    
    Boot.result<-community_membership_boot(cor_mat = A_sample,
                                           thresh = thresh,
                                           dats_met=dats_met_sample)
    
    
    modularity_score_boot[i]<-Boot.result$modularity_score
    cluster.track[[i+1]]<-Boot.result$cluster.output$cluster
    data_keep[[i+1]]<-dats_met_sample
    graph_keep[[i+1]]<-Boot.result$igraphstore
    adjacency_keep[[i+1]]<-Boot.result$adjacencymatrix
  }
  
  c.mat<-cluster.track
  m.mat<-c(modularity_score_boot,modularity_score_ori)
  result<-list(c.mat=c.mat,
               m.mat=m.mat,
               data_keep=data_keep,
               graph_keep=graph_keep,
               adjacency_keep=adjacency_keep)
  return(result)
}



scheme2.exp <- function(graph.input,
                        data.input,
                        clst.mat, thresh, 
                        PermuNo=10,B=20){
  cluster.output <-c()
  modularity_score<-c()
  degree.seq <- degree(graph.input) 
  kept.nodes<-V(graph.input)$name
  for(i in 1:PermuNo){
    
    random.graph <- sample_degseq(out.deg = degree.seq, method = "simple.no.multiple") 
    r1<-community_coexist(net.igraph=random.graph,dats_met=data.input,kept.nodes=kept.nodes)
    cluster.output<-rbind(cluster.output,r1$cluster.output$cluster)
    modularity_score[i]<-r1$modularity_score
  }
  
  
  permu.mat<-cluster.output
  
  
  B1 <- PermuNo
  agree.mat <- matrix(NA, nrow=B1, ncol=B1)
  diag(agree.mat) <- 1
  
  for (i in 1:(B1-1))
  {
    for (j in (i+1):B1)
    {
      agree.mat[i,j] <- mean(agreement(permu.mat[i,], permu.mat[j,]))
      agree.mat[j,i] <- agree.mat[i,j]
    }
  }
  
  results<-list()
  results$cluster.matrix <- permu.mat
  results$agree.matrix <- agree.mat
  results$mscore<-modularity_score
  
  
  return(results)
}

scheme2.module <- function(data.input, thresh, B,
                        large.size ,
                        cor.method){
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

grid_arrange_shared_legend <- function(plots) {
  
  g <- ggplotGrob(plots[[length(plots)]] + theme(legend.position = "bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$heights)
  plots[[length(plots)]] <- NULL
  
  p1<-arrangeGrob(grobs=lapply(plots, function(x)
    x + theme(legend.position = "none")),nrow = ceiling((length(plots)-1)/5))
  
  grid.arrange(
    p1,
    legend,
    nrow = 2,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
  
}
