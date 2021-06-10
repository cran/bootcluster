#' @name network.stability
#'
#' @title  Estimate of detect module stability
#'
#' @details This function estimates the modules' stability through bootstrapping 
#' approach for the given threshold. 
#' The approach to stability estimation is to compare the module composition of the reference correlation 
#' graph to the various bootstrapped correlation graphs, and to assess the stability at the 
#' (1) node-level, (2) module-level, and (3) overall. 
#'
#' @param data.input a \code{data.frame} of the data set where the rows are observations and columns are covariates 
#' @param threshold a \code{numeric} number of threshold for correlation matrix
#' @param B number of bootstrap re-samplings
#' @param cor.method the correlation method applied to the data set,three method are available: \code{"pearson", "kendall", "spearman"}.
#' @param large.size the smallest set of modules, the \code{large.size=0} is recommended to use right now. 
#' @param PermuNo number of random graphs for null 
#' @param scheme_2 \code{logical} \code{TRUE} if scheme 2 is used, \code{FASLE} if scheme 1 is used. 
#' Right now, only \code{FASLE} is recommended.
#' 
#'
#' @return
#' \describe{
#' \item{\code{stabilityresult}}{a \code{list} of result for nodes-wise stability}
#' \item{\code{modularityresult}}{\code{list} of modularity information with the given threshold}
#' \item{\code{jaccardresult}}{\code{list} estimated unconditional observed stability and 
#'      the estimates of expected stability under the null}
#' \item{\code{originalinformation}}{\code{list} information for original data,
#'       igraph object and adjacency matrix constructed with the given threshold}
#' }
#'
#' @author Mingmei Tian
#'
#' @references A framework for stability-based module detection in correlation graphs.
#' Mingmei Tian,Rachael Hageman Blair,Lina Mu, Matthew Bonner, Richard Browne and Han Yu.
#'
#' @importFrom igraph graph_from_adjacency_matrix fastgreedy.community induced.subgraph  V clusters degree sample_degseq 
#' @importFrom dplyr left_join
#' @importFrom compiler cmpfun
#' @importFrom doParallel registerDoParallel 
#' @importFrom parallel clusterExport makeCluster
#' @importFrom foreach foreach
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' data(wine)
#' x0 <- wine[1:50,]
#' 
#' mytest<-network.stability(data.input=x0,threshold=0.7, B=20, 
#' cor.method='pearson',large.size=0,
#' PermuNo = 10,
#' scheme_2 = FALSE)
#' }
#' @export

network.stability<-function(data.input,threshold, B=20, 
                            cor.method,large.size,
                            PermuNo ,
                            scheme_2 = FALSE){
  
  
  s2 <- scheme2.module(data.input=data.input, thresh=threshold,
                B=B,large.size=large.size,
                cor.method=cor.method)
  c.mat <- s2$cluster.matrix
  mscore<-s2$mscore
  
  ref <- 1
  if (scheme_2) {
    ref <- s2$ref.cluster
  }
  
  ref.cl <- c.mat[ref,]
  c.mat.1 <- c.mat[-ref,]
  min.agr <- c()
  for (i in 1:B)
  {
    min.agr[i] <- min.agreement(ref.cl, agreement(ref.cl, c.mat.1[i,]))
  }
  crit=mean(min.agr)
  
  
  datainputforexp<- s2$data[[ref]]
  graphinputforexp<-s2$graph[[ref]]
  
  
  s2.exp <- scheme2.exp(graph.input=graphinputforexp,
                        data.input=datainputforexp,
                        clst.mat=c.mat, 
                        thresh=threshold, 
                        PermuNo=PermuNo,B=B)
  p.mat <- s2.exp$cluster.matrix
  mscore.exp<-s2.exp$mscore
  
  
  clst.ref <- c.mat[ref, , drop = TRUE]
  
  stab.mat <- matrix(NA, nrow = (B+1), ncol =  dim(c.mat)[2])
  
  for (i in 1:(B+1)) {
    stab.mat[i,] <- agreement(clst.ref, c.mat[i,])
  }
  
  stabilityresult<-list()
  stabilityresult$ref<-ref
  stabilityresult$membership <- clst.ref
  stabilityresult$obs_wise <- colMeans(stab.mat)
  stabilityresult$overall <- mean(stabilityresult$obs_wise)
  
  modularityresult<-list()
  modularityresult$mscore=mscore
  modularityresult$mscore.exp=mscore.exp
  
  jaccardresult<-list()
  jaccardresult$min=crit
  jaccardresult$obsvalue=rowMeans(s2$agree.matrix)
  jaccardresult$expvalue=mean(s2.exp$agree.matrix,na.rm = TRUE)
  
  originalinformation<-list()
  originalinformation$data=s2$data[[1]]
  originalinformation$graph=s2$graph[[1]]
  originalinformation$adjacency=s2$adjacency[[1]]
  if(scheme_2==TRUE){
    originalinformation$data2=datainputforexp
    originalinformation$graph2=graphinputforexp
    originalinformation$adjacency2=s2$data[[ref]]
  }
  
  
  
  result<-list(jaccardresult=jaccardresult,
               stabilityresult=stabilityresult,
               modularityresult=modularityresult,
               originalinformation=originalinformation)
  
  return(result)
  
}

