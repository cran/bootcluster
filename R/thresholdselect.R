#' @name threshold.select
#'
#' @title Estimate of the overall Jaccard stability
#'
#' @details \code{threshold.select} is used to estimate of the overall Jaccard stability from 
#' a sequence of given threshold candidates, \code{threshold.seq}.
#' 
#' 
#' @param data.input a \code{data.frame} of the data set where the rows are observations and columns are covariates 
#' @param threshold.seq a \code{numeric} sequence of candidate threshold
#' @param B number of bootstrap re-samplings
#' @param cor.method the correlation method applied to the data set,three method are available: \code{"pearson", "kendall", "spearman"}.
#' @param large.size the smallest set of modules, the \code{large.size=0} is recommended to use right now. 
#' @param PermuNo number of random graphs for the estimation of expected stability
#' @param no_cores a \code{interger} number of CPU cores on the current host.
#'
#' @return
#' \describe{
#' \item{\code{stabilityresult}}{a \code{list} of result for nodes-wise stability}
#' \item{\code{modularityresult}}{a \code{list} of modularity information with each candidate threshold}
#' \item{\code{jaccardresult}}{a \code{list} estimated unconditional observed stability and 
#'      the estimates of expected stability under the nul}
#' \item{\code{originalinformation}}{a \code{list} information for original data,
#'       igraph object and adjacency matrix constructed with each candidate threshold}
#' \item{\code{threshold.seq}}{a \code{list} of candicate threshold given to the function}
#' }
#'
#' @author Mingmei Tian
#'
#' @references A framework for stability-based module detection in correlation graphs.
#' Mingmei Tian,Rachael Hageman Blair,Lina Mu, Matthew Bonner, Richard Browne and Han Yu.
#' @importFrom igraph graph_from_adjacency_matrix fastgreedy.community induced.subgraph  V clusters degree sample_degseq 
#' @importFrom dplyr left_join
#' @importFrom compiler cmpfun
#' @importFrom doParallel registerDoParallel 
#' @importFrom parallel clusterExport makeCluster stopCluster 
#' @importFrom foreach foreach %dopar%
#' @importFrom stats cor var
#' 
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' data(wine)
#' x0 <- wine[1:50,]
#' 
#' mytest<-threshold.select(data.input=x0,threshold.seq=seq(0.5,0.8,by=0.05), B=20, 
#' cor.method='pearson',large.size=0,
#' PermuNo = 10,
#' scheme_2 = FALSE)
#' }
#' @export


utils::globalVariables(c("k"))

threshold.select <- function(data.input,threshold.seq, B=20, 
                             cor.method,large.size,
                             no_cores=8,PermuNo = 10,
                             scheme_2 = FALSE){
  
  myFuncCmp <- cmpfun(network.stability)
  mcoptions <- list(preschedule = FALSE)
  # Initiate cluster
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  clusterExport(cl,'jaccard')
  clusterExport(cl,'min.agreement')
  clusterExport(cl,'community_coexist')
  clusterExport(cl,'community_membership')
  clusterExport(cl,'community_membership_boot')
  clusterExport(cl,'boost.community')
  clusterExport(cl,'scheme2.module')
  clusterExport(cl,'agreement')
  clusterExport(cl,'scheme2.exp')
  clusterExport(cl,'parelle.function')
  
  result<-foreach(k = threshold.seq ,
                  .packages=c("igraph", "reshape2", "plyr", "qvalue",
                              "wfg","geneplotter","dplyr","jaccard"),
                  .options.multicore = mcoptions
  ) %dopar% {
    
    myFuncCmp(data.input=data.input,threshold=k, B=B, 
              PermuNo = PermuNo,
              cor.method=cor.method,large.size=large.size,
              scheme_2 = scheme_2)
    
  }
  if(length(result)==0)print('Something wrong')
  
  stopCluster(cl)
  
  gc()
  
  
  clst.ref<-list()
  obs_wise<-list()
  overall<-c()
  mscore<-list()
  mscore.exp<-list()
  
  minvalue<-c()
  obsvalue<-list()
  expvalue<-c()
  ref<-c()
  data<-list()
  graph<-list()
  adjacency<-list()
  if(scheme_2==TRUE){
    data2<-list()
    graph2<-list()
    adjacency2<-list()
  }
  for(i in 1:length(result)){
    
    clst.ref[[i]]<-result[[i]]$stabilityresult$membership
    obs_wise[[i]]<-result[[i]]$stabilityresult$obs_wise
    overall[i]<-result[[i]]$stabilityresult$overall
    ref[i]<-result[[i]]$stabilityresult$ref
    
    mscore[[i]]<-result[[i]]$modularityresult$mscore
    mscore.exp[[i]]<-result[[i]]$modularityresult$mscore.exp
    
    minvalue[i]<-result[[i]]$jaccardresult$min
    obsvalue[[i]]<-result[[i]]$jaccardresult$obsvalue
    expvalue[i]<-result[[i]]$jaccardresult$expvalue
    
    data[[i]]<-result[[i]]$originalinformation$data
    graph[[i]]<-result[[i]]$originalinformation$graph
    adjacency[[i]]<-result[[i]]$originalinformation$adjacency
    
    if(scheme_2==TRUE){
      data2[[i]]<-result[[i]]$originalinformation$data2
      graph2[[i]]<-result[[i]]$originalinformation$graph2
      adjacency2[[i]]<-result[[i]]$originalinformation$adjacency2
      
    }
    
  }
  
  stabilityresult<-list()
  stabilityresult$membership <- clst.ref
  stabilityresult$obs_wise <- obs_wise
  stabilityresult$overall <- overall
  stabilityresult$ref<- ref
  
  modularityresult<-list()
  modularityresult$mscore=mscore
  modularityresult$mscore.exp=mscore.exp
  
  jaccardresult<-list()
  jaccardresult$min=minvalue
  jaccardresult$obsvalue=obsvalue
  jaccardresult$expvalue=expvalue
  
  originalinformation<-list()
  originalinformation$data=data
  originalinformation$graph=graph
  originalinformation$adjacency=adjacency
  if(scheme_2==TRUE){
    originalinformation$data2=data2
    originalinformation$graph2=graph2
    originalinformation$adjacency2=adjacency2
  }
  
  
  result<-list(jaccardresult=jaccardresult,
               stabilityresult=stabilityresult,
               modularityresult=modularityresult,
               originalinformation=originalinformation,
               threshold=threshold.seq)
  
  return(result)
  
}


