#' @name network.stability.output
#'
#' @title Plot method for objests from threshold.select
#'
#' @details \code{network.stability.output} is used to generate a series of network plots based on the given \code{threshold.seq},where the nodes are
#' colored by the level of stability. The network with optimal
#' threshold value selected by function \code{threshold.select} is colored as red.
#' 
#' @param input a \code{list} of results from function \code{threshold.select}
#' @param optimal.only a \code{logical} value indicating whether only plot the network with optimal threshold or not. The default is False, generating
#' all network figures with a large number of nodes could take some time.
#'
#' @return Plot of network figures
#'
#' @author Mingmei Tian
#'
#' @references A framework for stability-based module detection in correlation graphs.
#' Mingmei Tian,Rachael Hageman Blair,Lina Mu, Matthew Bonner, Richard Browne and Han Yu.
#' @importFrom grid unit.c
#' @import ggplot2 
#' @import gridExtra 
#' @import intergraph 
#' @import GGally
#' @importFrom network set.vertex.attribute
#' 
#'
#' @examples
#' \dontrun{
#' network.stability.output(result)
#' }
#' @export
network.stability.output<-function(input,optimal.only=FALSE){

  test.data<-input$originalinformation$data[[1]]
  test.value<-do.call(rbind,input$jaccardresult$obsvalue)
  test.mean<-rowMeans(test.value,na.rm = TRUE)
  reference.value<-input$jaccardresult$expvalue
  B<-ncol(test.value)-1
  substract.value<- test.value- matrix(rep(reference.value,B+1),ncol =(B+1))
  substract.mean<-rowMeans(substract.value,na.rm = TRUE)
  if(optimal.only){
    i=which.max(substract.mean)
    test.graph<-input$originalinformation$graph[[i]]
    test.adjacency<-input$originalinformation$adjacency[[i]]
    
    node_stability<-input$stabilityresult$obs_wise[[i]]
    node_color<- cut(node_stability ,
                     breaks=c(-0.1,0.2,0.4 ,0.6 ,0.8,1.1), 
                     labels=c("0.0-0.2","0.2-0.4","0.4-0.6",'0.6-0.8',"0.8-1.0"))
    
    net <- asNetwork(test.graph)
    set.vertex.attribute(net,"Node Stability",as.vector(node_color))
    #net %v% "Node Stability"<-as.vector(node_color)
    plot_list<-ggnet2(net,color = "Node Stability",mode = "kamadakawai",
                      size=3,edge.color = "grey",
                      palette =c("0.0-0.2"="#440154FF","0.2-0.4"="#3B528BFF", 
                                 "0.4-0.6"="#21908CFF" ,
                                 "0.6-0.8"="#5DC863FF","0.8-1.0"="#FDE725FF"))+
      ggtitle(paste0('Threshold=',input$threshold[i]))+
      theme(plot.title = element_text(color = "red"))
  }else{
    plot_list<-NULL
    for(i in 1:length(input$threshold)){
      
      test.graph<-input$originalinformation$graph[[i]]
      test.adjacency<-input$originalinformation$adjacency[[i]]
      
      node_stability<-input$stabilityresult$obs_wise[[i]]
      node_color<- cut(node_stability ,
                       breaks=c(-0.1,0.2,0.4 ,0.6 ,0.8,1.1), 
                       labels=c("0.0-0.2","0.2-0.4","0.4-0.6",'0.6-0.8',"0.8-1.0"))
      
      net <- asNetwork(test.graph)
      #net %v% "Node Stability"<-as.vector(node_color)
      set.vertex.attribute(net,"Node Stability",as.vector(node_color))
      
      
      if(i==which.max(substract.mean))  {
        plot_list[[i]]<-ggnet2(net,color = "Node Stability",mode = "kamadakawai",
                               size=3,edge.color = "grey",
                               palette =c("0.0-0.2"="#440154FF","0.2-0.4"="#3B528BFF", 
                                          "0.4-0.6"="#21908CFF" ,
                                          "0.6-0.8"="#5DC863FF","0.8-1.0"="#FDE725FF"))+
          ggtitle(paste0('Threshold=',input$threshold[i]))+
          theme(plot.title = element_text(color = "red"))
      }else{
        plot_list[[i]]<-ggnet2(net,color = "Node Stability",mode = "kamadakawai",size=3,edge.color = "grey",
                               palette =c("0.0-0.2"="#440154FF","0.2-0.4"="#3B528BFF", "0.4-0.6"="#21908CFF" ,
                                          "0.6-0.8"="#5DC863FF","0.8-1.0"="#FDE725FF"))+
          ggtitle(paste0('Threshold=',input$threshold[i]))
      }
      
      
      
    }
    
    
  }
  # Data for legend
  i=length(input$threshold)
  test.graph<-input$originalinformation$graph[[i]]
  test.adjacency<-input$originalinformation$adjacency[[i]]
  
  node_stability<-input$stabilityresult$obs_wise[[i]]
  #range(node_stability)
  node_color<- cut(node_stability ,
                   breaks=seq(range(node_stability)[1]-0.01,range(node_stability)[2]+0.01,length.out = 6), 
                   labels=c("0.0-0.2","0.2-0.4","0.4-0.6",'0.6-0.8',"0.8-1.0"))
  
  net <- asNetwork(test.graph)
  
  set.vertex.attribute(net,"Node Stability",as.vector(node_color))
  
  plot_list[[(length(input$threshold))+1]]<-ggnet2(net,color = "Node Stability",mode = "kamadakawai",size=3,edge.color = "grey",
                                                   palette =c("0.0-0.2"="#440154FF","0.2-0.4"="#3B528BFF", "0.4-0.6"="#21908CFF" ,
                                                              "0.6-0.8"="#5DC863FF","0.8-1.0"="#FDE725FF"))+
    ggtitle(paste0('Threshold=',input$threshold[i]))
  UseMethod("plot")
  
  # Final plot 
  
  return(grid_arrange_shared_legend(plot_list))
  
}