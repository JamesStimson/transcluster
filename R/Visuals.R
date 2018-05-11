
##################################
#
# James Stimson 10/04/2018
#
# Visualisations for beyond SNP thresholds
#
##################################

plotClusters <- function(clus, myModel, eWidth=2, vSize=30, vFontSize=2, vColor='cyan'){
  #cgraph <- graph(clus, directed=F)
  edges <- NULL
  wgts <- NULL
  isolates <- myModel$id
  for (i in seq(length(clus)-1)){
    group <- clus[[i]]
    for (j in seq(i+1, length(clus))){
      if (group == clus[[j]]){
        edges <- c(edges, myModel$id[[i]], myModel$id[[j]])
        wgts <- c(wgts, myModel$tcutoff[i,j])
        isolates <- isolates[isolates != myModel$id[[i]]]
        isolates <- isolates[isolates != myModel$id[[j]]]
      }
    }
  }
  #cgraph <- graph(edges=c("A","B"), isolates=myModel$id, directed=F)
  cgraph <- graph(edges=edges, isolates=isolates, directed=F)
  E(cgraph)$weight <- wgts
  plot(cgraph, edge.width=E(cgraph)$weight, vertex.size=vSize, vertex.color=vColor, frame=TRUE, main='Clusters', vertex.label.font=vFontSize)
  return(cgraph)
}

#' interactive plot of SNP clusters with weighted edges
#' @param clus clusters obtained from makeSNPClusters()
#' @param myModel cluster analyis model created using createModel()
#' @param eWidth width of edges
#' @param vSize size of vertices
#' @param vFontSize font size of vertex labels
#' @param vColor vertex colour
#' @param level SNP cut-off level
#' @param thick edge weight adjustment factor
#' @return graph object
#' @export
#' @examples
#' plotSNPClusters(fake, bcModel, level=4, vColor='orange', vSize=4, thick=1.25)
plotSNPClusters <- function(clus, myModel, eWidth=2, vSize=5, vFontSize=0.2, vColor='cyan', level=9, thick=1){
  #cgraph <- graph(clus, directed=F)
  edges <- NULL
  wgts <- NULL
  isolates <- myModel$id
  for (i in seq(length(clus)-1)){
    group <- clus[[i]]
    for (j in seq(i+1, length(clus))){
      if (myModel$snp[i,j]<=level){
        edges <- c(edges, myModel$id[[i]], myModel$id[[j]])
        wgts <- c(wgts, 1+thick*(level-myModel$snp[i,j]))
        isolates <- isolates[isolates != myModel$id[[i]]]
        isolates <- isolates[isolates != myModel$id[[j]]]
      }
    }
  }
  #cgraph <- graph(edges=c("A","B"), isolates=myModel$id, directed=F)
  gtitle <- paste0('SNP based clusters for S = ',level)
  cgraph <- graph(edges=edges, isolates=isolates, directed=F)
  E(cgraph)$weight <- wgts
  plot(cgraph, edge.width=E(cgraph)$weight, vertex.size=vSize, vertex.color=vColor, frame=TRUE, main=gtitle, vertex.label.font=vFontSize, vertex.label.dist=1)
  return(cgraph)
}

#' interactive plot of transmission clusters with weighted edges
#' @param clus clusters obtained from makeTransClusters()
#' @param myModel cluster analyis model created using createModel()
#' @param eWidth width of edges
#' @param vSize size of vertices
#' @param vFontSize font size of vertex labels
#' @param vColor vertex colour
#' @param level transmission SNP cut-off level
#' @param thick edge weight adjustment factor
#' @return graph object
#' @export
#' @examples
#' bcg <- plotTransClusters(fake, bcModel, level=3, vColor='lightblue',vSize=4, thick = 2, vFontSize=1)
plotTransClusters <- function(clus, myModel, eWidth=2, vSize=5, vFontSize=1, vColor='cyan', level=10, thick=1){
  #cgraph <- graph(clus, directed=F)
  edges <- NULL
  wgts <- NULL
  isolates <- myModel$id
  for (i in seq(length(clus)-1)){
    group <- clus[[i]]
    for (j in seq(i+1, length(clus))){
      if (myModel$tcutoff[i,j]<=level){
        edges <- c(edges, myModel$id[[i]], myModel$id[[j]])
        wgts <- c(wgts, 1+ thick*(level-myModel$tcutoff[i,j]))
        isolates <- isolates[isolates != myModel$id[[i]]]
        isolates <- isolates[isolates != myModel$id[[j]]]
      }
    }
  }
  #cgraph <- graph(edges=c("A","B"), isolates=myModel$id, directed=F)
  gtitle <- paste0('Transmission based clusters for T = ',level)
  cgraph <- graph(edges=edges, isolates=isolates, directed=F)
  E(cgraph)$weight <- wgts
  plot(cgraph, edge.width=E(cgraph)$weight, vertex.size=vSize, vertex.color=vColor, frame=TRUE, main=gtitle, vertex.label.font=vFontSize, vertex.label.dist=1)
  return(cgraph)
}

