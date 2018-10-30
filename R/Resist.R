# Incorporate resistance into model
# Instead of the usual work-flow:
# (1) Also run setSNPResistFromFile to load the resistance SNP numbers
# (2) Call makeTransClustersResist instead of makeTransClusters, which in turn calls setCutoffsResist instead of setCutoffs

#' set SNP resistant distance matrix from a file
#' @param thisModel Model object
#' @param filename Name of the input SNP data file
#' @keywords SNP
#' @export
#' @return the model
#' @examples
#' myModel <- setSNPResistFromFile(myModel, snpfile)
setSNPResistFromFile <- function(thisModel, filename){
  thisData <- read.csv(filename)
  thisModel$snpr <- as.matrix(thisData)
  return (thisModel)
}

#' set levels at which number of transmissions passes the threshold
#' @param thisModel the model
#' @param lambda_factor factor by which lambda is increased for resistance SNPs
#' @param maxK max number of transmissions looked for
#' @export
#' @return the model
setCutoffsResist <- function(thisModel, lambda_factor=5, maxK=25){
  thisModel$tcutoff <- matrix(0, nrow(thisModel$snp), ncol(thisModel$snp))
  for (i in seq(1, nrow(thisModel$snp)-1)){
    for (j in seq(i+1, ncol(thisModel$snp))){
      N_zero <- thisModel$snp[i,j] - thisModel$snpr[i,j]
      adjusted_lambda <- thisModel$lambda
      lambda_resist <- lambda_factor*thisModel$lambda
      if (thisModel$snp[i,j]>0 && thisModel$lambda>0 && lambda_resist>0){
        denom = (N_zero/thisModel$lambda) + (thisModel$snpr[i,j]/lambda_resist)
        adjusted_lambda = thisModel$snp[i,j]/denom
      }

      level = nTransCutoff(thisModel$snp[i,j], thisModel$date[i], thisModel$date[j], adjusted_lambda, thisModel$beta, thisModel$perc, maxK, maxK)
      thisModel$tcutoff[i,j] = level
      thisModel$tcutoff[j,i] = level
    }
  }
  return(thisModel)
}

#' make transmission-based clusters for each threshold set in the model using resistance SNPs
#' @param thisModel cluster analysis model created using createModel()
#' @param nameBase basis of name to be used for output cluster files
#' @param writeFile cluster files will be written if TRUE
#' @return list of clusters
#' @keywords cluster
#' @export
#' @examples
#' myTransClusters <- makeTransClustersResist(myModel, 'test')
makeTransClustersResist <- function(thisModel, nameBase='transcluster', writeFile=TRUE){
  allClusters <- NULL
  thisModel <- setCutoffsResist(thisModel)
  for (threshold in thisModel$thTrans){
    clusters = getClustersNonRecursive(threshold, thisModel$tcutoff, thisModel$id)
    if(writeFile){
      writeClusterFile(clusters, threshold, TRUE, nameBase)
    }
    clueClusters <- convertToClueFormat(thisModel$id, clusters)
    allClusters[[length(allClusters)+1]] <- clueClusters
  }
  names(allClusters) <- thisModel$thTrans

  return(allClusters)
}

