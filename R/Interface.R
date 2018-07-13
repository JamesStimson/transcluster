##################################
#
# James Stimson 05/03/18
# Interface for loading data into Model
#
##################################

#' set sample dates from a file
#' @param thisModel Model object
#' @param filename Name of the input date file
#' @return thisModel is returned with the dates added to it
#' @export
#' @examples
#' myModel <- setDatesFromFile(myModel, datefile)
setDatesFromFile <- function(thisModel, filename){
  #thisData <- read.csv(system.file("extdata", filename, package = "transcluster", mustWork = TRUE))
  thisData <- read.csv(filename)
  print(thisData[,1])
  thisModel$id <- as.vector(thisData[,1])
  thisModel$date <- as.vector(thisData[,2]) # date data is assumed to be as number, eg 2015.32

  # patch up missing dates
  maxdate = 0
  for (d in thisModel$date){
    if (is.na(d)){
      # do nothing
    }
    else if (d > maxdate){
      maxdate <- d
    }
  }
  thisModel$date[is.na(thisModel$date)] <- maxdate
  return (thisModel)
}


#' set SNP distance matrix from a file
#' @param thisModel Model object
#' @param filename Name of the input SNP data file
#' @keywords SNP
#' @export
#' @examples
#' myModel <- setSNPFromFile(myModel, snpfile)
setSNPFromFile <- function(thisModel, filename){
  #thisData <- read.csv(system.file("extdata", filename, package = "transcluster", mustWork = TRUE))
  thisData <- read.csv(filename)
  thisModel$snp <- as.matrix(thisData)
  return (thisModel)
}

#' set transmission distance matrix from a file
#' @param thisModel Model object
#' @param filename Name of the input transmission data file
#' @keywords transmission
#' @export
#' @examples
#' myModel <- setTransMatrixFromFile(myModel, transfile)
setTransMatrixFromFile <- function(thisModel, filename){
  thisData <- read.csv(filename)
  thisModel$trans <- as.matrix(thisData[,2:ncol(thisData)])
  return (thisModel)
}

#' compare similarity of clusters
#' @param SNPClusters List of clusters
#' @param transClusters List of clusters
#' @param roundTo Number of decimal places in return matrix
#' @export
compareClusters <- function(SNPClusters, transClusters, roundTo=3){
  # compare results using clue and return as a matrix
  compRes <- matrix(0, length(SNPClusters), length(transClusters))
  for (i in seq(length(SNPClusters))){
    for (j in seq(length(transClusters))){
      pSNP <- as.cl_partition(SNPClusters[[i]])
      pTrans <- as.cl_partition(transClusters[[j]])
      diss <- cl_dissimilarity(pSNP, pTrans, method="VI")
      compRes[i, j] <- round(diss, roundTo)
    }
  }
  return(compRes)
}

makeClustersFromSims <- function(thisModel, numSims, dateFileRoot, snpFileRoot, transFileRoot, subdir, writeFile=FALSE){

  sumSNPDiff = 0.0
  sumTransDiff = 0.0

  for (i in seq(1, numSims)){
    print (paste0('processing sim ', i))
    thisModel <- setDatesFromFile(thisModel, paste0(dateFileRoot, 'dates', i, '.csv'))
    #print (paste0(length(thisModel$date),' dates set'))
    thisModel <- setSNPFromFile(thisModel, paste0(snpFileRoot, i, '.csv'))
    #print (paste0(length(thisModel$snp),' snps set'))
    thisModel <- setTransMatrixFromFile(thisModel, paste0(transFileRoot, i, '.csv'))
    #print (' trans matrix set')
    SNPClusters <- makeSNPClusters(thisModel, paste0(subdir, i), writeFile)
    #print ('snp clusters made')
    knownTransClusters <- makeKnownTransClusters(thisModel, paste0(subdir, i), writeFile)
    transClusters <- makeTransClusters(thisModel, paste0(subdir, i, '_', thisModel$lambda, '_', thisModel$beta), writeFile)

    # compare results using clue
    # create a vector of clusters sizes so have baseline for comparison
    knownSizes <- NULL
    for (j in seq(length(knownTransClusters))){
      knownSizes[[j]] <- max(knownTransClusters[[j]])
    }
    level <- knownSizes[[1]]  # Direct transmissions
    SNP_index <- -1
    trans_index <- -1
    SNPSizes <- NULL
    for (j in seq(length(SNPClusters))){
      SNPSizes[[j]] <- max(SNPClusters[[j]])
      if(SNPSizes[[j]] <= level){
        SNP_index <- j
        break
      }
    }

    transSizes <- NULL
    for (j in seq(length(transClusters))){
      transSizes[[j]] <- max(transClusters[[j]])
      if(transSizes[[j]] <= level){
        trans_index <- j
        break
      }
    }

    pKnown <- as.cl_partition(knownTransClusters[[1]])
    pSNP <- as.cl_partition(SNPClusters[[SNP_index]])
    pTrans <- as.cl_partition(transClusters[[trans_index]])

    dissSNP <- cl_dissimilarity(pKnown, pSNP, method="VI")
    dissTrans <- cl_dissimilarity(pKnown, pTrans, method="VI")

    print(paste0(i, ' ',  dissSNP, ' ', dissTrans, ' ', SNP_index, ' ', trans_index))
    # low dissimilarity is good
    sumSNPDiff = sumSNPDiff + dissSNP
    sumTransDiff = sumTransDiff + dissTrans
  }
  print(sumSNPDiff)
  print(sumTransDiff)
}

# For testing makeKnownTransClusters:
makeTransClustersFromSims <- function(thisModel, numSims, dateFileRoot, snpFileRoot, transFileRoot){
  for (i in seq(1, numSims)){
    print (paste0('processing sim ', i))
    thisModel <- setDatesFromFile(thisModel, paste0(dateFileRoot, 'dates', i, '.csv'))
    thisModel <- setSNPFromFile(thisModel, paste0(snpFileRoot, i, '.csv'))
    thisModel <- setTransMatrixFromFile(thisModel, paste0(transFileRoot, i, '.csv'))
    knownTransClusters <- makeKnownTransClusters(thisModel, paste0(dateFileRoot, i))
  }
}
