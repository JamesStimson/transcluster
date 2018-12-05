# Incorporate spatial data into the model

#' set regions from a file
#' @param thisModel Model object
#' @param filename Name of the input date file
#' @return thisModel is returned with the dates added to it
#' @export
setRegionsFromFile <- function(thisModel, filename){

  thisData <- read.csv(filename)
  region_ids <- as.vector(thisData[,1])

  if(!all(region_ids==thisModel$id)){
    print("Warning: IDs in region file don't match those on the model.")
    print(region_ids)
  }

  thisModel$region <- as.vector(thisData[,2])

  # set defaults for colours and weights
  thisModel <- setRegionColours(thisModel)
  thisModel <- setRegionWeights(thisModel)

  thisModel <- setRegionLabels(thisModel, sort(unique(thisModel$region)))

  return (thisModel)
}

#' set region weightings
#' @param thisModel Model object
#' @param sameRegion Weighting for cases from the same region
#' @param diffRegion Weighting for cases from the different regions
#' @export
setRegionWeights <- function(thisModel, sameRegion=1.0, diffRegion=0.1)
{
  thisModel$regwgt <- sameRegion
  thisModel$diffwgt <- diffRegion
  return(thisModel)
}

#' set region colours
#' @param thisModel Model object
#' @param colours Vector of region colours
#' @export
setRegionColours <- function(thisModel, colours=c("orange", "lightblue", "red", "blue", "purple", "lightgreen", "yellow", "green", "brown")){
  thisModel$regcol <- colours
  return(thisModel)
}

#' set region labels
#' @param thisModel Model object
#' @param labels Vector of region labels
#' @export
setRegionLabels <- function(thisModel, labels=c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5", "Region 6")){
  thisModel$reglab <- labels
  return(thisModel)
}

#' return region weights
#' @param thisModel Model object
#' @param i index
#' @param j index
#' @return weight
regionWeight <- function(thisModel, i, j){
  if (thisModel$region[i] == thisModel$region[j]) return (thisModel$regwgt)
  return (thisModel$diffwgt)
}

#' return colour label
#' @param regionLabel region label
#' @return colour label
regionColour <- function(thisModel, regionLabel){
  index <- match(regionLabel, thisModel$reglab)
  if (is.na(index)) return ("gray")
  return(thisModel$regcol[[index]])
}

#' set levels at which number of transmissions passes the threshold, using spatial data
#' @param thisModel the model
#' @param maxK max number of transmissions looked for
#' @export
#' @return the model
setCutoffsSpatial <- function(thisModel, maxK=25){
  thisModel$tcutoff <- matrix(0, nrow(thisModel$snp), ncol(thisModel$snp))
  for (i in seq(1, nrow(thisModel$snp)-1)){
    for (j in seq(i+1, ncol(thisModel$snp))){
      level = nTransCutoff(thisModel$snp[i,j], thisModel$date[i], thisModel$date[j], thisModel$lambda, thisModel$beta, thisModel$perc, maxK, maxK, regionWeight(thisModel, i, j))
      thisModel$tcutoff[i,j] = level
      thisModel$tcutoff[j,i] = level
    }
  }
  return(thisModel)
}
