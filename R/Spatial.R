
#' set regions from a file
#' @param thisModel Model object
#' @param filename Name of the input date file
#' @return thisModel is returned with the dates added to it
#' @export
setRegionsFromFile <- function(thisModel, filename){
  # TO DO check order of IDs for data integrity
  thisData <- read.csv(filename)
  thisModel$region <- as.vector(thisData[,2])
  thisModel$subregion <- as.vector(thisData[,3])
  return (thisModel)
}

regionWeight <- function(regions, subregions, i, j){
  # regions <- thisModel$region
  if (subregions[i] == subregions[j]) return (1.0)
  if (regions[i] == regions[j]) return (0.2)
  return (0.1)
}

regionColour <- function(regionLabel){
  if (regionLabel == "Region 1A") return ("orange")
  if (regionLabel == "Region 1B") return ("lightblue")
  if (regionLabel == "Region 2") return ("red")
  if (regionLabel == "Region 3") return ("blue")
  if (regionLabel == "Region 4") return ("purple")
  return ("lightgreen")
}

#' set levels at which number of transmissions passes the threshold
#' @param thisModel the model
#' @export
#' @return the model
setCutoffsSpatial <- function(thisModel){
  print("new version")
  thisModel$tcutoff <- matrix(0, nrow(thisModel$snp), ncol(thisModel$snp))
  for (i in seq(1, nrow(thisModel$snp)-1)){
    for (j in seq(i+1, ncol(thisModel$snp))){
      #adjBeta <- thisModel$beta*regionWeight(thisModel$region, thisModel$subregion, i, j)
      level = nTransCutoff(thisModel$snp[i,j], thisModel$date[i], thisModel$date[j], thisModel$lambda, thisModel$beta, thisModel$perc, 25, 25, regionWeight(thisModel$region, thisModel$subregion, i, j))
      #print(paste0(i,':',j,':',level))
      thisModel$tcutoff[i,j] = level
      thisModel$tcutoff[j,i] = level
    }
  }
  return(thisModel)
}
