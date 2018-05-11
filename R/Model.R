##################################
#
# James Stimson 28/02/2018
#
# Model for beyond SNP thresholds
#
##################################

# Model list contains the data for the project
# 1. Create the model with at least minimal inputs
# 2. Run one of the functions to get outputs

# Inputs for thresholds can be vectors of integers or just one integer

# 21/03/18 return something in clue-friendly format when writing files
# TO consider: add clusters to model and return whole model when making clusters? A: return as list if the two things
#setwd("~/cluster")
#source('Methods.R')

ProbHeight = function(h, N, lambda, dd, dterm_list){
  psum <- 0.0

  for (i in seq(0, N))
  {
    hterm = lambda^{i+1} * h^i * exp(-lambda*h)/factorial(i)
    dterm = dterm_list[[i+1]]
    psum <- psum + hterm*dterm
  }
  return(psum)
}

# Gamma version of ProbHeight() -TO DO: Wire this up-
ProbHeightGamma = function(h, N, a, b, dd, dterm_list){
  psum <- 0.0

  for (i in seq(0, N))
  {
    hterm = b^{a*(i+1)} * h^{a*(i+1)-1} * exp(-b*h)/gamma(a*(i+1))
    dterm = dterm_list[[i+1]]
    psum <- psum + hterm*dterm
  }
  return(psum)
}

ProbTrans = function(t, k, beta){
  return({beta*t}^k * exp(-beta*t)/factorial(k))
}

Integrand = function(h, delta_time, N, k, lambda, beta, dterm_list){
  return(ProbHeight(h, N, lambda, delta_time, dterm_list) * ProbTrans(h+delta_time, k, beta))
}

ProbKTransmissions = function(N, k, t1, t2, lambda, beta){
  upper_limit = Inf
  delta_time = abs(t1-t2)
  dsum <- 0.0
  dterm_list = list()
  for (i in seq(0, N)){
    dterm = {lambda*delta_time}^{N-i} * exp(-lambda*delta_time)/factorial(N-i)
    dsum = dsum + dterm
    dterm_list[[i+1]] = dterm
  }
  result = tryCatch(
    stats::integrate(Integrand, lower=0, upper=upper_limit, delta_time,N,k,lambda,beta,dterm_list),
    error=function(error_message) {
      message("Integration error in ProbKTransmissions().")
      print(paste0(t1,':',t2))
      message(error_message)
      return(NA)
    }
  )
  return(result$val/dsum)
}

#' create levels at which number of transmissions passes the threshold
#' @param N SNP distance
#' @param t1 date of first sample
#' @param t2 date of second sample
#' @param lambda clock rate
#' @param beta transmission rate
#' @param perc_cutoff cut-off percentage
#' @param maxN for optimisation, SNP distance greater than this considered not in same cluster
#' @return number of transmissions reached
#' @export
nTransCutoff = function(N, t1, t2, lambda, beta, perc_cutoff, maxN=25, maxK=25){
  if (N > maxN) {return (maxK)}
  total_prob = 0.0
  for (k in seq(0, maxK)){
    this_prob = ProbKTransmissions(N, k, t1, t2, lambda, beta)
    total_prob = total_prob + this_prob
    if (total_prob > perc_cutoff) {return (k)}
  }
  return (maxK)
}

#' test whether two cases are within the same transmission cluster
#' @param N SNP distance
#' @param k_cutoff number of transmissions below this threshold considered to be in same cluster
#' @param t1 date of first sample
#' @param t2 date of second sample
#' @param lambda clock rate
#' @param beta transmission rate
#' @param perc_cutoff cut-off percentage
#' @param maxN for optimisation, SNP distance greater than this considered not in same cluster
#' @param N_resist number of resistance-conferring SNPs as a subset of N
#' @param lambda_resist clock rate for resistance-conferring SNPs
#' @return TRUE if within cluster, FALSE otherwise
#' @export
#' @examples
#' isWithinCluster <- WithinCluster(N=5, k_cutoff=3, t1=2012.13, t2=2015.76, lambda=0.7, beta=1.3, perc_cutoff=0.2)
withinCluster = function(N, k_cutoff, t1, t2, lambda, beta, perc_cutoff, maxN=25, N_resist=0, lambda_resist=0.0){
  if (N > maxN) {return (FALSE)}
  total_prob = 0.0
  for (k in seq(0, k_cutoff)){
    this_prob = ProbKTransmissions(N, k, t1, t2, lambda, beta)
    total_prob = total_prob + this_prob
    if (total_prob > perc_cutoff) {return (TRUE)}
  }
  return (FALSE)
}



#' create model to be used for cluster analysis
#' @param ids (optional) vector of unique sample IDs
#' @param dates (optional) vector of dates (formatted as numbers)
#' @param snpMatrix (optional) symmetric matrix of SNP distances supplied in the same order as ids
#' @param lambda (optional) clock rate
#' @param beta (optional) transmission rate
#' @param percCutOff (optional) cut-off percentage
#' @param thSNP (optional) SNP threshold values as a vector of integers
#' @param thTrans (optional) transmission threshold values as a vector of integers
#' @return model object for further analysis
#' @export
#' @examples
#' myModel <- createModel()
createModel <- function(ids=NULL, dates=NULL, snpMatrix=matrix(,0,0), lambda=1.0, beta=1.0, percCutOff=0.2, thSNP=seq(1,10), thTrans=seq(1,10)){

  if (!(length(ids)==length(dates) && length(ids)==nrow(snpMatrix) && length(ids)==ncol(snpMatrix))){
    print("Inconsistent inputs")
    return()
  }

  tcutoffMatrix = matrix(0, nrow=nrow(snpMatrix), ncol=ncol(snpMatrix))

  thisModel = list(id=ids, date=dates, snp=snpMatrix, tcutoff=tcutoffMatrix, lambda=lambda, beta=beta, perc=percCutOff, thSNP=thSNP, thTrans=thTrans)

  return(thisModel)
}

#' set parameters for existing cluster analysis model
#' @param thisModel cluster analyis model created using createModel()
#' @param lambda (optional) clock rate
#' @param beta (optional) transmission rate
#' @param percCutOff (optional) cut-off percentage
#' @return model object
#' @export
#' @examples
#' myModel <- setParams(myModel, lambda=0.5)
setParams <- function(thisModel, lambda=-1, beta=-1, percCutOff=-1){
  if (lambda>0)     {thisModel$lambda = lambda}
  if (beta>0)       {thisModel$beta = beta}
  if (percCutOff>0) {thisModel$perc = percCutOff}
  return(thisModel)
}

#' set SNP thresholds for existing cluster analysis model
#' @param thisModel cluster analsyis model created using createModel()
#' @param thSNP (optional) SNP threshold values as a vector of integers
#' @return model object
#' @export
#' @examples
#' myModel <- setSNPThresholds(myModel, seq(12))
setSNPThresholds <- function(thisModel, thSNP=seq(1,10)){
  if (length(thSNP)>0){
    thisModel$thSNP <- thSNP
  }
  else{
    print("Invalid input")
  }
  return(thisModel)
}

#' set transmission thresholds for existing cluster analysis model
#' @param thisModel cluster analysis model created using createModel()
#' @param thTrans (optional) transmission threshold values as a vector of integers
#' @return model object
#' @export
#' @examples
#' myModel <- setTransThresholds(myModel, seq(12))
setTransThresholds <- function(thisModel, thTrans=seq(1,10)){
  if (length(thTrans)>0){
    thisModel$thTrans <- thTrans
  }
  else{
    print("Invalid input")
  }
  return(thisModel)
}

# Requires index of ids, hence use of which
convertToClueFormat <- function(ids, clusters){
  grouping <- NULL
  for (i in seq(1,length(clusters))){
    for (j in seq(1, length(clusters[[i]]))){
      grouping[which(ids==clusters[[i]][[j]])] <- i
    }
  }
  return(grouping)
}

#' make SNP-based clusters for each threshold set in the model
#' @param thisModel cluster analysis model created using createModel()
#' @param nameRoot (optional) basis of name to be used for output cluster files
#' @param writeFile (optional) cluster files will be written if TRUE
#' @return list of clusters
#' @export
#' @examples
#' mySNPClusters <- makeSNPClusters(myModel, 'test')
makeSNPClusters <- function(thisModel, nameRoot='SNPcluster', writeFile=TRUE){
  allClusters <- NULL
  for (threshold in thisModel$thSNP){
    clusters <- getClustersNonRecursive(threshold, thisModel$snp, thisModel$id) #TESTING
    if(writeFile){
      writeClusterFile(clusters, threshold, FALSE, nameRoot)
    }
    clueClusters <- convertToClueFormat(thisModel$id, clusters)
    allClusters[[length(allClusters)+1]] <- clueClusters
  }
  names(allClusters) <- thisModel$thSNP
  return(allClusters)
}

makeKnownTransClusters <- function(thisModel, nameRoot, writeFile=TRUE){
  allClusters <- NULL
  for (threshold in thisModel$thTrans){
    # Need to cut the trans matrix down to snp matrix size
    clusters = getClustersNonRecursive(threshold, thisModel$trans[1:ncol(thisModel$snp), 1:ncol(thisModel$snp)], thisModel$id) #TESTING
    if(writeFile){
      writeClusterFile(clusters, threshold, FALSE, nameRoot, TRUE)
    }
    clueClusters <- convertToClueFormat(thisModel$id, clusters)
    allClusters[[length(allClusters)+1]] <- clueClusters
  }
  names(allClusters) <- thisModel$thTrans
  return(allClusters)
}

#' set levels at which number of transmissions passes the threshold
#' @param thisModel the model
#' @return the model
setCutoffs <- function(thisModel){
  thisModel$tcutoff <- matrix(0, nrow(thisModel$snp), ncol(thisModel$snp))
  for (i in seq(1, nrow(thisModel$snp)-1)){
    for (j in seq(i+1, ncol(thisModel$snp))){
      level = nTransCutoff(thisModel$snp[i,j], thisModel$date[i], thisModel$date[j], thisModel$lambda, thisModel$beta, thisModel$perc, 25, 25)
      #print(paste0(i,':',j,':',level))
      thisModel$tcutoff[i,j] = level
      thisModel$tcutoff[j,i] = level
    }
  }
  return(thisModel)
}

#' make transmission-based clusters for each threshold set in the model
#' @param thisModel cluster analysis model created using createModel()
#' @param nameRoot (optional) basis of name to be used for output cluster files
#' @param writeFile (optional) cluster files will be written if TRUE
#' @return list of clusters
#' @keywords cluster
#' @export
#' @examples
#' myTransClusters <- makeTransClusters(myModel, 'test')
makeTransClusters <- function(thisModel, nameRoot='transcluster', writeFile=TRUE){
  allClusters <- NULL
  thisModel <- setCutoffs(thisModel)
  for (threshold in thisModel$thTrans){
    clusters = getClustersNonRecursive(threshold, thisModel$tcutoff, thisModel$id)
    if(writeFile){
      writeClusterFile(clusters, threshold, TRUE, nameRoot)
    }
    clueClusters <- convertToClueFormat(thisModel$id, clusters)
    allClusters[[length(allClusters)+1]] <- clueClusters
  }
  names(allClusters) <- thisModel$thTrans
  #####################
  # new code 19/04/2018

  #####################

  return(allClusters)
}

#################################
#
# Test code, vignette
#
#################################

DOTEST = TRUE

if(DOTEST){
  ids = c('A', 'B', 'C', 'D')
  dates = c(2018, 2014, 2016, 2016)
  snpMatrix = matrix(c(0,5,7,7,5,0,8,8,7,8,0,6,7,8,6,0),nrow=4,ncol=4)

  myModel <- createModel(ids, dates, snpMatrix)
  myModel <- setParams(myModel, lambda=1.5)
  myModel <- setParams(myModel, beta=2.3)
  myModel <- setSNPThresholds(myModel, c(5))
  myModel <- setTransThresholds(myModel, c(7,8,9))

  mySNPClusters = makeSNPClusters(myModel, 'test2')
  myTransClusters = makeTransClusters(myModel, 'test2')
}

