##################################
#
# James Stimson 01/03/2018
#
# Methods for beyond SNP thresholds
#
##################################

# ported from python code in SNPMethod.py

#source('Probability.R')

#TESTPATH = "/Users/jamesstimson/cluster/sim/test/"
# TO DO test non-recursive version

writeClusterFile <- function(clusters, threshold, isTrans=FALSE, nameRoot='default', isReal=FALSE){
  label = 'SNP'
  if(isTrans){
    label = 'trans'
  }
  if(isReal){
    label = 'true'
  }
  fileName = paste0(nameRoot, '_cluster_', label, '_', threshold, '.csv')
  cat(paste0(label, ' threshold,', threshold, ',\n'), file=fileName, append=F)
  for (i in seq(1,length(clusters))){
    for (j in seq(1, length(clusters[[i]]))){
      cat(clusters[[i]][[j]], file=fileName, append=T)
      cat(',', file=fileName, append=T)
    }
    cat('\n', file=fileName, append=T)
  }
}

isConnected <- function(threshold, dMatrix, i, j){
  return((dMatrix[i, j] <= threshold)&&(i!=j))
}

addCase <- function (thiscase, thiscluster, dMatrix, mylabels, threshold){
  # Add any other cases connected to this case to the cluster it's in
  for (other_case in seq(1, ncol(dMatrix))){
    if (isConnected(threshold, dMatrix, thiscase, other_case) && !(mylabels[[other_case]] %in% thiscluster)){
      thiscluster <- c(thiscluster, mylabels[[other_case]])
      thiscluster <- addCase(other_case, thiscluster, dMatrix, mylabels, threshold)
    }
  }
  return(thiscluster)
}

addCaseNonRecursive <- function (thiscase, thiscluster, dMatrix, mylabels, threshold){
  # Add any other cases connected to this case to the cluster it's in
  for (other_case in seq(1, ncol(dMatrix))){
    if (isConnected(threshold, dMatrix, thiscase, other_case) && !(mylabels[[other_case]] %in% thiscluster)){
      #print(paste0('adding',mylabels[[other_case]]))
      thiscluster <- c(thiscluster, mylabels[[other_case]])
    }
  }
  return(thiscluster)
}

getClustersNonRecursive <- function(threshold, dMatrix, mylabels, isTrans=FALSE){
  theseClusters <- NULL
  #print('non-recursive test')
  #for (case in seq(1, length(mylabels))){  # this assumed mylabels same as ncol of dMatrix
  for (case in seq(1, ncol(dMatrix))){      # assumes length(mylabels) >= ncol(dMatrix)
    #print(mylabels[[case]])
    found = FALSE
    for (c in theseClusters){
      if (mylabels[case] %in% c){
        found = TRUE
        #print('found')
        break
      }
    }
    #
    if (!found){
      # NEW extra stage IN PROG look for connections to all existing clusters
      ccount <- 0
      tc_index = 0
      for (tc in theseClusters){
        tc_index <- tc_index + 1
        for (element in tc){
          # if match found, add to this cluster
          if (isConnected(threshold, dMatrix, case, which(mylabels==element))){
            #tc <- c(tc, mylabels[[case]])    #add this element to this cluster FAILING, MUST BE A COPY
            theseClusters[[tc_index]] <- c(tc, mylabels[[case]]) # fix, go back to indexed main object
            ccount <- ccount+1
            #print(paste0('adding ',mylabels[[case]],' ', element))
            break
          }
        }
      }
      # if ccount==1 all is well
      if(ccount>1){ # merge clusters
        mindex = -1
        cindex <- 0
        for (c in theseClusters){
          cindex <- cindex+1
          if (mylabels[case] %in% c){
            if(mindex<0) {mindex <- cindex}
            else{
              # merge clusters
              print('Warning: untested code: cluster merging')
              theseClusters[[mindex]] <- c(theseClusters[[mindex]], c)
            }
          }
        }
      }#ccount>1

      if(ccount==0){
        new_cluster = list(mylabels[[case]])
        #print(paste0('new ',mylabels[[case]]))
        new_cluster = addCaseNonRecursive(case, new_cluster, dMatrix, mylabels, threshold)
        # one of the new cases might be in an existing cluster
        new_found = FALSE
        for(newel in new_cluster){
          for (tc in theseClusters){
            if (newel %in% tc){
              new_found = TRUE
              tc <- new_cluster[!newel]
              print('Warning: untested code: new cluster has existing case')
              break
            }
          }
        }
        if(!new_found){
          theseClusters[[length(theseClusters)+1]] <- new_cluster
        }
      }

    }
  }
  #print(theseClusters)
  return(theseClusters)
}


getClusters <- function(threshold, dMatrix, mylabels, isTrans=FALSE){
  theseClusters <- NULL
  #for (case in seq(1, length(mylabels))){  # this assumed mylabels same as ncol of dMatrix
  for (case in seq(1, ncol(dMatrix))){      # assumes length(mylabels) >= ncol(dMatrix)
    found = FALSE
    for (c in theseClusters){
      if (mylabels[case] %in% c){
        found = TRUE
        break
      }
    }
    if (!found){
      new_cluster = list(mylabels[[case]])
      new_cluster = addCase(case, new_cluster, dMatrix, mylabels, threshold)
      theseClusters[[length(theseClusters)+1]] <- new_cluster
    }
  }
  return(theseClusters)
}

#################################
#
# Test code, vignette
#
#################################

#threshold = 7
#x = getClusters(threshold, myModel$tcutoff, mylabels)
#writeClusterFile(x, threshold)

