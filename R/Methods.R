##################################
#
# James Stimson 01/03/2018
#
# Methods for beyond SNP thresholds
#
##################################

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
      thiscluster <- c(thiscluster, mylabels[[other_case]])
    }
  }
  return(thiscluster)
}

getClustersNonRecursive <- function(threshold, dMatrix, mylabels, isTrans=FALSE){
  theseClusters <- NULL
  for (case in seq(1, ncol(dMatrix))){      # assumes length(mylabels) >= ncol(dMatrix)
    found = FALSE
    for (c in theseClusters){
      if (mylabels[case] %in% c){
        found = TRUE
        break
      }
    }
    if (!found){
      # look for connections to all existing clusters
      ccount <- 0
      tc_index = 0
      for (tc in theseClusters){
        tc_index <- tc_index + 1
        for (element in tc){
          # if match found, add to this cluster
          if (isConnected(threshold, dMatrix, case, which(mylabels==element))){
            theseClusters[[tc_index]] <- c(tc, mylabels[[case]])
            ccount <- ccount+1
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
              # Need to (a) remove this case from c (b) add c to mindex cluster (c) delete c, by setting to NULL
              c[[which(c==mylabels[case])]] <- NULL
              theseClusters[[mindex]] <- c(theseClusters[[mindex]], c)
              theseClusters[[cindex]] <- NULL
            }
          }
        }
      }#ccount>1

      if(ccount==0){
        new_cluster = list(mylabels[[case]])
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
  return(theseClusters)
}


getClusters <- function(threshold, dMatrix, mylabels, isTrans=FALSE){
  theseClusters <- NULL
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


