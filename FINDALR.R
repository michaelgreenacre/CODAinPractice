FINDALR <- function(data, weight = FALSE) { 
  
  ### -------------------------------------------------------------------
  ### function to identify the best reference for a set of ALRs
  ### various statistics are computed for each reference to assist
  ### the choice of best reference
  ### data is a normalized data matrix
  ### equal weighting is default for the logratio geometry here
  ### row (sample) weighting not implemented in this version
  
  if(sum(data[1,]!=1)) data <- data/rowSums(data)
  
  ### first compute the exact logratio geometry
  data.lra <- LRA(data, weight=weight)
  data.lra.rpc <- data.lra$rowpcoord 
  tot.var <- sum(data.lra$sv^2)
  
  ### loop on all the potential references, computing Procrustes correlation
  ### of each set of ALRs with the exact geometry
  procrust.cor <- rep(0, ncol(data))
  dim <- min(nrow(data), ncol(data)) - 1
  for(j in 1:ncol(data)) {
    # ALR transformation
    alr <- ALR(data, denom=j, weight=weight)
    # ALR geometry using PCA 'by hand' using SVD, without or with weighting
    if(!weight) {
      alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) * sqrt(1/ncol(alr$LR)))
      alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
    }
    if(weight) {
      c <- colMeans(data)
      cc <- c*c[j]
      cc <- cc[-j]
      alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) %*%  diag(sqrt(cc)))
      alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
    }
    procrust.cor[j] <- protest(alr.rpc[,1:dim],data.lra.rpc, permutations=0)$t0
  }

  ### the variances of the log-transformed parts
  var.log <- as.numeric(apply(log(data), 2, var))

  ### highest Procrustes correlation
  procrust.max <- max(procrust.cor)

  ### which reference gives maximum Procrustes
  procrust.ref <- which(procrust.cor==procrust.max)
 
  ### lowest log variance
  var.min <- min(var.log)
  
  ### which reference gives lowest log variance
  var.ref <- which(var.log==var.min)
  
  return(list(tot.var=tot.var, procrust.cor=procrust.cor, 
              procrust.max=procrust.max, procrust.ref=procrust.ref,
              var.log=var.log, var.min=var.min, var.ref=var.ref))
}

