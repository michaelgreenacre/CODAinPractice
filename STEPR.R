STEPR <- function(data, y, method = NA, family = "gaussian", 
                  nsteps = ncol(data)-1, top = 1, previous = NA, 
                  criterion = "Bonferroni", alpha = 0.05, 
                  previousparts=NA, denom=NA)
### data = compositional data matrix, supposed normalized (but will check)
### y = response variable, can be (i) 0/1 binary (ii) continuous (iii) count
### method = logratio selection method, can be 1=unconstrained, 2=no overlap, 3=ALR/subcomposition
### family = corresponds to y, can be "binomial", "gaussian", "poisson")
### nsteps = number of stepwise selections, default is number of parts minus 1
### top = number of alternative logratios in last step, default = 1. If top > 1, usually nsteps = 1 but not necessarily
### previous = logratios, or other predictors, already in (for one-step-at-a-time selection)
### criterion = stopping criterion, default NA (no stopping), "AIC", "BIC", "Bonferroni"
### alpha = overall significance level for Bonferroni (default 0.05)
### previousparts = the index numbers of the parts already included (in Method 2)
### denom = the denominator to be used in Method 3

{
# preliminaries and error checks
  set.seed(123456789)
  data <- as.matrix(data)
  if (!is.numeric(y) & !is.factor(y)) 
    stop("Response variable neither numeric nor a factor")
  if (!is.factor(y) & family=="binomial")
    stop("Response for binomial family should be a factor")
  if (is.factor(y) & nlevels(y)!=2)
    stop("Response factor must have two levels")
  if (sum(data < 0) > 0) 
    stop("Negative values in compositional data are not allowed")  
  if (sum(data == 0) > 0) 
    stop("Zero values in matrix data on which logratios constructed -- please replace")
  BonValue <- qchisq(1-alpha/(ncol(data)-1), 1)
  if(method == 1) {
    nratios <- nsteps
    deviances   <- matrix(999999, ncol(data), ncol(data))
    logLiks    <- matrix(999999, ncol(data), ncol(data))
    for (j in 2:ncol(data)) {
        for (i in 1:(j - 1)) {
            foo <- as.data.frame(list(logratio=log(data[, i]/data[, j])))
            if (!is.na(previous[[1]][1])) {
                foo <-  as.data.frame(list(previous=previous, logratios=foo))
            }
            foo.glm <- glm(y ~ ., family=family, data=foo)
            deviances[i, j] <- deviance(foo.glm)
            logLiks[i,j]    <- -2*logLik(foo.glm)
        }
    }
    logLik.min <- min(logLiks)
    ratios <- as.matrix(which(logLiks == logLik.min, arr.ind = TRUE))
    logratios <- as.matrix(log(data[, ratios[1, 1]]/data[, ratios[1, 2]]))
    rationames <- paste(colnames(data)[ratios[1, 1]], colnames(data)[ratios[1, 2]], sep = "/")
    rownames(ratios) <- rationames
    colnames(logratios) <- rationames
    predictors <- as.data.frame(list(logratios=logratios))
    npar <- 1
    if (!is.na(previous[[1]][1])) {
        predictors <- as.data.frame(list(previous=previous, logratios=logratios))
        npar <- ncol(predictors)
    }
    if (family == "gaussian") npar <- npar +1
    ratio.glm <- glm(y ~ ., family = family, data = predictors)
    logLik <- -2*logLik(ratio.glm)
    deviance <- deviance(ratio.glm)
    AIC <- logLik + 2 * (npar+1)
    BIC <- logLik+ log(nrow(data)) * (npar+1)
    Bonferroni <- logLik + BonValue * (npar+1)

    if (nratios == 1 & top == 1) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, 
                            null.deviance = ratio.glm$null.deviance))
    }
    if (nratios == 1 & top > 1) {
        logLik.order <- order(logLiks)
        deviance.top <- deviances[logLik.order[1:top]]
        logLik.top  <- logLiks[logLik.order[1:top]]
        ratios.top <- which(logLiks <= logLik.top[top], arr.ind=TRUE)
        ratios.top <- ratios.top[order(logLiks[ratios.top]),] 
        rationames.top <- paste(colnames(data)[ratios.top[, 1]], 
            colnames(data)[ratios.top[, 2]], sep = "/")
        rownames(ratios.top) <- rationames.top
        logratios.top <- log(data[, ratios.top[, 1]]/data[, ratios.top[, 2]])
        colnames(logratios.top) <- rationames.top
        AIC.top <- logLik.top + 2*(npar+1)
        BIC.top <- logLik.top + log(nrow(data)) * (npar+1)
        Bonferroni.top <- logLik.top + BonValue * (npar+1)
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, 
                            null.deviance = ratio.glm$null.deviance,
                            ratios.top = ratios.top, logratios.top = logratios.top, 
                            logLik.top = logLik.top, deviance.top = deviance.top,
                            AIC.top = AIC.top, BIC.top = BIC.top, 
                            Bonferroni.top = Bonferroni.top))
    }
### there are some more steps, nsteps-1 remaining, also test for stopping if specified
### ----------------------------------------------------------------------------------
    for(jratio in 2:nratios) {
        deviances   <- matrix(999999, ncol(data), ncol(data))
        logLiks    <- matrix(999999, ncol(data), ncol(data))
        for (j in 2:ncol(data)) {
            for (i in 1:(j - 1)) {
                foo <- as.data.frame(list(logratios = logratios, logratio = log(data[, i]/data[, j])))
                if (!is.na(previous[[1]][1])) {
                    foo <- as.data.frame(list(previous=previous, logratios=foo))
                }
                foo.glm <- glm(y ~ ., family=family, data = foo)
                deviances[i, j] <- deviance(foo.glm)
                logLiks[i, j]   <- -2*logLik(foo.glm)            }
        }
        logLik.min <- min(logLiks)
### test for multiple solutions
        foo <- as.matrix(which(logLiks == logLik.min, arr.ind = TRUE))
        if(nrow(foo)>1) foo <- foo[sample(1:nrow(foo))[1],]
        ratios <- rbind(ratios, foo)
        rationames <- c(rationames, paste(colnames(data)[ratios[jratio, 1]], colnames(data)[ratios[jratio, 2]], sep = "/"))
        rownames(ratios) <- rationames
        logratios <- cbind(logratios, log(data[, ratios[jratio, 1]]/data[, ratios[jratio, 2]]))
        colnames(logratios) <- rationames
        predictors <- as.data.frame(list(logratios=logratios))
        npar <- jratio
        if (!is.na(previous[[1]][1])) {
            predictors <- as.data.frame(list(previous=previous, logratios=logratios))
            npar <- ncol(predictors)
        }
        if(family == "gaussian") npar<-jratio+1
        ratio.glm <- glm(y ~ ., family = family, data = predictors)
        deviance <- c(deviance, deviance(ratio.glm))
        logLik   <- c(logLik, -2*logLik(ratio.glm))
        AIC <- c(AIC, -2*logLik(ratio.glm)+2*(npar+1))
        BIC <- c(BIC, -2*logLik(ratio.glm)+log(nrow(data))*(npar+1))
        Bonferroni <- c(Bonferroni, -2*logLik(ratio.glm) + BonValue * (npar+1))
### test for stopping
        if(!is.na(criterion)) {
            deciding <- 9999
            if(criterion=="AIC") deciding <- AIC[jratio] - AIC[jratio-1]
            if(criterion=="BIC") deciding <- BIC[jratio] - BIC[jratio-1]
            if(criterion=="Bonferroni") deciding <- Bonferroni[jratio] - Bonferroni[jratio-1]  
            if(deciding == 9999) stop("Stopping criterion has to be one of AIC, BIC, or Bonferroni, in quotation marks")
            if(deciding > 0) {
                print(paste("Criterion increases when ", jratio, "-th ratio enters", sep=""))
            break
            }
        }
    }
    if (top == 1 & is.na(criterion)) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))
    }
    if (top == 1 & !is.na(criterion) & deciding < 0) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))    
    }
    if (top == 1 & !is.na(criterion) & deciding > 0) {
        return(list(names = rationames[-jratio], ratios = ratios[-jratio,], logratios = logratios[,-jratio], 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))    
    }
  }
##-----------------------------------------------------------------------------------------------------------
  if(method == 2) {
    nratios <- min(nsteps, floor(ncol(data)/2))
    rationames <- rep("", nratios)
    deviances   <- matrix(999999, ncol(data), ncol(data))
    logLiks     <- matrix(999999, ncol(data), ncol(data))
    for (j in 2:ncol(data)) {
      if(j %in% previousparts) next
        for (i in 1:(j - 1)) {
          if(i %in% previousparts) next
            foo <- as.data.frame(list(logratio=log(data[, i]/data[, j])))
            if (!is.na(previous[[1]][1])) {
                foo <-  as.data.frame(list(previous=previous, logratios=foo))
            }
            foo.glm <- glm(y ~ ., family=family, data=foo)
            deviances[i, j] <- deviance(foo.glm)
            logLiks[i, j]   <- -2*logLik(foo.glm)
        }
    }
    logLik.min <- min(logLiks)
    ratios <- as.matrix(which(logLiks == logLik.min, arr.ind = TRUE))
    logratios <- as.matrix(log(data[, ratios[1, 1]]/data[, ratios[1, 2]]))
    rationames <- paste(colnames(data)[ratios[1, 1]], colnames(data)[ratios[1, 2]], sep = "/")
    colnames(logratios) <- rationames
    predictors <- as.data.frame(list(logratios=logratios))
    npar <- 1
    if (!is.na(previous[[1]][1])) {
        predictors <- as.data.frame(list(previous=previous, logratios=logratios))
        npar <- ncol(predictors)
    }
    if (family == "gaussian") npar <- npar + 1
    ratio.glm <- glm(y ~ ., family = family, data = predictors)
    logLik   <- -2*logLik(ratio.glm)
    deviance <- deviance(ratio.glm)
    AIC <- logLik + 2 * (npar+1)
    BIC <- logLik + log(nrow(data)) * (npar+1)
    Bonferroni <- logLik + BonValue * (npar+1)
    if (nratios == 1 & top == 1) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, 
                            null.deviance = ratio.glm$null.deviance))
    }
    if (nratios == 1 & top > 1) {
        logLik.order <- order(logLiks)
        logLik.top <- logLiks[logLik.order[1:top]]
        deviance.top  <- deviances[logLik.order[1:top]]
        ratios.top <- which(logLiks <= logLik.top[top], arr.ind=TRUE)
        ratios.top <- ratios.top[order(logLik[ratios.top]),] 
        rationames.top <- paste(colnames(data)[ratios.top[, 1]], 
            colnames(data)[ratios.top[, 2]], sep = "/")
        rownames(ratios.top) <- rationames.top
        logratios.top <- log(data[, ratios.top[, 1]]/data[, ratios.top[, 2]])
        colnames(logratios.top) <- rationames.top
        AIC.top <- logLik.top + 2*npar
        BIC.top <- logLik.top + log(nrow(data)) * npar
        Bonferroni.top <- logLik.top + BonValue * npar
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, 
                            null.deviance = ratio.glm$null.deviance,
                            ratios.top = ratios.top, logratios.top = logratios.top, 
                            logLik.top = logLik.top, deviance.top = deviance.top,
                            AIC.top = AIC.top, BIC.top = BIC.top, 
                            Bonferroni.top = Bonferroni.top))
    }
### there are some more steps, nsteps-1 remaining, also test for stopping if specified
### ----------------------------------------------------------------------------------
    for(jratio in 2:nratios) {
        deviances   <- matrix(999999, ncol(data), ncol(data))
        logLiks    <- matrix(999999, ncol(data), ncol(data))
        iindex <- jindex <- (1:ncol(data))[-c(ratios[,1], ratios[,2])]
        for (j in jindex[2:length(jindex)]) {
          if(j %in% previousparts) next
            for (i in iindex[1:(length(jindex)-1)]) {
              if(i %in% previousparts) next
                foo <- as.data.frame(list(logratios = logratios, logratio = log(data[, i]/data[, j])))
                if (!is.na(previous[[1]][1])) {
                    foo <- as.data.frame(list(previous=previous, logratios=foo))
                }
                foo.glm <- glm(y ~ ., family=family, data = foo)
                deviances[i, j] <- deviance(foo.glm)
                logLiks[i, j]   <- -2*logLik(foo.glm)
            }
        }
        logLik.min <- min(logLiks)
### test for multiple solutions
        foo <- as.matrix(which(logLiks == logLik.min, arr.ind = TRUE))
        if(nrow(foo)>1) foo <- foo[sample(1:nrow(foo))[1],]
        ratios <- rbind(ratios, foo)
        rationames <- c(rationames, paste(colnames(data)[ratios[jratio, 1]], colnames(data)[ratios[jratio, 2]], sep = "/"))
        rownames(ratios) <- rationames
        logratios <- cbind(logratios, log(data[, ratios[jratio, 1]]/data[, ratios[jratio, 2]]))
        colnames(logratios) <- rationames
        predictors <- as.data.frame(list(logratios=logratios))
        npar <- jratio
        if (!is.na(previous[[1]][1])) {
            predictors <- as.data.frame(list(previous=previous, logratios=logratios))
            npar <- ncol(predictors)
        }
        if (family == "gaussian") npar <- jratio + 1
        ratio.glm <- glm(y ~ ., family = family, data = predictors)
        deviance <- c(deviance, deviance(ratio.glm))
        logLik   <- c(logLik, -2*logLik(ratio.glm))
        AIC <- c(AIC, -2*logLik(ratio.glm) + 2 * (npar+1))
        BIC <- c(BIC, -2*logLik(ratio.glm) + log(nrow(data)) * (npar+1))
        Bonferroni <- c(Bonferroni, -2*logLik(ratio.glm) + BonValue * (npar+1))
### test for stopping
        if(!is.na(criterion)) {
            deciding <- 9999
            if(criterion=="AIC") deciding <- AIC[jratio] - AIC[jratio-1]
            if(criterion=="BIC") deciding <- BIC[jratio] - BIC[jratio-1]
            if(criterion=="Bonferroni") deciding <- Bonferroni[jratio] - Bonferroni[jratio-1]  
            if(deciding == 9999) stop("Stopping criterion has to be one of AIC, BIC, or Bonferroni, in quotation marks")
            if(deciding > 0) {
                print(paste("Criterion increases when ", jratio, "-th ratio enters", sep=""))
            break
            }
        }
    }
    if (top == 1 & is.na(criterion)) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))
    }
    if (top == 1 & !is.na(criterion) & deciding < 0) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))    
    }
    if (top == 1 & !is.na(criterion) & deciding > 0) {
        return(list(names = rationames[-jratio], ratios = ratios[-jratio,], logratios = logratios[,-jratio], 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))    
    }

  
  }
##--------------------------------------------------------------------------------------------------------
  if(method == 3) {
    nratios <- min(nsteps, floor(ncol(data)/2))
    rationames <- rep("", nratios)
    deviances   <- matrix(999999, ncol(data), ncol(data))
    logLiks     <- matrix(999999, ncol(data), ncol(data))
    if(is.na(denom)) {
        for (j in 2:ncol(data)) {
            for (i in 1:(j - 1)) {
                foo <- as.data.frame(list(logratio=log(data[, i]/data[, j])))
                if (!is.na(previous[[1]][1])) {
                    foo <-  as.data.frame(list(previous=previous, logratios=foo))
                }
                foo.glm <- glm(y ~ ., family=family, data=foo)
                deviances[i, j] <- deviance(foo.glm)
                logLiks[i, j]    <- -2*logLik(foo.glm)
            }
        }
    }
    if(!is.na(denom)) {
      
        for (i in (1:ncol(data))[-denom]) {
          foo <- as.data.frame(list(logratio=log(data[, i]/data[, denom])))
          if (!is.na(previous[[1]][1])) {
            foo <-  as.data.frame(list(previous=previous, logratios=foo))
          }
          foo.glm <- glm(y ~ ., family=family, data=foo)
          deviances[i, denom] <- deviance(foo.glm)
          logLiks[i, denom]    <- -2*logLik(foo.glm)
        }
      
     }
    logLik.min <- min(logLiks)
    ratios <- as.matrix(which(logLiks == logLik.min, arr.ind = TRUE))
    logratios <- as.matrix(log(data[, ratios[1, 1]]/data[, ratios[1, 2]]))
    rationames <- paste(colnames(data)[ratios[1, 1]], colnames(data)[ratios[1, 2]], sep = "/")
    colnames(logratios) <- rationames
    predictors <- as.data.frame(list(logratios=logratios))
    npar <- 1
    if (!is.na(previous[[1]][1])) {
        predictors <- as.data.frame(list(previous=previous, logratios=logratios))
        npar <- ncol(predictors)
    }
    if (family == "gaussian") npar <- npar + 1
    ratio.glm <- glm(y ~ ., family = family, data = predictors)
    logLik   <- -2*logLik(ratio.glm)
    deviance <- deviance(ratio.glm)
    AIC <- -2*logLik(ratio.glm) + 2 * (npar+1)
    BIC <- -2*logLik(ratio.glm) + log(nrow(data)) * (npar+1)
    Bonferroni <- -2*logLik(ratio.glm) + BonValue * (npar+1)
    if (nratios == 1 & top == 1) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni,
                            null.deviance = ratio.glm$null.deviance))
    }
    if (nratios == 1 & top > 1) {
        logLik.order <- order(logLiks)
        logLik.top  <- logLiks[logLik.order[1:top]]
        deviance.top <- deviances[logLik.order[1:top]]
        ratios.top <- which(logLiks <= logLik.top[top], arr.ind=TRUE)
        ratios.top <- ratios.top[order(logLiks[ratios.top]),] 
        rationames.top <- paste(colnames(data)[ratios.top[, 1]], 
            colnames(data)[ratios.top[, 2]], sep = "/")
        rownames(ratios.top) <- rationames.top
        logratios.top <- log(data[, ratios.top[, 1]]/data[, ratios.top[, 2]])
        colnames(logratios.top) <- rationames.top
        AIC.top <- logLik.top + 2 * (npar+1)
        BIC.top <- logLik.top + log(nrow(data)) * (npar+1)
        Bonferroni.top <- logLik.top + BonValue * (npar+1)
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, 
                            null.deviance = ratio.glm$null.deviance, 
                            ratios.top = ratios.top, logratios.top = logratios.top, 
                            logLik.top = logLik.top, deviance.top = deviance.top,
                            AIC.top = AIC.top, BIC.top = BIC.top, 
                            Bonferroni.top = Bonferroni.top))
    }
### there are some more steps, nsteps-1 remaining, also test for stopping if specified
### ----------------------------------------------------------------------------------
    for(jratio in 2:nratios) {
        deviances   <- matrix(999999, ncol(data), ncol(data))
        logLiks     <- matrix(999999, ncol(data), ncol(data))
        iindex <- (1:ncol(data))[-c(ratios[,1], ratios[1,2])]
        for (i in iindex) {
                foo <- as.data.frame(list(logratios = logratios, logratio = log(data[, i]/data[, ratios[1,2]])))
                if (!is.na(previous[[1]][1])) {
                    foo <- as.data.frame(list(previous=previous, logratios=foo))
                }
                foo.glm <- glm(y ~ ., family=family, data = foo)
            deviances[i, ratios[1,2]] <- deviance(foo.glm)
            logLiks[i, ratios[1,2]]   <- -2*logLik(foo.glm)
        }
        logLik.min <- min(logLiks)
### test for multiple solutions
        foo <- as.matrix(which(logLiks == logLik.min, arr.ind = TRUE))
        if(nrow(foo)>1) foo <- foo[sample(1:nrow(foo))[1],]
        ratios <- rbind(ratios, foo)
        rationames <- c(rationames, paste(colnames(data)[ratios[jratio, 1]], colnames(data)[ratios[jratio, 2]], sep = "/"))
        rownames(ratios) <- rationames
        logratios <- cbind(logratios, log(data[, ratios[jratio, 1]]/data[, ratios[jratio, 2]]))
        colnames(logratios) <- rationames
        predictors <- as.data.frame(list(logratios=logratios))
        npar <- jratio
        if (!is.na(previous[[1]][1])) {
            predictors <- as.data.frame(list(previous=previous, logratios=logratios))
            npar <- ncol(predictors)
        }
        if (family == "gaussian") npar <- npar + 1
        ratio.glm <- glm(y ~ ., family = family, data=predictors)
        logLik   <- c(logLik, -2*logLik(ratio.glm))
        deviance <- c(deviance, deviance(ratio.glm))
        AIC <- c(AIC, -2*logLik(ratio.glm) + 2*(npar+1))
        BIC <- c(BIC, -2*logLik(ratio.glm) + log(nrow(data)) * (npar+1))
        Bonferroni <- c(Bonferroni, -2*logLik(ratio.glm) + BonValue * (npar+1))
### test for stopping
        if(!is.na(criterion)) {
            deciding <- 9999
            if(criterion=="AIC") deciding <- AIC[jratio] - AIC[jratio-1]
            if(criterion=="BIC") deciding <- BIC[jratio] - BIC[jratio-1]
            if(criterion=="Bonferroni") deciding <- Bonferroni[jratio] - Bonferroni[jratio-1]  
            if(deciding == 9999) stop("Stopping criterion has to be one of AIC, BIC, or Bonferroni, in quotation marks")
            if(deciding > 0) {
                print(paste("Criterion increases when ", jratio, "-th ratio enters", sep=""))
            break
            }
        }
    }
    if (top == 1 & is.na(criterion)) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))
    }
    if (top == 1 & !is.na(criterion) & deciding < 0) {
        return(list(names = rationames, ratios = ratios, logratios = logratios, 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))    
    }
    if (top == 1 & !is.na(criterion) & deciding > 0) {
        return(list(names = rationames[-jratio], ratios = ratios[-jratio,], logratios = logratios[,-jratio], 
                            logLik = logLik, deviance = deviance, AIC = AIC, BIC = BIC, 
                            Bonferroni = Bonferroni, null.deviance = ratio.glm$null.deviance))    
    }
  
  }

}


