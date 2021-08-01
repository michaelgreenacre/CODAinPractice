### RNAseq data of Jovanovic et al. (2015)
### Files obtained from Quinn et al. (2019), repository http://doi.org/10.5281/zenodo.3270954
mice <- read.csv("rnaseq-x.csv", row.names=1)
mice.annot <- read.csv("rnaseq-y.csv", row.names=1)
mice <- t(mice)
dim(mice)
# [1]  28 3147

### this data set has some zeros that are not exactly zero!
sum(mice<1e-10)
# [1] 34
### so redefine them as zeros
mice[which(mice<1e-10, arr.ind=TRUE)] <- 0
### replace zeros
require(zCompositions)
mice.no0 <- cmultRepl(mice, output = "p-counts")

### normalize data
mice.pro <- mice.no0 / rowSums(mice.no0)
mice.pro <- as.matrix(mice.pro)

### -------------------------------------------------------------------
### illustrating logratio analysis and total logratio variance
### several equivalent ways to obtain logratio variance are illustrated
### make sure 'easyCODA' package is installed
### skip this section if you just want to find best ALR transformation

require(easyCODA)
## first way: directly from sum of squared singular values of logratio analysis 
mice.lra <- LRA(mice.pro, weight=FALSE)
sum(mice.lra$sv^2)
# [1] 0.2099017

## second way: computed from the variances of the centered logratios
## note the averaging (/3147) as well as the adjustment (n-1)/n to get variances as average sum of squares
mice.clr <- CLR(mice.pro, weight=FALSE)
sum(apply(mice.clr$LR, 2, var))*(27/28)/3147
# [1] 0.2099017

## third way: using easyCODA's LR.VAR function
LR.VAR(mice.clr)
# [1] 0.2099017

## fourth way: same as third way but on transposed matrix
## note the renormalization
### on transposed matrix
tmice.no0 <- t(mice.no0)
tmice.pro <- tmice.no0 / rowSums(tmice.no0)
LR.VAR(CLR(tmice.pro, weight=FALSE))
# [1] 0.2099017

### -------------------------------------------------------------------

### use function FINDALR to identify best references
### note that FINDALR has default weight=FALSE (unweighted)

starttime <- Sys.time()
alr.refs <- FINDALR(mice.pro)
endtime   <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 82.66236 secs

alr.refs$procrust.max # [1] 0.9977376
alr.refs$procrust.ref # 1] 1318

alr.refs$var.min # [1] 0.004145442
alr.refs$var.ref # 1] 1557

starttime <- Sys.time()
rabbits.alr.refs <- FINDALR(Rabbits.pro)
endtime   <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 82.66236 secs

### quartiles of transcript #1557
quantile(log(mice.pro[,1557]), c(0,0.25,0.5,0.75,1))
#        0%       25%       50%       75%      100% 
# -8.322311 -8.217493 -8.177022 -8.138795 -8.034504 

### top 10 for correlation
mice.alr.results[order(mice.alr.results[,3], decreasing=TRUE),3][1:10]
#         Vcp    Pafah1b1   Arf1;Arf3       Rplp0        Aco2        Canx 
#   0.9977376   0.9974278   0.9970546   0.9968020   0.9965415   0.9964751 
# Vamp3;Vamp2        Wdr1 Cyfip1;Shyc       Rpl19 
#   0.9961483   0.9960804   0.9958606   0.9957104 

### top 10 for log variance
mice.alr.results[order(mice.alr.results[,4]),4][1:10]
#          Bag1      Pafah1b1        Ube2l3 Rab11a;Rab11b   Vamp3;Vamp2 
#   0.004145442   0.006259947   0.006903511   0.007030669   0.007209450 
#        Tmsb4x          Gnb1         Arpc3   Med15;Pcqap         Anxa2 
#   0.007711698   0.008519982   0.008785561   0.008980469   0.009315313 

### reference "Pafah1b1" chosen, which comes second on both criteria

### number of the chosen reference "Pafah1b1"
which(colnames(mice.pro)=="Pafah1b1")
# [1] 1179

### rank of the chosen reference "Pafah1b1" in decreasing order of average relative abundance
which(sort(colMeans(mice.pro), decreasing=TRUE)==sort(colMeans(mice.pro), decreasing=TRUE)["Pafah1b1"])
# Pafah1b1 
#     1617

### quartiles of transcript #1179
quantile(log(mice.pro[,1179]), c(0,0.25,0.5,0.75,1))
#        0%       25%       50%       75%      100% 
# -9.685366 -9.620854 -9.569166 -9.500670 -9.371171 

### -------------------------------------------------------------------
### logratio distances
mice.dist <- dist(mice.clr$LR)
alr1179.dist <- dist(ALR(mice.pro, denom=1179, weight=FALSE)$LR)
## FIGURE 4 of article
par(mar=c(4.2,4,1,1), mgp=c(2,0.7,0), las=1)
plot(mice.dist, alr1179.dist, xlab="Exact logratio distances",
     ylab="Distances using ALRs", font.lab=2)

### -------------------------------------------------------------------
### LRA of all logratios and PCA of ALRs w.r.t. #1179

## %s of variance explained by dimensions
round(100*mice.lra$sv^2/sum(mice.lra$sv^2),2)
#  [1] 59.24 15.54  7.62  3.97  1.93  1.85  1.23  1.15  0.98  0.75  0.68  0.62
# [13]  0.60  0.54  0.51  0.45  0.42  0.35  0.33  0.29  0.27  0.22  0.18  0.15
# [25]  0.14  0.00  0.00

## simplified sample names
mice.names <- paste(rep(c("M","M","L","L"),7), rep(c(0,1,2,4,6,9,12), each=4))

## FIGURE 5 of article
par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2, mfrow=c(2,1))
## LRA
plot(mice.lra$rowpcoord, type="n", asp=1, main="LRA of transcripts", 
     xlab="LRA dimension 1 (59.2%)", ylab="LRA dimension 2 (15.5%)")
abline(v=0, h=0, col="gray", lty=2)
text(mice.lra$rowpcoord, labels=mice.names, cex=0.8)

##PCA
alr1179 <- ALR(mice.pro, denom=1179, weight=FALSE)$LR
alr1179.pca <- PCA(alr1179, weight=FALSE)

## %s of explained variance
round(100*alr1179.pca$sv^2/sum(alr1179.pca$sv^2),2)
#  [1] 61.96 14.07  6.92  3.68  2.01  1.75  1.39  1.07  1.02  0.74  0.66  0.57
# [13]  0.55  0.49  0.46  0.43  0.39  0.36  0.30  0.27  0.25  0.21  0.18  0.14
# [25]  0.13  0.00  0.00  0.00

alr1179.rpc <- alr1179.pca$rowpcoord
plot(alr1179.pca$rowpcoord, type="n", asp=1, main="PCA of ALRs w.r.t. 1179", 
     xlab="PCA dimension 1 (62.0%)", ylab="PCA dimension 2 (14.1%)")
abline(v=0, h=0, col="gray", lty=2)
text(alr1179.pca$rowpcoord, labels=mice.names, cex=0.8)
### ----------------------------------------------------------------------------
### Finding ALRs for different sized data sets, from 500 to 3000 in steps of 500
procrustes.sim <- matrix(0, 100, 7)
colnames(procrustes.sim) <- c(100, seq(500, 3000, 500))
procrustes.ref <- matrix(0, 100, 7)
colnames(procrustes.ref) <- c(100, seq(500, 3000, 500))

set.seed(123)
for(i in 1:100) {
  for(j in 1:7) {
    ncomp <- c(100, seq(500, 3000, 500))[j]
    simul <- sample(1:ncol(mice.no0))[1:ncomp]
    mice.sub <- mice.pro[,simul]
    mice.sub.pro <- mice.sub / rowSums(mice.sub)
    results <- FINDALR(mice.sub.pro)
    procrustes.sim[i,j] <- results$procrust.max
    procrustes.ref[i,j] <- simul[results$procrust.ref]
  }
}

boxplot(as.numeric(procrustes.sim) ~ factor(rep(1:7, each=100)), ylim=c(0.987,0.999), bty="n", xaxt="n", yaxt="n",
        xlab="Number of components", ylab="Procrustes correlations", pars=list(boxwex=0.5))
axis(1, at=1:7, labels=c(100, seq(500, 3000, 500)), cex.axis=0.9)
axis(2, at=seq(0.987,0.999, 0.002), labels=seq(0.987,0.999, 0.002), cex.axis=0.8)

procrustes.ref
