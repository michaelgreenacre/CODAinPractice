### -------------------------------------------------------------------------------
### Install easyCODA package from CRAN in usual way or latest version from R-Forge:
### (Note: version 0.32 of easyCODA was used here)

install.packages("easyCODA", repos="http://R-Forge.R-project.org")

### -------------------------------------------
### Aar Massif data set is in data object 'Aar'
### in the package 'compositions'

library(compositions)
data(Aar)
### columns 3-12 are the oxide data
aar <- Aar[,3:12]
### close the data set
aar <- aar / rowSums(aar)
rownames(aar) <- 1:nrow(aar)

### ---------------------------------------------------------------
### compute the (unweighted) logratio variance three different ways
library(easyCODA)
### (1) as the sum of all pairwise logratio variances divided by 100 (=J squared)
### option weight=FALSE implies equal weights, each logratio gets weight (1/10)x(1/10)
### notice that the weights are always passed on in the logratio object
aar.LR <- LR(aar, weight=FALSE)
LR.VAR(aar.LR)
# [1] 0.1522912
### (2) as the sum of all CLR variances divided by 10 (=J)
aar.CLR <- CLR(aar, weight=FALSE)
LR.VAR(aar.CLR)
# [1] 0.1522912
### (3) as the average of all squared elements of the double-centred log-transformed data matrix
aarlog <- log(aar)
aarlog.DC <- sweep(aarlog, 1, rowMeans(aarlog))
aarlog.DC <- sweep(aarlog.DC, 2, colMeans(aarlog.DC))
sum(aarlog.DC^2) /(nrow(aar)*ncol(aar))
# [1] 0.1522912

### ---------------------------------------------------------------
### Stepwise selection of ratios, including the three amalgamations

### Define the amalgamations and add them to the set of 10 parts
mafic      <- rowSums(aar[,c(5,10,4)])
felsic     <- rowSums(aar[,c(7,1,3,8)])
carb_apat  <- rowSums(aar[,c(6,9)])
aar.amalg  <- cbind(aar, mafic, felsic, carb_apat)

### Perform the stepwise analysis
aar.step   <- STEP(aar.amalg, aar, weight=FALSE)

### Table 1 of article
cbind(aar.step$ratios, round(100*aar.step$R2max,1), round(aar.step$pro.cor,3))

### -------------------------------------------------------------------------
### LRA of original data and PCA of the 9 selected LRs
### (abbreviation LR = pairwise logratio, in this case including SLR balance)

par(mar=c(4.2,4,3,3), font.lab=2, cex.lab=1.2, mgp=c(2.7,0.7,0), las=1, mfrow=c(1,2), cex.axis=0.8)

### LRA (logratio analysis, = PCA of the CLRs) (Fig 2(a) of article)
aar.lra <- LRA(aar, weight=FALSE)
PLOT.LRA(aar.lra, map="contribution")

### PCA of the selected logratios (Fig 2(b) of article)
rownames(aar.step$logratios) <- 1:nrow(aar)
aar.ratios.pca <- PCA(aar.step$logratios, weight=FALSE)
PLOT.PCA(aar.ratios.pca, map="contribution", axes.inv=c(1,-1), rescale=2)

### --------------------------------------------------------
### Procrustes correlation of full-space geometry of samples
### (the element $rowpcoord means row principal coordinates)
protest(aar.ratios.pca$rowpcoord, aar.lra$rowpcoord)$t0
# [1] 0.993197

### Procrustes correlation of two configurations of samples in two dimensions
protest(aar.ratios.pca$rowpcoord[,1:2], aar.lra$rowpcoord[,1:2])$t0
# [1] 0.9971076

### ----------------------------------------------------------------------------
### Ratio selection after SLR balance of mafic/felsic is forced in at first step 
### First check the top 10 ratios at the first step
foo.step <- STEP(aar.amalg, aar, weight=FALSE, nsteps=1, top=10)
cbind(foo.step$ratios.top, round(100*foo.step$R2.top,1))
#               row col     
# MgO/Na2O        5   7 69.1
# Na2O/mafic      7  11 69.0
# MnO/felsic      4  12 68.9
# mafic/felsic   11  12 68.8
# Al2O3/mafic     3  11 68.6
# .....
### (mafic/felsic is 4th in the list, with an explained variance of 68.8%)

### The remaining ratios 
aar.step2 <- STEP(aar.amalg, aar, weight=FALSE, previous=log(aar.amalg[,"mafic"]/aar.amalg[,"felsic"]), previous.wt=0.3*0.4)
cbind(aar.step2$ratios, round(100*aar.step2$R2max,1))

### PCA of the new selection of logratios (Fig 3 of article)
### First add the mafic/felsic logratio at the start with the others
aar.ratios2 <- cbind(log(aar.amalg[,"mafic"]/aar.amalg[,"felsic"]), aar.step2$logratios)
colnames(aar.ratios2)[1] <- "mafic/felsic"
rownames(aar.ratios2) <- 1:nrow(aar)
aar.ratios2.pca <- PCA(aar.ratios2, weight=FALSE)
par(mar=c(4.2,4,3,3), font.lab=2, cex.lab=1.2, mgp=c(2.7,0.7,0), las=1, mfrow=c(1,1), cex.axis=0.8)
PLOT.PCA(aar.ratios2.pca, map="contribution", axes.inv=c(1,-1), rescale=2)

### ----------------------------------------------------------
### Procrustes correlation of full-space geometry of samples
### (for the new set of ratios, which is sub-"optimal" because
###  of forcing mafic/felsic in at step 1)
protest(aar.ratios2.pca$rowpcoord, aar.lra$rowpcoord)$t0
# [1] 0.9824111

### Procrustes correlation of two configurations of samples in two dimensions
protest(aar.ratios2.pca$rowpcoord[,1:2], aar.lra$rowpcoord[,1:2])$t0
# [1] 0.9817307
