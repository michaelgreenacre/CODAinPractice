### read data, 40 fatty acids (FAs) in 42 copepods sampled 
### in spring, summer and winter
FA0 <- read.csv("https://raw.githubusercontent.com/michaelgreenacre/CODAinPractice/master/copepods_TL.csv", header=TRUE, row.names=1, check.names=FALSE)

### season and total lipids (TL) in columns 1 and 2
FA_season <- FA0[,1]
table(FA_season)
# spring summer winter 
#     12     22      8 
FA_TL <- FA0[,2]
summary(FA_TL)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9425  1.3700  2.8300  3.7591  5.4363  8.8650 


### remove 'season' and 'TL' columns (just do this once!)
FA0 <- FA0[,-c(1,2)]

sum(FA0==0)
# [1] 187

### proportion of zeros
sum(FA0==0)/(nrow(FA0)*ncol(FA0))
# [1] 0.1113095

### for first analyses, use all FAs with mean over 1%
FA <- FA0[,colMeans(FA0)>1]
dim(FA)
# [1] 42 13

### close (normalize) to a subcomposition
FA <- FA / rowSums(FA)

### graph of ratios: all ratios, additive w.r.t. 14:0
### need to install package 'igraph'
require(igraph)
### all pairwise ratios
ratios <- matrix(0, 0.5*ncol(FA)*(ncol(FA)-1), 2)
ii <- 1
for(j in 2:ncol(FA)) {
  for(i in 1:(j-1)) {
    ratios[ii,] <- c(i,j)
    ii <- ii+1
  }
}
ratios.g <- graph(t(ratios), directed=FALSE)
V(ratios.g)$name <- colnames(FA)
plot(ratios.g, vertex.color="gray90", vertex.size=33, vertex.label.color="black", 
     vertex.label.cex=0.8, vertex.label.font=4, vertex.label.family="sans-serif",
     layout=layout_in_circle)

### additive ratios w.r.t. 14:0
ratios <- matrix(0, ncol(FA)-1, 2)
for(j in 2:ncol(FA)) {
  ratios[j-1,] <- c(1,j)
}
ratios.g <- graph(t(ratios), directed=TRUE)
V(ratios.g)$name <- colnames(FA)
plot(ratios.g, vertex.color="gray90", vertex.size=33, vertex.label.color="black", 
     vertex.label.cex=0.8, vertex.label.font=4, vertex.label.family="sans-serif",
     layout=layout_in_circle, edge.arrow.size=0.5)

### many ways of computing total logratio variance (TotVar)...

### data set in FA
I <- nrow(FA)
J <- ncol(FA)

### ... using variances of all pairwise logratios (LRs)
FA.LR <- LRtemp(FA, weight=FALSE)$LR
sum(apply(FA.LR, 2, VAR)) / J^2
# [1] 0.2461614

### ... using variances of centered logratios (CLRs)
FA.CLR <- CLR(FA, weight=FALSE)$LR
sum(apply(FA.CLR, 2, VAR)) / J
# [1] 0.2461614

### ordered contributions
round(sort(100*apply(FA.CLR, 2, VAR)/sum(apply(FA.CLR, 2, VAR))),1)
#  22:1(n-9) 22:1(n-11)  20:1(n-9)       16:0  18:1(n-9)       14:0  22:6(n-3)  20:1(n-7) 
#        1.1        1.3        1.9        2.0        2.0        2.1        2.6        2.6 
# 20:5(n-3)  18:2(n-6)  16:1(n-7)       18:0  18:4(n-3) 
#       3.5        3.6       13.2       18.4       45.7

### ... using the variation matrix
tau <- matrix(0, J, J)
for(j in 2:J) {
  for(k in 1:(j-1)) {
    tau[j,k] <- tau[k,j] <- VAR(log(FA[,j]/FA[,k]))
  }
}
### averaging the variation matrix gives twice the total variance
sum(tau)/ J^2
# [1] 0.4923229

### ... using the sample (between-row) logratio distances
FA.rowdist <- dist(FA.CLR)/sqrt(J)
### TotVar = sum of squared sample logratio distances divided by I^2
sum(FA.rowdist^2) / I^2
# [1] 0.2461614

### ... using the part (between-column) logratio distances
### first have to compute CLRs on transposed data
tFA <- t(FA) / rowSums(t(FA))
tFA.CLR <- CLR(tFA, weight=FALSE)$LR
FA.coldist <- dist(tFA.CLR)/sqrt(I)
### TotVar = sum of squared part logratio distances divided by J^2
sum(FA.coldist^2) / J^2
# [1] 0.2461614

### ... using double-centered log(FA), divided by sqrt(I*J)
dcFA <- sweep(log(FA), 1, rowMeans(log(FA)))
dcFA <- sweep(dcFA, 2, colMeans(dcFA))
sum(dcFA^2) / (I*J)
# [1] 0.2461614

### ... using logratio distances (between-row or between-column)
### obtained from the double-centered matrix
FA.rowdist <- dist(dcFA / sqrt(J))
FA.coldist <- dist(t(dcFA) / sqrt(I))
sum(FA.rowdist^2) / I^2
# [1] 0.2461614
sum(FA.coldist^2) / J^2
# [1] 0.2461614

### the variances of LRs are equal to squared part (between-column) distances
LR.VARs <- apply(FA.LR, 2, VAR)
head(LR.VARs,20)
head(FA.coldist^2,20)

### logratio analysis (LRA) = PCA of the LRs = PCA of the CLRs
FA.LRA <- LRA(FA, weight=FALSE)
### sum of singular values = TotVar
sum(FA.LRA$sv^2)
# [1] 0.2461614
round(FA.LRA$sv^2, 4)
# [1] 0.1629 0.0465 0.0181 0.0073 0.0041 0.0032 0.0014 0.0011 0.0008 0.0003 0.0002 0.0002

### rpc=row principal coordinates
### csc= column standard coordinates
### cpc=column principal coordinates
### rsc= row standard coordinates
FA.rpc <- FA.LRA$rowpcoord
FA.csc <- FA.LRA$colcoord
FA.cpc <- FA.LRA$colpcoord
FA.rsc <- FA.LRA$rowcoord
season.col <- c("chocolate","red","blue")
season.num <- as.numeric(FA_season)

### biplot and cluster analysis of samples
par(mar=c(4.2,4,2,1), mgp=c(1.8,0.6,0), font.lab=2, cex.lab=0.8, mfrow=c(1,2), cex.axis=0.7)
rescale <- 0.5
### biplot
plot(rbind(FA.rpc, rescale*FA.csc), type="n", asp=1, 
     xlab="LRA dimension 1 (66.2%)", ylab="LRA dimension 2 (18.9%)", las=1)
abline(v=0, h=0, col="gray", lty=2)
text(rescale*FA.csc, labels=colnames(FA), col="black", font=4, cex=0.8)
points(FA.rpc, pch=19, col=season.col[season.num], cex=0.8)
legend("bottomright", legend=c("spring","summer","winter"), 
       pch=19, col=season.col, text.col=season.col,
       title="Seasons", title.col="black", bty="n", cex=0.8)
### clustering, using WARD in easyCODA
FA.clust <- WARD(FA.CLR, weight=FALSE)
plot(FA.clust, hang=-1, labels=rep("", nrow(FA)), ylab="Cluster height", main="", xlab="Samples", 
     ylim=c(0, 0.45), las=1)
sample.order <- FA.clust$order
points((1:42), rep(-0.0002,42), pch=19, col=season.col[season.num][sample.order], cex=0.8)

### sum of squared cluster node heights = TotVar
sum(FA.clust$height^2)
# [1] 0.2461614

### clustering of samples using hclust in R
FA.clust2 <- hclust(dist(FA.rpc), method="ward.D2")
sum(FA.clust2$height^2)
# [1] 20.67756
### sum of squared heights = 2I times TotVar
sum(FA.clust2$height^2) / (2*I)
# [1] 0.2461614

### biplot and cluster analysis of the parts
par(mar=c(4.2,4,1,1), mgp=c(1.8,0.6,0), font.lab=2, cex.lab=0.8, mfrow=c(1,2), cex.axis=0.7, las=1)
rescale <- 0.3
plot(rbind(rescale*FA.rsc, FA.cpc), type="n", asp=1, main="",
     xlab="LRA dimension 1 (66.2%)", ylab="LRA dimension 2 (18.9%)")
abline(v=0, h=0, col="gray", lty=2)
points(rescale*FA.rsc, pch=19, col=season.col[season.num], cex=0.8)
text(FA.cpc, labels=colnames(FA), col="black", font=4, cex=0.8)
legend("bottomright", legend=c("spring","summer","winter"), 
       pch=19, col=season.col, text.col=season.col,
       title="Seasons", title.col="black", bty="n", cex=0.8)
### clustering the parts uding WARD
FA.clust3 <- WARD(tFA.CLR, weight=FALSE)
plot(FA.clust3, hang=-1, main="", labels=colnames(FA), font=4, cex=0.8, font.lab=4, yaxt="n",
      xlab="Parts (fatty acids)",  ylab="Cluster height", ylim=c(0,0.4), yaxt="n")
axis(2, at=seq(0,0.3,0.1), labels=seq(0,0.3,0.1))

### sum of squared cluster node heights = TotVar
sum(FA.clust$height^2)
# [1] 0.2461614

### decomposition of TotVar along nodes or along dimensions
rev(round(100*FA.clust$height^2/sum(FA.clust$height^2),2))
#  [1] 59.57 15.81  4.73  3.13  2.59  2.13  1.71  1.38 ...
round(100*FA.LRA$sv^2 / sum(FA.LRA$sv^2), 2)
#  [1] 66.19 18.91  7.34  2.98  1.67  1.28  0.58  0.44 ...

### ---------------------------------------------------------------
### LR selection

require(easyCODA)

### notice we use unweighted LRs here

FA.STEP <- STEP(FA, weight=FALSE)
FA.STEP$names
#  [1] "16:0/18:4(n-3)"      "16:1(n-7)/18:0"      "16:1(n-7)/20:1(n-9)"
#  [4] "14:0/22:6(n-3)"      "14:0/20:1(n-9)"      "18:2(n-6)/20:5(n-3)"
#  [7] "18:0/18:2(n-6)"      "20:1(n-9)/20:1(n-7)" "18:1(n-9)/20:1(n-7)"
# [10] "18:4(n-3)/22:1(n-9)" "16:0/18:1(n-9)"      "18:0/22:1(n-11)" 
FA.STEP$R2max
# [1] 0.6555195 0.8376898 0.9089140 0.9413417 0.9592644 0.9760015 0.9840515 0.9898978
# [9] 0.9951828 0.9975647 0.9989557 1.0000000

### first 3 explain more than 90% of the variance
cbind(FA.STEP$names[1:5], round(100*FA.STEP$R2max[1:5],2))
# [1,] "16:0/18:4(n-3)"      "65.55"
# [2,] "16:1(n-7)/18:0"      "83.77"
# [3,] "16:1(n-7)/20:1(n-9)" "90.89"
# [4,] "14:0/22:6(n-3)"      "94.13"
# [5,] "14:0/20:1(n-9)"      "95.93"

### use the 5 parts from the first 3 LRs as a subcomposition
FAsub5 <- cbind( FA[,"16:0"], FA[,"16:1(n-7)"], FA[,"18:0"], 
                FA["18:4(n-3)"], FA[,"20:1(n-9)"])
colnames(FAsub5) <- c("16:0","16:1(n-7)","18:0","18:4(n-3)","20:1(n-9)")
FAsub5 <- FAsub5 / rowSums(FAsub5)

### the CLRs of these 5 parts explain slightly more % of variance 
### because the 3 logratios are not connected
FAsub5.CLR <- CLR(FAsub5, weight=FALSE)$LR
FAsub5.LR  <- LR(FAsub5, weight=FALSE)$LR

rda(FA.CLR ~ FA.STEP$logratios[,1:3])
#               Inertia Proportion Rank
# Total         3.27815    1.00000     
# Constrained   2.97956    0.90891    3
# Unconstrained 0.29859    0.09109    9

rda(FA.CLR ~ FAsub5.LR)
#               Inertia Proportion Rank
# Total         3.27815    1.00000     
# Constrained   3.03961    0.92723    4
# Unconstrained 0.23854    0.07277    8

rda(FA.CLR ~ FAsub5.CLR)
#               Inertia Proportion Rank
# Total         3.27815    1.00000     
# Constrained   3.03961    0.92723    4
# Unconstrained 0.23854    0.07277    8


### PCA of just the three selected logratios
FA3.PCA <- PCA(FA.STEP$logratios[,1:3], weight=FALSE)
summary(FA3.PCA)
#   dim    value      %   cum%   scree plot               
#   1      0.751469  62.7  62.7  ****************         
#   2      0.409952  34.2  96.9  *********                
#   3      0.036857   3.1 100.0  *                        
#   -------- -----                                 
#   Total: 1.198278 100.0 

FA3.rpc <- FA3.PCA$rowpcoord
FA3.csc <- FA3.PCA$colcoord
### reverse second axis
FA3.rpc[,2] <- -FA3.rpc[,2]
FA3.csc[,2] <- -FA3.csc[,2]
par(mar=c(4.2,4,1,1), mgp=c(1.8,0.6,0), font.lab=2, cex.lab=0.8, mfrow=c(1,2), cex.axis=0.7)
rescale <- 0.5
plot(rbind(FA3.rpc, rescale*FA3.csc), type="n", asp=1, 
     xlab="PCA dimension 1 (62.7%)", ylab="PCA dimension 2 (34.2%)")
abline(v=0, h=0, col="gray", lty=2)
arrows(0, 0, 0.93*rescale*FA3.csc[,1], 0.93*rescale*FA3.csc[,2], col="gray", lwd=2, length=0.1, angle=10)
text(rescale*FA3.csc, labels=FA.STEP$names[1:3], col="black", font=4, cex=0.8)
points(FA3.rpc, pch=19, col=season.col[season.num], cex=0.8)
legend("bottomleft", legend=c("spring","summer","winter"), 
       pch=19, col=season.col, text.col=season.col,
       title="Seasons", title.col="black", bty="n", cex=0.8)

FA3.clust <- WARD(FA.STEP$logratios[,1:3], weight=FALSE)
sample3.order <- FA3.clust$order
plot(FA3.clust, hang=-1, labels=rep("", nrow(FA)),  main="", ylab="Cluster height",
     xlab="Samples", yaxt="n")
axis(2, at=seq(0,0.8,0.1), labels=seq(0,0.8,0.1))
points((1:42), rep(-0.0002,42), pch=19, col=season.col[season.num][sample3.order], cex=0.8)

### two subcompositions: five-part one above, and the other 8-part one
index.sub5 <- c(2,3,4,7,8)
index.sub8 <- (1:13)[-index.sub5]
FAsub5 <- FA[, index.sub5]
FAsub5 <- FAsub5 / rowSums(FAsub5)
FAsub8 <- FA[, index.sub8]
FAsub8 <- FAsub8 / rowSums(FAsub8)

### variance explained: 5-part and 8-part subcompositions
rda(FA.CLR ~ CLR(FAsub5, weight=FALSE)$LR)
#               Inertia Proportion Rank
# Total         3.27815    1.00000     
# Constrained   3.03961    0.92723    4   92.7% explained
# Unconstrained 0.23854    0.07277    8

rda(FA.CLR ~ CLR(FAsub8, weight=FALSE)$LR)
#               Inertia Proportion Rank
# Total          3.2781     1.0000     
# Constrained    2.8232     0.8612    7   86.1% explained
# Unconstrained  0.4549     0.1388    5

### conditioning on one of the subcompositions at a time
rda(FA.CLR ~ CLR(FAsub8, weight=FALSE)$LR + Condition(CLR(FAsub5, weight=FALSE)$LR))
#                Inertia Proportion Rank
# Total         3.278150   1.000000     
# Conditional   3.039610   0.927233    4
# Constrained   0.229247   0.069932    7
# Unconstrained 0.009293   0.002835    1

rda(FA.CLR ~ CLR(FAsub5, weight=FALSE)$LR + Condition(CLR(FAsub8, weight=FALSE)$LR))
#                Inertia Proportion Rank
# Total         3.278150   1.000000     
# Conditional   2.823202   0.861218    7
# Constrained   0.445655   0.135947    4
# Unconstrained 0.009293   0.002835    1

### including both subcompositions (thus only one cross logratio missing)
rda(FA.CLR ~ cbind(CLR(FAsub5, weight=FALSE)$LR, CLR(FAsub8, weight=FALSE)$LR))
#                Inertia Proportion Rank
# Total         3.278150   1.000000     
# Constrained   3.268857   0.997165   11
# Unconstrained 0.009293   0.002835    1

### variances contained
### "full" 13-part composition  (13*12/2 = 78 LRs)
LR.VAR(FA.CLR, weight=FALSE)
# [1] 0.2461614
sum(apply(FA.LR, 2, VAR))/13^2
# [1] 0.2461614

### 5-part subcomposition (5*4/2 = 10 LRs)
sum(apply(FAsub5.LR, 2, VAR))/13^2
# [1] 0.07564679

### 8-part subcomposition (8*7/2 = 28 LRs)
sum(apply(FAsub8.LR, 2, VAR))/13^2
# [1] 0.0272609

### remaining logratios (40 LRs)
sum(apply(FA.LR, 2, VAR))/13^2 - sum(apply(FAsub5.LR, 2, VAR))/13^2 - sum(apply(FAsub8.LR, 2, VAR))/13^2
# [1] 0.1432537

### decomposition 0.24616 = 0.07565 + 0.02726 + 0.14325
###  #LR pairs       78   =    10   +    28    +   40


### ---------------------------------------------------------------
### modeling

### case 1: compositional data are responses
FA.rda <- rda(FA.CLR ~ factor(FA_season))
summary(FA.rda)
# Partitioning of variance:
#               Inertia Proportion
# Total          3.2781     1.0000
# Constrained    2.4996     0.7625  (76.3% explained)
# Unconstrained  0.7786     0.2375

### notice that rda's "inertia" has to be averaged by J 
### and adjusted by (I-1)/I to get TotVar
(FA.rda$tot.chi / J) * (I-1)/I
# [1] 0.2461614
### the part of explained variance
(FA.rda$CCA$tot.chi / J) * (I-1)/I
# [1] 0.1876975

### case 2: 5-part subcomposition are responses
summary(rda(FAsub5.CLR ~ factor(FA_season)))
# Partitioning of variance:
#               Inertia Proportion
# Total          2.6192     1.0000
# Constrained    2.2067     0.8425  (84.3% explained)
# Unconstrained  0.4125     0.1575

### case 3: 8-part subcomposition are responses
summary(rda(FAsub8.CLR ~ factor(FA_season)))
# Partitioning of variance:
#               Inertia Proportion
# Total          0.5899     1.0000
# Constrained    0.2754     0.4669  (46.7% explained)
# Unconstrained  0.3145     0.5331

### Total Lipids (TL) as response

### do a stepwise search for explanatory variables of log(FA_TL)
R2 <- rep(0, 13*12/2)
ratios <- matrix("", 13*12/2, 2)
jk <- 1
for(k in 2:13) {
  for(j in 1:(k-1)) {
    foo.lm <- lm(log(FA_TL) ~ log(FA[,j]/FA[,k]))
    R2[jk] <- summary(foo.lm)$r.squared
    ratios[jk,] <- colnames(FA)[c(j,k)]
    jk <- jk+1
  }
}
max(R2)
# [1] 0.7423907

### order the FA ratios
ratios.order <- order(R2, decreasing=TRUE)
cbind(ratios[ratios.order,],round(R2[ratios.order],3))[1:10,]
#  [1,] "18:0"      "20:1(n-9)"  "0.742"
#  [2,] "18:0"      "18:4(n-3)"  "0.726"
#  [3,] "18:0"      "22:1(n-11)" "0.71"
#  [4,] "18:0"      "22:1(n-9)"  "0.706"
#  [5,] "18:4(n-3)" "20:1(n-7)"  "0.702"
#  [6,] "18:0"      "18:2(n-6)"  "0.664"
#  [7,] "18:0"      "18:1(n-9)"  "0.662"
#  [8,] "18:2(n-6)" "18:4(n-3)"  "0.651"
#  [9,] "16:0"      "18:4(n-3)"  "0.65" 
# [10,] "18:4(n-3)" "22:6(n-3)"  "0.648"

### second step 
R2 <- rep(0, 13*12/2)
ratios <- matrix("", 13*12/2, 2)
jk <- 1
for(k in 2:13) {
  for(j in 1:(k-1)) {
    foo.lm <- lm(log(FA_TL) ~ log(FA[,"18:0"]/FA[,"20:1(n-9)"]) + log(FA[,j]/FA[,k]))
    R2[jk] <- summary(foo.lm)$r.squared
    ratios[jk,] <- colnames(FA)[c(j,k)]
    jk <- jk+1
  }
}
max(R2)
# [1] 0.8326866

### order the FA ratios
ratios.order <- order(R2, decreasing=TRUE)
cbind(ratios[ratios.order,],round(R2[ratios.order],3))[1:10,]
# [1,] "18:4(n-3)" "20:1(n-7)"  "0.833"
# [2,] "18:4(n-3)" "22:6(n-3)"  "0.819"
# [3,] "16:0"      "18:4(n-3)"  "0.816"
# [4,] "18:4(n-3)" "22:1(n-11)" "0.806"
# [5,] "14:0"      "18:4(n-3)"  "0.805"
# [6,] "18:1(n-9)" "18:4(n-3)"  "0.804"
# [7,] "18:2(n-6)" "18:4(n-3)"  "0.804"
# [8,] "18:4(n-3)" "22:1(n-9)"  "0.801"
# [9,] "18:4(n-3)" "20:5(n-3)"  "0.801"
#[10,] "18:0"      "18:4(n-3)"  "0.799"  

### the best two LRs are (8,4) and (7,9)
### estimate coefficients using logs to base 2 
summary(lm(log(FA_TL,2) ~ log(FA[,8]/FA[,4],2)+log(FA[,7]/FA[,9],2)))
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -0.34947    0.18034  -1.938   0.0599 .  
# log(FA[, 8]/FA[, 4], 2)  0.49537    0.08980   5.517 2.43e-06 ***
# log(FA[, 7]/FA[, 9], 2)  0.32302    0.07041   4.588 4.56e-05 ***
# Multiple R-squared:  0.8327,    Adjusted R-squared:  0.8241 
# F-statistic: 97.05 on 2 and 39 DF,  p-value: 7.225e-16


### multiplicative effects of the LR predictors
2^0.49537; 2^0.32302
# [1] 1.409682
# [1] 1.250946


### adding a link 8/7 to make a four-part subcomposition
summary(lm(log(FA_TL,2) ~ log(FA[,8]/FA[,4], 2)+log(FA[,7]/FA[,9] ,2) + log(FA[,8]/FA[,7] ,2)))
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -1.19162    0.48410  -2.462 0.018487 *  
# log(FA[, 8]/FA[, 4], 2)  0.41450    0.09727   4.261 0.000129 ***
# log(FA[, 7]/FA[, 9], 2)  0.65633    0.19126   3.432 0.001461 ** 
# log(FA[, 8]/FA[, 7], 2)  0.24065    0.12900   1.866 0.069833 .  
# Multiple R-squared:  0.8467,    Adjusted R-squared:  0.8346 
# F-statistic: 69.97 on 3 and 38 DF,  p-value: 1.549e-15

### SLR suggested...
summary(lm(log(FA_TL,2) ~ log((FA[,7]+FA[,8])/(FA[,4]+FA[,9]), 2)))
# Coefficients:
#                                                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                     -0.70365    0.20508  -3.431  0.00141 ** 
# log((FA[, 7] + FA[, 8])/(FA[, 4] + FA[, 9]), 2)  0.90892    0.07702  11.802 1.32e-14 ***
# Multiple R-squared:  0.7769,    Adjusted R-squared:  0.7713 
# F-statistic: 139.3 on 1 and 40 DF,  p-value: 1.321e-14


### ILR suggested...
summary(lm(log(FA_TL, 2) ~ log(sqrt((FA[,7]*FA[,8]))/sqrt((FA[,4]*FA[,9])), 2)))
# Coefficients:
#                                                             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                 -0.22455    0.14478  -1.551    0.129    
# log(sqrt((FA[, 7] * FA[, 8]))/sqrt((FA[, 4] * FA[, 9])), 2)  0.79442    0.05745  13.827   <2e-16 ***
# Multiple R-squared:  0.827,     Adjusted R-squared:  0.8227 
# F-statistic: 191.2 on 1 and 40 DF,  p-value: < 2.2e-16

### LR choice with expert knowledge
### Martin chose 18:4(n-3)/20:1(n-7) and 22:1(n-9)/16:0  (numbers 7/9 & 12/2) 
summary(lm(log(FA_TL, 2) ~ log(FA[,7]/FA[,9],2)+log(FA[,12]/FA[,2],2)))
# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               3.50509    0.70318   4.985 1.31e-05 ***
# log(FA[, 7]/FA[, 9], 2)   0.48767    0.05865   8.314 3.62e-10 ***
# log(FA[, 12]/FA[, 2], 2)  0.97723    0.22081   4.426 7.52e-05 ***
# Multiple R-squared:  0.8017,    Adjusted R-squared:  0.7915 
# F-statistic: 78.84 on 2 and 39 DF,  p-value: 1.983e-14


### ----------------------------------------------------------------
### the coefficients of the log-contrast of c(4,7,8,9)
FAsub4 <- lm(log(FA_TL,2) ~ log(FA[,8]/FA[,4], 2)+log(FA[,7]/FA[,9] ,2) + log(FA[,8]/FA[,7] ,2))$coefficients[2:4]
FAsub4.a <- rep(0,4)
FAsub4.a[1] <- -FAsub4[1]
FAsub4.a[2] <- FAsub4[2]-FAsub4[3]
FAsub4.a[3] <- FAsub4[1]+FAsub4[3]
FAsub4.a[4] <- -FAsub4[2]
sum(FAsub4.a) # checks = 0

##########################################################################################
### Analysis of full composition FA0
apply(FA0==0, 2, sum)
#     14:0  14:1(n-5)     i-15:0     a-15:0       15:0  15:1(n-6)     i-16:0       16:0 
#        0          6          0          4          0          8         40          0 
# 16:1(n-9)  16:1(n-7)  16:1(n-5)     i-17:0     a-17:0  16:2(n-4)       17:0  16:3(n-4) 
#       39          0          0          4          5          0         15          0 
# 16:4(n-1)       18:0  18:1(n-9)  18:1(n-7)  18:2(n-6)  18:3(n-6)  18:3(n-3)  18:4(n-3) 
#         7          0          0          0          0         11          0          0 
#      20:0 20:1(n-11)  20:1(n-9)  20:1(n-7)  20:2(n-6)  20:3(n-6)  20:4(n-6)  20:3(n-3) 
#         0         30          0          0          1          0          5          8 
# 20:4(n-3)  20:5(n-3) 22:1(n-11)  22:1(n-9)  22:1(n-7)  22:5(n-3)  22:6(n-3)  24:1(n-9) 
#         0          0          0          0          0          0          0          4 

### smallest positive values
FAmin <- min(FA0[FA0[,1]>0,1])
for(j in 2:ncol(FA0)) FAmin <- c(FAmin, min(FA0[FA0[,j]>0,j]))

#--------
### (2/3)min method, 
FA1.min <- FA0
for(j in 1:ncol(FA0)) FA1.min[FA0[,j]==0,j] <- (2/3) * FAmin[j]
FA1.min <- FA1.min/rowSums(FA1.min)
### (FAs 7 and 9 have to be omitted, because some methods broke down with so many zeros)
sum(LRA(FA1.min[,-c(7,9)], weight=FALSE)$sv^2)
# [1] 0.3506684

#--------
require(zCompositions)
### lrDA: Imputation by single simulated values
FA1.lrDA <- lrDA(FA0[,-c(7,9)], label=0, ini.cov="multRepl", dl=FAmin[-c(7,9)], n.iters=150)
FA1.lrDA <- FA1.lrDA/rowSums(FA1.lrDA)
sum(LRA(FA1.lrDA, weight=FALSE)$sv^2)
# [1] 0.7987203

protest(LRA(FA1.lrDA, weight=FALSE)$rowpcoord, LRA(FA1.min[,-c(7,9)], weight=FALSE)$rowpcoord)$t0
# [1] 0.8253652

#--------
require(robCompositions)
### lrEM method
FA1.lrEM <- lrEM(FA0[,-c(7,9)], label=0, ini.cov="multRepl", dl=FAmin[-c(7,9)])
FA1.lrEM <- FA1.lrEM/rowSums(FA1.lrEM)
sum(LRA(FA1.lrEM, weight=FALSE)$sv^2)
# [1] 0.4550464

protest(LRA(FA1.lrEM, weight=FALSE)$rowpcoord, LRA(FA1.min[,-c(7,9)], weight=FALSE)$rowpcoord)$t0
# [1] 0.9409831
protest(LRA(FA1.lrEM, weight=FALSE)$rowpcoord, LRA(FA1.lrDA, weight=FALSE)$rowpcoord)$t0
# [1] 0.8512292

#--------------------
### BDLs model-based
FA1.BDLs <- imputeBDLs(FA0[,-c(7,9)], method="lmrob")
FA1.BDLs$x <- FA1.BDLs$x/rowSums(FAimp.BDLs$x)
sum(LRA(FA1.BDLs$x, weight=FALSE)$sv^2)
# [1] 0.818633

protest(LRA(FA1.BDLs$x, weight=FALSE)$rowpcoord, LRA(FA1.min[,-c(7,9)], weight=FALSE)$rowpcoord)$t0
# [1] 0.8342288
protest(LRA(FA1.BDLs$x, weight=FALSE)$rowpcoord, LRA(FA1.lrDA, weight=FALSE)$rowpcoord)$t0
# [1] 0.7984514
protest(LRA(FAimp.BDLs$x, weight=FALSE)$rowpcoord, LRA(FA1.lrEM, weight=FALSE)$rowpcoord)$t0
# [1] 0.8285052

### accumulate all 108 substituted values in FA.zeros
sum(FA0[-c(7,9)]==0)
# 108

FA.zeros <- matrix(0, 108, 4)
colnames(FA.zeros) <- c("2/3","lrDA","lrEM","BDLs")
zero.locs <- which(FA0[-c(7,9)]==0, arr.ind=TRUE)
for(i in 1:108) FA.zeros[i,1] <- 100*FA1.23[zero.locs[i,1], zero.locs[i,2]]
for(i in 1:108) FA.zeros[i,2] <- 100*FA1.lrDA[zero.locs[i,1], zero.locs[i,2]]    
for(i in 1:108) FA.zeros[i,3] <- 100*FA1.lrEM[zero.locs[i,1], zero.locs[i,2]]    
for(i in 1:108) FA.zeros[i,4] <- 100*FAimp.BDLs$x[zero.locs[i,1], zero.locs[i,2]]    
apply(FA.zeros, 2, quantile, c(0.25,0.5,0.75))
#            2/3        lrDA       lrEM         BDLs
# 25% 0.03534661 0.007181469 0.01757676 0.0001253048
# 50% 0.04695294 0.021474773 0.04247181 0.0302875937
# 75% 0.06191269 0.059849171 0.06257556 0.0333122465

### many small substituted values for BDLs:
sum(FA.zeros[,4]<0.001)
# [1] 29

### plot the zeros
require(colorspace)
col.zeros <- rainbow_hcl(5)

### put substitutes into categories for bar-plotting
zero.cats <- (seq(0, max(FA.zeros), length=51))
zero.tab  <- matrix(0, nrow=50, ncol=4)
for(j in 1:4) {
  for(itab in 1:50) {
    zero.tab[itab,j] <- sum((FA.zeros[,j]>zero.cats[itab]) & (FA.zeros[,j]<zero.cats[itab+1]))
  }
}
lows <- 160 - cumsum(apply(zero.tab, 2, max)+10) - 35
# PDF file: width=6, height=4
par(mar=c(4.2,4,1,1), mgp=c(1.5,0,0), font.lab=2, cex.lab=1, mfrow=c(1,1), cex.axis=0.8)
plot(0,0,type="n", xlab="Substituted values (%)", ylab="Counts for each method         ", xlim=c(-7,51), ylim=c(0,120),
     bty="n", xaxt="n", yaxt="n")
axis(1, at=c(1,16,32,48), labels=c("0", "0.05", "0.10", "0.15"), tick=FALSE)
for(j in 1:4) segments(0.6, lows[j], 50, lows[j], col="gray")
for(j in 1:4) {
  low <- lows[j]
  for(itab in 1:50) segments(itab, low, itab, low+zero.tab[itab,j], col=col.zeros[j], 
                          lwd=6, lend=1)
}
text(rep(0.5,4), lows+2.5, labels=rev(c("BDLs","lrEM","lrDA","(2/3)min")), 
     col=col.zeros, font=2, pos=2)

###########################################
# microbiome data set (Baxter etal, 2016) #
###########################################
# --------------------------------------------------------------------------------------------
#### read the Baxter microbiome data
baxter <- read.table("https://raw.githubusercontent.com/michaelgreenacre/CODAinPractice/master/Baxter_OTU_table.txt", header=TRUE)
dim(baxter)
# [1] 490 338  (490 samples, 335 bacteria)

#### remove first 3 columns (only do this once!)
baxter <- baxter[,-c(1:3)]

### percentage of data zeros
100 * sum(baxter==0) / (nrow(baxter)*ncol(baxter))
# [1] 58.70911

### -----------------------------------------------------------------------------------------
### read meta data of Baxter microbiome study
meta <- read.table("https://raw.githubusercontent.com/michaelgreenacre/CODAinPractice/master/Baxter_Metadata.txt", header=TRUE)
dim(meta)
# [1] 490  27
colnames(meta)
#  [1] "sample"       "fit_result"   "Site"         "Dx_Bin"       "dx"           "Hx_Prev"      "Hx_of_Polyps"
#  [8] "Age"          "Gender"       "Smoke"        "Diabetic"     "Hx_Fam_CRC"   "Height"       "Weight"      
# [15] "BMI"          "White"        "Native"       "Black"        "Pacific"      "Asian"        "Other"       
# [22] "Ethnic"       "NSAID"        "Abx"          "Diabetes_Med" "stage"        "Location" 

### the group labels, also convert to numbers
dx <- meta[,"dx"]
table(dx)
# adenoma  cancer  normal 
#     198     120     172 
dx.num <- as.numeric(as.factor(dx))
table(dx.num)
#   1   2   3 
# 198 120 172 
sex <- meta[,"Gender"]
table(sex)
#   f   m 
# 243 247 
sex.num <- as.numeric(as.factor(sex))
table(sex.num)
#   1   2 
# 243 247 

### CCA of the Baxter data, with groups and age as constraining variables
### using fourth-root transformation
baxter.pro <-100 * baxter/rowSums(baxter)
baxter.pro <- baxter.pro^0.25
baxter.pro <- baxter.pro / rowSums(baxter.pro)
require(easyCODA)
dx.Z <- DUMMY(dx.num)
colnames(dx.Z) <- c("A","C","N")
### there's an outlier, sample 371
baxter.cm <- colMeans(baxter.pro)
baxter.cca <- cca(baxter.pro ~ dx.Z + scale(meta[,"Age"]))

baxter.cm <- colMeans(baxter.pro[-371,])
baxter.cca <- cca(baxter.pro[-371,] ~ dx.Z[-371,] + scale(meta[-371,"Age"]))
#               Inertia Proportion Rank
# Total         1.62263    1.00000     
# Constrained   0.02174    0.01340    3
# Unconstrained 1.60089    0.98660  334
# (1.34% of the variance explained by group differences)
anova(baxter.cca, by="terms")
#                       Df ChiSquare      F Pr(>F)    
# dx.Z                   2   0.01434 2.1065  0.001 ***
# scale(meta[, "Age"])   1   0.00808 2.3735  0.001 ***
# Residual             486   1.65414   
anova(baxter.cca, by="axis")
#           Df ChiSquare      F Pr(>F)    
# CCA1       1   0.01081 3.1752  0.001 ***
# CCA2       1   0.00842 2.4727  0.001 ***
# CCA3       1   0.00319 0.9386  0.613    
# Residual 486   1.65414  

### plot the results in the constrained space of the group centroids
### and only show the high contributing OTUs
baxter.rpc <- baxter.cca$CCA$wa %*% diag(1/sqrt(apply(baxter.cca$CCA$wa^2, 2, sum))) %*% diag(sqrt(baxter.cca$CCA$eig))
baxter.ccc <- baxter.cca$CCA$v %*% diag(1/sqrt(apply(baxter.cca$CCA$v^2, 2, sum)))
baxter.cpc <- (baxter.cca$CCA$v / sqrt(490)) %*% diag(sqrt(baxter.cca$CCA$eig))

apply( baxter.cca$CCA$v^2 * baxter.cm, 2, sum)
# CCA1 CCA2 CCA3 
#    1    1    1 

## reverse the coordinates
baxter.rpc <- -baxter.rpc
baxter.ccc <- -baxter.ccc
baxter.crd <- (baxter.ccc[,1]^2 * baxter.cca$CCA$eig[1] + baxter.ccc[,2]^2 * baxter.cca$CCA$eig[2]) / (baxter.cca$CCA$eig[1]+baxter.cca$CCA$eig[2])
baxter.ctr <- baxter.crd > mean(baxter.crd)  # note: mean(baxter.crd) = 1/length(baxter.crd)
## how many have contributions higher than average
sum(baxter.ctr)
# [1] 72
dx.aggr <- aggregate(baxter.rpc ~ as.factor(dx.num[-371]), FUN=mean)
colnames(dx.aggr) <- c("Dim","CCA1","CCA2","CCA3")
dx.aggr
#   Dim         CCA1          CCA2          CCA3
# 1   1 -0.001388531  0.0004392753 -1.099845e-03
# 2   2  0.004851585 -0.0001914329  8.809265e-05
# 3   3 -0.001794475 -0.0003695656  1.198246e-03

summary(lm(meta[-371,"Age"] ~ baxter.rpc[,1]+baxter.rpc[,2]))
age.coef <- lm(scale(meta[-371,"Age"]) ~ scale(baxter.rpc[,1])+scale(baxter.rpc[,2])+scale(baxter.rpc[,3]))$coefficients[2:4]
age.coef/30

## plot "by hand"
require(ellipse)
require(colorspace)
## needs my function CIplot_biv which can be downloaded here:
# source("https://www.fbbva.es/microsite/multivariate-statistics/static/CIplot_biv.R")
dx.col <- c("forestgreen","chocolate", "blue")
OTUshort <- substr(colnames(baxter.pro),7,9)
rescale <- 0.1
# pdf(file="Figure9new.pdf", width=8, family="ArialMT", height=6, useDingbats=FALSE)
par(mar=c(4,4,1,1), mgp=c(2.5,0.8,0), font.lab=2, cex.axis=0.8)
plot(rbind(rescale*baxter.ccc), type="n", asp=1, xlab="CCA dimension 1",
     ylab="CCA dimension 2") #, ylim=c(-0.013,0.02))
abline(v=0, h=0, col="gray", lty=2)
arrows(0,0,0.95*rescale*baxter.ccc[baxter.ctr,1],0.95*rescale*baxter.ccc[baxter.ctr,2], length=0.1, angle=10, col="pink", lwd=1)
arrows(0,0,0.95*age.coef[1]/30,0.95*age.coef[2]/30, length=0.1, angle=10, col="gray", lwd=2)
points(baxter.rpc, pch=1, col=dx.col[dx.num[-371]], cex=0.6)
text(rescale*baxter.ccc[baxter.ctr,], labels=OTUshort[baxter.ctr], col="red", cex=0.7, font=4)
text(age.coef[1]/30,age.coef[2]/30, label="Age", font=4, cex=0.8)
legend("topright", legend=c("Adenoma","Cancer","Normal"), pch=1, bty="n", col=dx.col, 
       text.col=dx.col, pt.cex=0.8)
for(i in 1:5) {
  set.seed(123)
  CIplot_biv(baxter.rpc[,1], baxter.rpc[,2], group=dx.num[-371], shownames=FALSE, groupcol=rep("white",6), shade=TRUE, add=TRUE, alpha=0.999) 
}
set.seed(123)
CIplot_biv(baxter.rpc[,1], baxter.rpc[,2], group=dx.num[-371], shownames=FALSE, groupcol=dx.col, shade=TRUE, add=TRUE, alpha=0.99) 
set.seed(123)
CIplot_biv(baxter.rpc[,1], baxter.rpc[,2], group=dx.num[-371], groupnames=c("A","C","N"), groupcol=dx.col,  add=TRUE, alpha=0.99) 

# dev.off()

## top 10 contributors to the first dimension
OTUshort[order(baxter.ccc[,1]^2, decreasing=TRUE)][1:10]
[1] "260" "310" "105" "281" "264" "297" "057" "288" "298" "340"


#-----------------------------------------------
### Subcompositional incoherence stress exercise
### Two functions: STRESS and chidist
STRESS <- function(D, Dhat) {
# Function to compute Kruskal stress 1 between 
# two distance matrices
# D and Dhat are distance objects
  foo.num <- 0
  foo.den <- 0
  D <- as.numeric(D)
  Dhat <- as.numeric(Dhat)
  for(ij in i:length(D)) {
    foo.num <- foo.num+(D[ij]-Dhat[ij])^2
    foo.den <- foo.den + D[ij]^2
  }
  sqrt(foo.num/foo.den)
}

chidist <- function(mat,rowcol=1) {
# function to calculate chi-square distances between row or column
# profiles of a matrix
# e.g. chidist(N,1) calculates the chi-square distances between row profiles
#      (for row profiles, chidist(N) is sufficient)
#      chidist(N,2) calculates the chi-square distances between column profiles
  mat <- as.matrix(mat)
  if(rowcol==1) {
    prof<-mat/apply(mat,1,sum)
    rootaveprof<-sqrt(apply(mat,2,sum)/sum(mat))
  }
  if(rowcol==2) {
    prof<-t(mat)/apply(mat,2,sum)
    rootaveprof<-sqrt(apply(mat,1,sum)/sum(mat))
  }
  dist(scale(prof,F,rootaveprof))  
}


### exclude outlier
baxter.pro1 <- baxter[-371,] / rowSums(baxter[-371,])

incoh.CA   <- matrix(0, nrow=100, ncol=9)
baxter.cm1 <- colMeans(baxter.pro1)
D.chi      <- as.matrix(chidist(baxter.pro1, 2))
set.seed(123)
j <- 1
nparts <- 10*j

set.seed(123)
for(j in seq(9,1-1)) {
  nparts <- 10*j
  for(i in 1:100) {
# find the subcompositional parts  
    jparts <- sample(1:335, ceiling(335*nparts/100), prob=baxter.cm1)
    foo <- baxter.pro1[,jparts]
# remove parts all zeros
    allzero <- which(colSums(foo)==0)
    if(length(allzero)>0) {
      jparts <- jparts[-allzero]
      foo <- baxter.pro[,jparts]
    }
# incoherence in CA  
    D <- as.dist(D.chi[jparts, jparts])
# remove samples all zeros
    allzero <- which(rowSums(foo)==0)
    if(length(allzero)>0) {
      foo <- foo[-allzero,]
    }
    D2 <- chidist(foo, 2)
    incoh.CA[i,j] <- STRESS(D,D2)
  }
}

incoh.CA.quants <- apply(incoh.CA, 2, quantile, c(0.025,0.975))
round(incoh.CA.quants,4)
#         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]
#  2.5% 0.0300 0.0274 0.0229 0.0205 0.0166 0.0128 0.0092 0.0050 0.0015
# 97.5% 0.1485 0.0990 0.0950 0.0792 0.0516 0.0382 0.0311 0.0202 0.0132


### Figure 10
plot(rep(1:9, each=2), as.numeric(incoh.CA.quants), xlab="Percentage in subcomposition",
     ylab="Stress", bty="n", xaxt="n", ylim=c(0, 0.20), type="n", font.lab=2)
for(j in 1:9) segments(j, incoh.CA.quants[1,j], j, incoh.CA.quants[2,j], col="blue", lwd=2)
eps <- 0.06
for(j in 1:9) segments(j-eps, incoh.CA.quants[1,j], j+eps, incoh.CA.quants[1,j], col="blue", lwd=2, lend=2)
for(j in 1:9) segments(j-eps, incoh.CA.quants[2,j], j+eps, incoh.CA.quants[2,j], col="blue", lwd=2, lend=2)
axis(1, at=1:9, labels=c("10%","20%","30%","40%","50%","60%","70%","80%","90%"))
points(1:9, apply(incoh.CA, 2, median), pch=21, col="blue", bg="white", cex=0.9)


### top 10 OTUs, numerical order, not OTU number 50% subcompositions
top10 <- order(baxter.ccc[,1]^2, decreasing=TRUE)[1:10]
# [1] 239 279 105 256 243 269  57 263 270 301
OTUshort[top10]
# [1] "260" "310" "105" "281" "264" "297" "057" "288" "298" "340"
cancer <- dx[-371]=="cancer"
table(cancer)
# FALSE  TRUE 
#   369   120 
cancer.preds <- as.matrix(baxter.pro1)^0.5
cancer.preds <- cancer.preds / rowSums(cancer.preds)
cancer.preds <- cancer.preds %*% diag(1/apply(cancer.preds, 2, max)) 
colnames(cancer.preds) <- colnames(baxter)
cancer.preds <- as.data.frame(cancer.preds)
foo.glm <- glm(factor(cancer) ~ Otu000310+Otu000105, family="binomial", data=as.data.frame(cancer.preds))
summary(foo.glm)
Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -1.671      0.134 -12.466  < 2e-16 ***
# Otu000310      6.918      1.851   3.738 0.000186 ***
# Otu000105     19.157      4.786   4.002 6.27e-05 ***
  
cancer.regr <- matrix(0, 100, 2)
cancer.regr[1,] <- foo.glm$coefficients[2:3]
cancer.se <- summary(foo.glm)$coefficients[2:3,2]
cancer.ci <- cbind(cancer.regr[1,]-2*cancer.se, cancer.regr[1,]+2*cancer.se)
cancer.pvals <- matrix(0, 100, 2)
cancer.pvals[1,] <- summary(foo.glm)$coefficients[2:3,4]
cancer.predict <- matrix(0,100,4)
foo.predict <- predict(foo.glm)
table(foo.predict>0, cancer)
#           cancer
#         FALSE TRUE
#   FALSE   364   83
#   TRUE      5   37

cancer.predict[1,] <- as.numeric(table(foo.predict>0, cancer))
[1] 357  13  84  36

### check 50% subcompositions including the top10
set.seed(123)
for(iter in 2:100) {
  baxter.sub <- baxter.pro1[,top10]
  colnames(baxter.sub) <- colnames(baxter)[top10]
  baxter.sub <- cbind(baxter.sub, baxter.pro1[, sample((1:335)[-top10], 158, prob=baxter.cm1[-top10])])
  baxter.sub.pro <- baxter.sub/rowSums(baxter.sub)
  cancer.preds2 <- as.matrix(baxter.sub.pro)^0.5
  cancer.preds2 <- cancer.preds2 / rowSums(cancer.preds2)
  cancer.preds2 <- cancer.preds2 %*% diag(1/apply(cancer.preds2, 2, max)) 
  colnames(cancer.preds2) <- colnames(baxter.sub)
  cancer.preds2 <- as.data.frame(cancer.preds2)
  foo.glm <- glm(factor(cancer) ~ Otu000310+Otu000105, family="binomial", data=cancer.preds2)
  cancer.regr[iter,] <- foo.glm$coefficients[2:3]
  cancer.pvals[iter,] <- summary(foo.glm)$coefficients[2:3,4]
  foo.predict <- predict(foo.glm)
  cancer.predict[iter,] <- as.numeric(table(foo.predict>0, cancer))
}

apply(cancer.pvals, 2, range)

# Figure window: width=8, height=5
par(mar=c(2.2,4,2,1), mgp=c(2.5,0.7,0), font.lab=2, mfrow=c(1,2), xpd=TRUE)
plot(0,0,type="n",xlim=c(0.5, 3),ylim=c(0,30),
     xlab="", ylab="Logistic regression coefficients", bty="n", xaxt="n")
for(i in 1:100) points(rep(1,100), cancer.regr[,1], pch=21, col="blue", cex=0.8)
for(i in 1:100) points(rep(2,100), cancer.regr[,2], pch=21, col="red", cex=0.8)
eps <- 0.1
segments(1, cancer.ci[1,1], 1, cancer.ci[1,2], col="blue", lend=2 )
segments(2, cancer.ci[2,1], 2, cancer.ci[2,2], col="red", lend=2 )
segments(1-0.5*eps, cancer.ci[1,1], 1+0.5*eps, cancer.ci[1,1], col="blue", lend=2 )
segments(2-0.5*eps, cancer.ci[2,1], 2+0.5*eps, cancer.ci[2,1], col="red", lend=2 )
segments(1-0.5*eps, cancer.ci[1,2], 1+0.5*eps, cancer.ci[1,2], col="blue", lend=2 )
segments(2-0.5*eps, cancer.ci[2,2], 2+0.5*eps, cancer.ci[2,2], col="red", lend=2 )
segments(1-eps, cancer.regr[1,1], 1+eps, cancer.regr[1,1], col="blue", lwd=2, lend=2 )
segments(2-eps, cancer.regr[1,2], 2+eps, cancer.regr[1,2], col="red", lwd=2, lend=2 )
text(1:2, cancer.ci[,1]-1, labels=c("OTU310","OTU105"),
     col=c("blue", "red"), font=4, cex=0.9)
text(0, 30, label="(a)", pos=3, offset=1.5, cex=1.2)

plot(0.1,0.1,type="n",xlim=c(0.5, 3), ylim=c(10^-5, 10^-3),log="y",
     xlab="", ylab="p-values", bty="n", xaxt="n", yaxt="n", las=1)
axis(2, at=c(10^-5, 10^-4, 10^-3), labels=c("10^-5","10^-4","10^-3"))
for(i in 1:100) points(rep(1,100), cancer.pvals[,1], pch=21, col="blue", cex=0.8)
for(i in 1:100) points(rep(2,100), cancer.pvals[,2], pch=21, col="red", cex=0.8)
eps <- 0.1
segments(1-eps, cancer.pvals[1,1], 1+eps, cancer.pvals[1,1], col="blue", lwd=2, lend=2 )
segments(2-eps, cancer.pvals[1,2], 2+eps, cancer.pvals[1,2], col="red", lwd=2, lend=2 )
text(1:2, apply(cancer.pvals, 2, min), labels=c("OTU310","OTU105"),
     col=c("blue", "red"), font=4, cex=0.9, pos=1)
text(0, 10^-3, label="(b)", pos=3, offset=1.5, cex=1.2)


