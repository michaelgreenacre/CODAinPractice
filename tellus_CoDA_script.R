### -----------------------------------------------------------------------
### "Aitchison's Compositional Data Analysis 40 years On: A Reappraisal"
### This is the analysis of Tellus cation data with zeros

### Set your working directory where the data have been donwloaded

### Software note:
### package easyCODA required
library(easyCODA)
### For the FINDALR function, either make sure you have the latest version of easyCODA
### (presently version 0.35.1 on RForge)
### or if you have installed from CRAN include the function directly from GitHub as follows:
### source("https://raw.githubusercontent.com/michaelgreenacre/CODAinPractice/master/FINDALR.R")

TELLUS <- read.table("tellus.xrf.a.cation.txt", header=TRUE)
dim(TELLUS)
# [1] 6799   61
colnames(TELLUS)
#  [1] "Sample"     "Easting"    "Northing"   "Unitname"   "AgeBracket" "col"       
#  [7] "ab"         "pH"         "LOI"        "Si"         "Al"         "Fe"        
# [13] "Mg"         "Mn"         "Ca"         "Na"         "K"          "P"         
# [19] "Ti"         "S"          "Ag"         "As"         "Ba"         "Bi"        
# [25] "Br"         "Cd"         "Ce"         "Cl"         "Co"         "Cr"        
# [31] "Cs"         "Cu"         "Ga"         "Ge"         "Hf"         "I"         
# [37] "In"         "La"         "Mo"         "Nb"         "Nd"         "Ni"        
# [43] "Pb"         "Rb"         "Sb"         "Sc"         "Se"         "Sm"        
# [49] "Sn"         "Sr"         "Ta"         "Te"         "Th"         "Tl"        
# [55] "U"          "V"          "W"          "Y"          "Yb"         "Zn"        
# [61] "Zr"        
  
### Age Brackets (AB)    
table(TELLUS[,7])
# CzOl CzPl  Mes NeoP   Pg   Pl PlCr PlDv PlOr PlSi 
#  154 1691  330 1013  124  170 1534  304 1311  168 
AB <- as.numeric(factor(TELLUS[,7]))
table(AB)
#    1    2    3    4    5    6    7    8    9   10 
#  154 1691  330 1013  124  170 1534  304 1311  168 
ABnames <- unique(TELLUS[,7])
ABnames <- sort(ABnames)

### Tellus cation data, with zeros
tellus0 <- TELLUS[,10:61]
dim(tellus0)
# [1] 6799   52

### Number of zeros and percentage of total
sum(tellus0==0) 
# [1] 3883
100*sum(tellus0==0)/(nrow(tellus0)*ncol(tellus0))
# [1} 1.098295
colSums(tellus0==0)
#   Si   Al   Fe   Mg   Mn   Ca   Na    K    P   Ti    S   Ag   As   Ba   Bi   Br   Cd   Ce 
#    0    0    0    0    0    0    0    0    0    0 2535    0    1    0  480    0    0    0 
#   Cl   Co   Cr   Cs   Cu   Ga   Ge   Hf    I   In   La   Mo   Nb   Nd   Ni   Pb   Rb   Sb 
#    0    0    0    0    2   48    0    0    0    0    0   92    0  581    0    0    0    0 
#   Sc   Se   Sm   Sn   Sr   Ta   Te   Th   Tl    U    V    W    Y   Yb   Zn   Zr 
#   41    0    0    0    0   59    0    1    0   10    0   33    0    0    0    0 

### Number of samples with some zeros
length(which(rowSums(tellus0==0)>0)) 
# [1] 3303
table(rowSums(tellus0==0))
#    0    1    2    3    4    5 
# 3496 2851  332  113    6    1 

### minimum positive values observed
tellus.min <- min(tellus0[tellus0[,1]>0,1])
for(j in 2:52) tellus.min <- c(tellus.min, min(tellus0[tellus0[,j]>0,j]))

### replace by 2/3 minimum positive value to form matrix tellus
tellus <- tellus0
for(j in 1:52) {
  if(sum(tellus0[,j]==0) > 0) {
    for(i in 1:nrow(tellus)) tellus[tellus0[,j]==0,j] <- tellus.min[j]*2/3
  }
}
sum(tellus==0)
# [1] 0

### Data matrices are tellus: with replaced values; tellus0: with zeros
### .pro are the normalized/closed profiles
tellus.pro  <- tellus/rowSums(tellus)
tellus0.pro <- tellus0/rowSums(tellus0)

### Negative correlations in the original Tellus data and closed data
tellus0.cor <- cor(tellus0)
sum(as.dist(tellus0.cor)<0)/length(as.dist(tellus0.cor))
# [1] 0.387632
tellus0.pro.cor <- cor(tellus0.pro)
sum(as.dist(tellus0.pro.cor)<0)/length(as.dist(tellus0.pro.cor))
# [1] 0.3989442

### Table of positive and negative correlations in original and closed data
table(as.dist(tellus0.cor)>0, as.dist(tellus0.pro.cor)>0)
#         FALSE TRUE
#   FALSE   307  207
#   TRUE    222  590
 	
### average percentages of elements
round(100*colMeans(tellus.pro),5)
#       Si       Al       Fe       Mg       Mn       Ca       Na        K        P 
# 66.71287 17.22345  4.76249  3.06781  0.09769  1.83961  2.36898  2.30087  0.31112 
#       Ti        S       Ag       As       Ba       Bi       Br       Cd       Ce 
#  0.75834  0.35673  0.00003  0.00133  0.02041  0.00002  0.00855  0.00005  0.00223 
#       Cl       Co       Cr       Cs       Cu       Ga       Ge       Hf        I 
#  0.06613  0.00209  0.02024  0.00021  0.00530  0.00142  0.00014  0.00024  0.00080 
#       In       La       Mo       Nb       Nd       Ni       Pb       Rb       Sb 
#  0.00002  0.00118  0.00007  0.00104  0.00076  0.00640  0.00237  0.00460  0.00009 
#       Sc       Se       Sm       Sn       Sr       Ta       Te       Th       Tl 
#  0.00216  0.00012  0.00029  0.00021  0.00758  0.00004  0.00001  0.00018  0.00003 
#        U        V        W        Y       Yb       Zn       Zr 

### rare earth correlations vs. in full composition
rare <- c(29, 18, 32, 37, 39, 49, 50, 44, 46)
tellus0.rare.pro <- CLOSE(tellus0.pro[,rare])
round(cor(tellus0.rare.pro), 3)
#        La     Ce     Nd     Sc     Sm      Y     Yb     Th      U
# La  1.000  0.920  0.218 -0.910  0.328 -0.344  0.292  0.625  0.412
# Ce  0.920  1.000  0.176 -0.910  0.353 -0.368  0.319  0.618  0.398
# Nd  0.218  0.176  1.000 -0.396 -0.694  0.123 -0.722  0.328 -0.249
# Sc -0.910 -0.910 -0.396  1.000 -0.224  0.069 -0.186 -0.724 -0.464
# Sm  0.328  0.353 -0.694 -0.224  1.000 -0.312  0.990 -0.032  0.603
# Y  -0.344 -0.368  0.123  0.069 -0.312  1.000 -0.300  0.056 -0.028
# Yb  0.292  0.319 -0.722 -0.186  0.990 -0.300  1.000 -0.046  0.588
# Th  0.625  0.618  0.328 -0.724 -0.032  0.056 -0.046  1.000  0.311
# U   0.412  0.398 -0.249 -0.464  0.603 -0.028  0.588  0.311  1.000
round(cor(tellus0.pro[,rare]), 3)
#        La     Ce     Nd     Sc     Sm     Y     Yb     Th      U
# La  1.000  0.892  0.830 -0.208  0.175 0.479  0.117  0.612  0.332
# Ce  0.892  1.000  0.741 -0.167  0.173 0.403  0.129  0.551  0.251
# Nd  0.830  0.741  1.000  0.033 -0.207 0.627 -0.263  0.546  0.233
# Sc -0.208 -0.167  0.033  1.000 -0.074 0.313 -0.044 -0.388 -0.149
# Sm  0.175  0.173 -0.207 -0.074  1.000 0.032  0.985 -0.009  0.322
# Y   0.479  0.403  0.627  0.313  0.032 1.000  0.024  0.459  0.543
# Yb  0.117  0.129 -0.263 -0.044  0.985 0.024  1.000 -0.036  0.312
# Th  0.612  0.551  0.546 -0.388 -0.009 0.459 -0.036  1.000  0.429
# U   0.332  0.251  0.233 -0.149  0.322 0.543  0.312  0.429  1.000

### should we weight the parts?
tellus.clr.unw <- CLR(tellus.pro, weight=FALSE)$LR
tellus.clr.unw.var <- apply(tellus.clr.unw, 2, var)
tellus.pro.cm <- colMeans(tellus.pro)
### ---------------------------------------------------
### Figure 4: plotting CLR variances against part means
### (pdf and dev.off functions commented out, can be used for saving PDFs)
# pdf(file="Fig_4.pdf", width=5, height=5, useDingbats=FALSE, family="ArialMT")
par(mar=c(4.2,4,1,1), font.lab=2, las=1, mfrow=c(1,1))
plot(tellus.pro.cm, tellus.clr.unw.var, log="xy", type="n", 
     xlab="Average compositional values (log-scale)", ylab="Variance of CLR (log-scale)")
text(tellus.pro.cm, tellus.clr.unw.var, labels=colnames(tellus), col="red", font=4, cex=0.8)
# dev.off()
### ---------------------------------------------------
### seems like unweighted OK: no tendency for rarer elements to have much higher variance
### in fact Al, the second highest component, has the lowest CLR variance


### variances and their contributions to total
TotVar <- mean(tellus.clr.unw.var)
TotVar
# [1] 0.3446613

### sort the part contributions to variance
sort(100*tellus.clr.unw.var/sum(tellus.clr.unw.var), decreasing=TRUE)
#         Nd         Br          S         Cl         Mn         Ni          I         Cr 
# 10.9354250  5.5526814  5.5490654  5.3977073 	 5.1673774  4.2466589  3.8078010  3.4451737 
#         Rb         Cu         Co         Ga         Pb         Sc         Se         Cd 
#  3.4010389  2.9310840  2.8993020  2.6291876  2.1762499  2.0443790  2.0355740  2.0174516 
#         :          :          :          :          :          :          :          :
#         :          :          :          :          :          :          :          :        Sm         Nb         Ge         Al 
#         Sm         Nb         Ge         Al 
#  0.3686688  0.3586418  0.3582653  0.1494869 

### find the denominator part for an ALR transformation
(tellus.findalr <- FINDALR(tellus.pro))
# $tot.var
# [1] 0.3446107 

# $procrust.cor
#  [1] 0.9662964 0.9907692 0.9339259 0.9545244 0.8336251 0.8938565 0.8945049 0.9181235 0.9289943
# [10] 0.9633254 0.8580249 0.9636980 0.8731239 0.9593046 0.8155351 0.8670293 0.8849160 0.9726034
# [19] 0.8794430 0.8981615 0.8811368 0.9300113 0.8417604 0.8948470 0.9663585 0.9619381 0.8150313
# [28] 0.9465207 0.9696384 0.8760541 0.9826409 0.7984858 0.8619365 0.8771585 0.8727109 0.9585282
# [37] 0.8797559 0.9212887 0.9892701 0.9690770 0.8738061 0.8632258 0.9623934 0.9128071 0.9812157
# [46] 0.9062143 0.9344361 0.9084048 0.9631968 0.9889152 0.8965161 0.9197674

# $procrust.max
# [1] 0.9907692

# $procrust.ref
# [1] 2

# $var.log
#  [1] 0.01105921 0.03038162 0.33140143 0.17374619 0.95729658 0.39101702 0.20187174 0.21395846 0.14493670
# [10] 0.13661518 1.33979835 0.18759453 0.43182290 0.07191086 0.34766988 1.34610010 0.53016945 0.08702211
# [19] 1.24307113 0.60417549 0.63413672 0.24444501 0.63733602 0.41873633 0.12317329 0.05033178 0.97267331
# [28] 0.27404725 0.09099433 0.20412732 0.05161335 1.79752728 0.87174489 0.58182783 0.46392747 0.24809741
# [37] 0.41898994 0.59088562 0.12563614 0.18092795 0.24537927 0.27541449 0.16340938 0.26874333 0.16899488
# [46] 0.30902500 0.36712962 0.35735229 0.10523465 0.12860962 0.34451145 0.17111440

# $var.min
# [1] 0.01105921

# $var.ref
# [1] 1

### (notice that Si has lowest variance of log, but Al had lowest variance of CLR
### Al has highest Procrustes correlation and the second lowest variance of log

### 5-number summary of log(Al)
quantile(log(tellus.pro[,"Al"]), c(0.025, 0.25, 0.5, 0.75, 0.975))
#      2.5%       25%       50%       75%     97.5% 
# -2.179636 -1.863984 -1.758703 -1.658929 -1.480048 


### -------------------------------------------------
### Figure 5: diagnosis of the reference part for ALR
### (pdf and dev.off functions commented out, can be used for saving PDFs)
# pdf(file="Fig_4.pdf", width=5, height=5, useDingbats=FALSE, family="ArialMT")
par(mar=c(4.2,4,1,1), font.lab=2, las=1, mfrow=c(1,1))
plot(tellus.findalr$var.log, tellus.findalr$procrust.cor, ylim=c(0.8,1), 
     xlab="Variance of log", ylab="Procrustes correlation")
points(tellus.findalr$var.log[2], tellus.findalr$procrust.cor[2], col="red", pch=21, bg="red")
points(tellus.findalr$var.log[1], tellus.findalr$procrust.cor[1], col="blue", pch=21, bg="blue")
text(tellus.findalr$var.log[2], tellus.findalr$procrust.cor[2], col="red", label="Al", pos=3, font=4)
text(tellus.findalr$var.log[1], tellus.findalr$procrust.cor[1], col="blue", label="Si", pos=3, font=4)
# dev.off()
### -------------------------------------------------

### ordering in terms of ALR variances (ALRs w.r.t. Al)
tellus.alr <- ALR(tellus.pro, denom=2, weight=FALSE)$LR
tellus.order.alr <- order(apply(tellus.alr, 2, var), decreasing=TRUE)
tellus.log <- log(tellus.pro[,c(1,3:52)])

### -------------------------------------------------------
### Figure S1: all log-transforms versus ALRs w.r.t. ref Al
### (pdf and dev.off functions commented out, can be used for saving PDFs)
# pdf(file="Fig_S1.pdf", width=6, height=11, useDingbats=FALSE, family="ArialMT")
### use a very tall vertical window to fit in the 51 plots in a 9-by-6 grid
par(mar=c(1,0.5,2,0.5), mgp=c(2,0.7,0), cex.axis=0.8, mfrow=c(9,6))
for(j in 1:51) plot(tellus.alr[,tellus.order.alr[j]], tellus.log[,tellus.order.alr[j]], 
                    main=colnames(tellus)[-2][tellus.order.alr[j]], cex=0.4, 
                    ylab="", xlab="", xaxt="n",yaxt="n", col="lightblue")
# dev.off()
### -------------------------------------------------------

### comparing distances between two samples using all logratios and using the ALRs
tellus.clr.unw <- CLR(tellus.pro, weight=FALSE)$LR
tellus.alr     <- ALR(tellus.pro, denom=2, weight=FALSE)$LR
### using 10000 random distance pairs
foo <- matrix(0, nrow=10000, ncol=2)
k <- 1
set.seed(123)
sample1 <- sample(1:6799, 10000, replace=TRUE)
sample2 <- sample(1:6799, 10000, replace=TRUE)
for(i in 1:10000) {
  if(sample1[i]==sample2[i]) next
  foo[k,1] <- sqrt(sum((tellus.clr.unw[sample1[i],] - tellus.clr.unw[sample2[i],])^2)) / sqrt(52)
  foo[k,2] <- sqrt(sum((tellus.alr[sample1[i],] - tellus.alr[sample2[i],])^2)) / sqrt(51)  
  k <- k+1
}
sum(foo[,2]<foo[,1])
# [1] 0   (all distances based on ALRs below the corresponding ones based on CLRs)

### ---------------------------------------------------------
### Figure 6: Scatterplot of distances based on CLRs and ALRs
### (pdf and dev.off functions commented out, can be used for saving PDFs)
# pdf(file="Fig_6.pdf", width=5, height=5, useDingbats=FALSE, family="ArialMT")
par(mar=c(4.2,4,1,1), font.lab=2, las=1)
plot(foo[,1], foo[,2], xlab="Distance based on ALRs w.r.t. Al", ylab="Logratio distance based on all LRs", 
     xlim=c(0,2), ylim=c(0,2), asp=1, col="lightblue", cex=0.5)
abline(a=0, b=1, col="red", lty=2)
# dev.off()
### ---------------------------------------------------------

### total variance in logratio analysis (the CLRs, equivalently all LRs)
tellus.lra <- LRA(tellus.pro, weight=FALSE)
sum(tellus.lra$sv^2)  
# [1] 0.3446107
### total variance in PCA of ALRs (ref: Al)
tellus.pca <- PCA(tellus.alr, weight=FALSE)
### (remember that all total variances are averaged, not summed)
sum(tellus.pca$sv^2)  
# [1] 0.3786807

### percentages of variance for CLRs in LRA
round(100*tellus.lra$sv^2/sum(tellus.lra$sv^2),3)
#  [1] 45.188 23.813  4.934  3.131  2.555  2.088  1.847  1.790  1.352  1.148  1.105  1.057  0.950  0.841
# [15]  0.812  0.675  0.634  0.560  0.477  0.428  0.404  0.374  0.358  0.335  0.277  0.249  0.237  0.224
# [29]  0.204  0.195  0.171  0.167  0.157  0.146  0.141  0.132  0.119  0.110  0.105  0.088  0.083  0.073
# [43]  0.054  0.052  0.040  0.031  0.024  0.022  0.015  0.013  0.013

### percentages of variance for ALRs (ref: Al)
round(100*tellus.pca$sv^2/sum(tellus.pca$sv^2),3)
#  [1] 44.893 22.101  7.330  2.914  2.606  2.083  1.801  1.713  1.461  1.231  1.065  1.013  0.912  0.807
# [15]  0.778  0.708  0.601  0.539  0.519  0.435  0.375  0.362  0.334  0.316  0.263  0.250  0.231  0.212
# [29]  0.207  0.189  0.170  0.159  0.154  0.146  0.135  0.129  0.118  0.110  0.102  0.091  0.082  0.070
# [43]  0.065  0.049  0.048  0.031  0.025  0.022  0.018  0.013  0.012


### WARD clustering of parts (i.e., on transposed matrix) using 10% of the cases 
tellus.10 <- tellus.pro[seq(1,nrow(tellus), 10),]
dim(tellus.10)
# [1] 680  52
tellus.clus <- WARD(CLR(CLOSE(t(tellus.10)), weight=FALSE), weight=FALSE)
### -----------------------------------------------------------------------
### Figure 2: Ward clustering of parts
### (pdf and dev.off functions commented out, can be used for saving PDFs)
# pdf(file="Fig_2.pdf", width=10, height=5, useDingbats=FALSE, family="ArialMT")
par(mar=c(4.2,4,1,1), mgp=c(2,0.7,0), font.lab=2)
plot(tellus.clus, labels=colnames(tellus), xlab="Elements", ylab="Height", main="")
# dev.off()
### -----------------------------------------------------------------------

### Amalgamation clustering of parts (matrix not transposed) using 10% of the cases 
### (this stepwise algorithm takes some time --- needs optimizing)
tellus.aclus <- ACLUST(tellus.10, weight=FALSE)
### Alternative vertical scale 
tellus.aclus.alt <- tellus.aclus
tellus.aclus.alt$height <- 100*tellus.aclus$height/tellus.aclus$height[51]
### -----------------------------------------------------------------------
### Figure 3: Amalgamation clustering of parts
### (pdf and dev.off functions commented out, can be used for saving PDFs)
# pdf(file="Fig_3.pdf", width=10, height=5, useDingbats=FALSE, family="ArialMT")
par(mar=c(4.2,4,1,1), mgp=c(2,0.7,0), font.lab=2)
plot(tellus.aclus.alt, labels=colnames(tellus), xlab="Elements",
     ylab="Percentage variance (%)", main="")
# dev.off()
### -----------------------------------------------------------------------

### LRA of tellus and row principal coordinates 
### .rpc = row principal coordinates
tellus.lra <- LRA(tellus.pro, weight=FALSE)
tellus.lra.rpc <- tellus.lra$rowpcoord

### PCA of tellus ALRs w.r.t. Al and row principal coordinates
tellus.alr.al <- ALR(tellus.pro, denom=2, weight=FALSE)$LR
tellus.pca <- PCA(tellus.alr.al, weight=FALSE)
tellus.pca.rpc <- tellus.pca$rowpcoord

### Procrustes correlations between all-logratios and ALRs... 
### ...in full space 
protest(tellus.lra.rpc, tellus.pca.rpc, permutations=0)$t0
# [1] 0.9907692
### ... and in reduced 2-D space
protest(tellus.lra.rpc[,1:2], tellus.pca.rpc[,1:2], permutations=0)$t0
# [1] 0.9971027

### contribution coordinates (in LRA equal weights are 1/52)
tellus.lra.ccc <- tellus.lra$colcoord * sqrt(1/52)
### high contributors
tellus.lra.ctr <- (tellus.lra.ccc[,1]^2 > 1/ncol(tellus)) | (tellus.lra.ccc[,2]^2 > 1/ncol(tellus)) 
sum(tellus.lra.ctr)
[1] 25

### colours for Age Bracket groups
require(colorspace)
tellus.col <- rainbow_hcl(10, l=50, c=70)

### function add.alpha for colour transparency
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

tellus.col.alpha <- add.alpha(rainbow_hcl(10, l=50, c=70), 0.2)
col <- c("blue","red")  # colours for possible use in graphics
tellus.pch <- c(5,3,1,2,4,5,3,1,2,4)

### -----------------------------------------------------------------------
### Figure 7: left and right figures of the LRA biplot, and Age Brackets
### save as png insert into PPT, and eventually save all as png

# invert 2nd axis for this figure (only do this once)
tellus.lra.rpc[,2] <- -tellus.lra.rpc[,2]
tellus.lra.ccc[,2] <- -tellus.lra.ccc[,2]

### use horizontal rectangular window
rescale <- 2 # for points
dim <- c(1,2)
perc.hor <- 45.2; perc.ver <- 23.8
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * rbind(tellus.lra.rpc, rescale*tellus.lra.ccc), type = "n", asp = 1, 
     xlab = paste("LRA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("LRA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = 
        round(axTicks(3)/rescale, 2), col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = 
        round(axTicks(4)/rescale, 2), col = "black", col.ticks = col[2], col.axis = col[2])
arrows(0, 0, 0.92 * rescale * tellus.lra.ccc[tellus.lra.ctr, 1], 
       0.92 * rescale * tellus.lra.ccc[tellus.lra.ctr, 2], length = 0.1, angle = 10, col = "pink")
points(tellus.lra.rpc, pch = tellus.pch[AB], col = tellus.col.alpha[AB], font = 1, cex = 0.5)
text(rescale * tellus.lra.ccc[tellus.lra.ctr,], labels = colnames(tellus.pro)[tellus.lra.ctr], 
     col = "red", cex = 0.9, font = 4)    
legend("bottomleft", legend=ABnames, pch=tellus.pch, 
       col=tellus.col, text.col=tellus.col, pt.cex=0.6, cex=0.8)
# dev.off()

require(ellipse)

# png(file="Fig_7_right.png",width=7,height=5.5,units="in",res=144)
rescale <- 3  # for ellipses
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * 0.44*tellus.lra.rpc, type = "n", asp = 1, 
     xlab = paste("LRA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("LRA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)

set.seed(123)
CIplot_biv(tellus.lra.rpc[,1], tellus.lra.rpc[,2], group=AB, groupcols=tellus.col, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           shownames=FALSE)
set.seed(123)
CIplot_biv(tellus.lra.rpc[,1], tellus.lra.rpc[,2], group=AB, groupcols=tellus.col, 
           add=TRUE, shade=FALSE, groupnames=ABnames, alpha=0.99)

# dev.off()

### same for dimension reduction of ALRs (both axes reversed here, 
### weights here = 1/51 for 51 ALRs
tellus.pca.rpc <- -tellus.pca$rowpcoord
tellus.pca.ccc <- -tellus.pca$colcoord * sqrt(1/51)
tellus.pca.ctr <- (tellus.pca.ccc[,1]^2 > 1/ncol(tellus)) | (tellus.pca.ccc[,2]^2 > 1/ncol(tellus)) 
sum(tellus.pca.ctr)
[1] 28

### -----------------------------------------------------------------------
### Figure 8: left and right figures of the ALR biplot, and Age Brackets
# save as png insert into PPT, and eventually save all as png
# png(file="Fig_8_left.png",width=7,height=5.5,units="in",res=144)
rescale <- 2 # for points
dim <- c(1,2)
perc.hor <- 44.9; perc.ver <- 22.1
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * rbind(tellus.pca.rpc, rescale*tellus.pca.ccc), type = "n", asp = 1, 
     xlab = paste("PCA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("PCA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
arrows(0, 0, 0.92 * rescale * tellus.pca.ccc[tellus.pca.ctr, 1], 
       0.92 * rescale * tellus.pca.ccc[tellus.pca.ctr, 2], 
       length = 0.1, angle = 10, col = "pink")
points(tellus.pca.rpc, pch = tellus.pch[AB], col = tellus.col.alpha[AB], font = 1, cex = 0.5)
text(rescale * tellus.pca.ccc[tellus.pca.ctr,], labels = colnames(tellus.alr)[tellus.pca.ctr], 
     col = "red", cex = 0.9, font = 4)    
legend("bottomleft", legend=ABnames, pch=tellus.pch, 
       col=tellus.col, text.col=tellus.col, pt.cex=0.6, cex=0.8)
# dev.off()

# png(file="Fig_8_right.png",width=7,height=5.5,units="in",res=144)
rescale <- 3 # for ellipses
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * 0.44*tellus.lra.rpc, type = "n", asp = 1, 
     xlab = paste("PCA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("PCA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
require(ellipse)

set.seed(123)
CIplot_biv(tellus.pca.rpc[,1], tellus.pca.rpc[,2], group=AB, groupcols=tellus.col, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           shownames=FALSE)
set.seed(123)
CIplot_biv(tellus.pca.rpc[,1], tellus.pca.rpc[,2], group=AB, groupcols=tellus.col, 
           add=TRUE, shade=FALSE, groupnames=ABnames, alpha=0.99)
# dev.off()

### study of the Box-Cox transformation in CA on the geometry of the parts
### should work with the columns: a 't' before 'tellus' indicates transposed
ttellus <- t(tellus.pro)
ttellus.pro <- CLOSE(ttellus)
# ttellus.clr <- CLR(ttellus.pro, weight=FALSE)
ttellus.lra <- LRA(ttellus.pro, weight=FALSE)
ttellus.lra.rpc <- ttellus.lra$rowpcoord 

### for original CA/chi-square
ttellus.ca <- CA(ttellus.pro)
ttellus.ca.rpc <- ttellus.ca$rowpcoord 
protest(ttellus.lra.rpc, ttellus.ca.rpc, permutations=0)$t0
# [1] 0.8663466

### Now CA/chi-square with power transformation on data with zeros replaced
### Sequence of powers down to the smallest 0.0001 
BoxCox <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- ttellus.pro^alpha
  foo.ca <- CA(foo)
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox[k] <- protest(ttellus.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

### Repeat with original data zeros not replaced
ttellus0 <- t(tellus0)
ttellus0.pro <- CLOSE(ttellus0)

BoxCox0 <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- ttellus0.pro^alpha
  foo.ca <- CA(foo)
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox0[k] <- protest(ttellus.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

### What is maximum Procrustes correlation
max(BoxCox0)
# [1] 0.9431539 

### For which power?
c(seq(1,0.01,-0.01),0.0001)[which(BoxCox0 == max(BoxCox0))]
# [1] 0.5

### -----------------------------------------------------------------------
### Figure 9: plots of Procrustes correlations for Box-Cox transformation 
###           CA for data with zeros replaced and data with original zeros
# pdf(file="Fig_9.pdf", width=5, height=5, useDingbats=FALSE, family="ArialMT")
par(mar=c(4.2,4,1,1), font.lab=2, las=1)
plot(c(seq(1,0.01,-0.01),0.0001), BoxCox, xlab="Power of Box-Cox transformation", 
     ylab="Procrustes correlation", type="l", lwd=2, col="blue", ylim=c(0.35,1),
     bty="n", xaxt="n", yaxt="n")
axis(1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(2)
lines(c(seq(1,0.01,-0.01),0.0001), BoxCox0, lwd=2, col="red", lty=3)
segments(0.5,0,0.5,BoxCox0[51], col="pink", lwd=2, lty=2)
legend("bottomright", legend=c("zeros replaced","with zeros"), 
       bty="n",
       col=c("blue","red"), 
       lwd=c(2,2), lty=c(1,3), cex=0.8)

# dev.off()
### -----------------------------------------------------------------------

### Dimension reduction with CA of square-root transformed compositions
### (doing it on columns as for Box-Cox, although makes no difference)
ttellus0.ca <- CA(CLOSE(ttellus0.pro^0.5))
ttellus0.ca$sv <- ttellus0.ca$sv/0.5
round(100*ttellus0.ca$sv^2/sum(ttellus0.ca$sv^2),3)
#  [1] 43.384 22.532  5.194  4.729  3.068  2.649  2.119  1.832  1.591  1.436  1.152  1.116  0.987  0.813
# [15]  0.729 

### note again: rows are columns, and columns are rows, 
### and axes are reversed to agree with previous biplots
ttellus0.ca.cpc <- -ttellus0.ca$colpcoord
ttellus0.ca.rcc <- -ttellus0.ca$rowcoord * sqrt(ttellus0.ca$rowmass)
ttellus0.ca.ctr <- (ttellus0.ca.rcc[,1]^2 > 1/nrow(ttellus0)) | (ttellus0.ca.rcc[,2]^2 > 1/nrow(ttellus0)) 
sum(ttellus0.ca.ctr)
[1] 22

### -----------------------------------------------------------------------
### Figure 10: left and right figures of the CA biplot, and Age Brackets
# save as png insert into PPT, and eventually save all as png
# png(file="Fig_10_left.png",width=7,height=5.5,units="in",res=144)
rescale <- 1.5 # for points
dim <- c(1,2)
perc.hor <- 43.4; perc.ver <- 22.5
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * rbind(ttellus0.ca.cpc, rescale*ttellus0.ca.rcc), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
arrows(0, 0, 0.95 * rescale * ttellus0.ca.rcc[ttellus0.ca.ctr, 1], 
       0.95 * rescale * ttellus0.ca.rcc[ttellus0.ca.ctr, 2], 
       length = 0.1, angle = 10, col = "pink")
points(ttellus0.ca.cpc, pch = tellus.pch[AB], col = tellus.col.alpha[AB], font = 1, cex = 0.5)
text(rescale * ttellus0.ca.rcc[ttellus0.ca.ctr,], labels = rownames(ttellus0)[ttellus0.ca.ctr], col = "red", 
     cex = 0.9, font = 4)    
legend("bottomleft", legend=ABnames, 
       pch=tellus.pch, col=tellus.col, text.col=tellus.col, pt.cex=0.6, cex=0.8)

# dev.off()


# png(file="Fig_10_right.png",width=7,height=5.5,units="in",res=144)
rescale <- 3  # for ellipses
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * 0.5*ttellus0.ca.cpc, type = "n", 
     asp = 1, xlab = paste("CA dimension ", dim[1], " (", 
     round(perc.hor, 1), "%)", sep = ""), ylab = paste("CA dimension ", dim[2], " (", 
     round(perc.ver, 1), "%)", sep = ""), main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
require(ellipse)
set.seed(123)
CIplot_biv(ttellus0.ca.cpc[,1], ttellus0.ca.cpc[,2], group=AB, groupcols=tellus.col, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           shownames=FALSE)
set.seed(123)
CIplot_biv(ttellus0.ca.cpc[,1], ttellus0.ca.cpc[,2], group=AB, groupcols=tellus.col, 
           add=TRUE, shade=FALSE, groupnames=ABnames, alpha=0.99)
# dev.off()

### Procrustes between case coordinates in LRA and corresponding ones in CA
protest(tellus.lra.rpc, ttellus0.ca.cpc, permutations=0)$t0
# [1] 0.9569627


### ---------------------------------------------------------------------------------------
### k-means clustering of LRA and ALR and CA coodinates and comparison (3-cluster solution)

### for LRA
set.seed(123)
lra.km3 <- kmeans(tellus.lra.rpc, centers=3, nstart=50, iter.max=200)
# cluster sizes
lra.km3$size
## [1]  798 4513 1488

### for PCA of ALRs (ref:Al)
set.seed(123)
pca.km3 <- kmeans(tellus.pca.rpc, centers=3, nstart=50, iter.max=200)
# cluster sizes
pca.km3$size
## [1]  842 4478 1479

# for CA of power transformed (square root)
set.seed(123)
ca.km3 <- kmeans(ttellus0.ca.cpc, centers=3, nstart=50, iter.max=200)
# cluster sizes
ca.km3$size
##  914 4396 1489

### tables of agreements
table(lra.km3$cluster, pca.km3$cluster)[c(2,3,1),c(2,3,1)]
#        2    3    1
#   2 4467    9   37
#   3    7 1470   11
#   1    4    0  794

(4467+1470+794) / 6799
# 0.9899985

table(lra.km3$cluster, ca.km3$cluster)[c(2,3,1),c(2,3,1)]
#        3    1    2
#   2 4369   45   99
#   3   17 1444   27
#   1   10    0  788

(4369+1444+788) / 6799
# 0.9708781

### adjusted Rand index
require(pdfCluster)
adj.rand.index(lra.km3$cluster, pca.km3$cluster) 
# [1] 0.9708841
adj.rand.index(lra.km3$cluster, ca.km3$cluster) 
# [1] 0.9142874


### Quasi-coherence

### Subcompositional incoherence exercise for regular CA using chi-square distance
### with and without square-root transformation

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
  dist(scale(prof,center=FALSE,scale=rootaveprof))  
}

procr.CA    <- matrix(0, nrow=100, ncol=44)
procr.CA.05 <- matrix(0, nrow=100, ncol=44)
 
tellus.cm  <- colMeans(tellus0.pro)
D.chi      <- as.matrix(chidist(tellus0.pro, 2))
D.chi.05   <- as.matrix(chidist(tellus0.pro^0.5, 2))

set.seed(1234567)
for(j in seq(44,4,-2)) {
  nparts <- j
  for(i in 1:100) {
# find the subcompositional parts  
    jparts <- sample(1:52, nparts)
    foo <- tellus0.pro[,jparts]
# remove parts all zeros
    allzero <- which(colSums(foo)==0)
    if(length(allzero)>0) {
      jparts <- jparts[-allzero]
      foo <- tellus.pro[,jparts]
    }
# incoherence in CA via MDS of distances  
    D <- as.dist(D.chi[jparts, jparts])
    D.05 <- as.dist(D.chi.05[jparts, jparts])
    D.rpc <- cmdscale(D, eig=TRUE, k=length(jparts)-1)$points
    D.rpc.05 <- cmdscale(D.05, eig=TRUE, k=length(jparts)-1)$points
# remove samples that may have all zeros
    allzero <- which(rowSums(foo)==0)
    if(length(allzero)>0) {
      foo <- foo[-allzero,]
    }
    D2 <- chidist(foo, 2)
    D2.05 <- chidist(foo^0.5, 2)
    D2.rpc <- cmdscale(D2, eig=TRUE, k=length(jparts)-1)$points
    D2.rpc.05 <- cmdscale(D2.05, eig=TRUE, k=length(jparts)-1)$points 
    procr.CA[i,j] <- protest(D2.rpc, D.rpc, permutations=0)$t0
    procr.CA.05[i,j] <- protest(D2.rpc.05, D.rpc.05, permutations=0)$t0
  }
}

procr.CA.quants <- apply(procr.CA, 2, quantile, c(0.025,0.975), na.rm=TRUE)
round(procr.CA.quants[,seq(4,44,2)],4)
#         [,1]  [,2]   [,3]  [,4]  [,5]   [,6]   [,7]  [,8]   [,9]  [,10]  [,11]  [,12]  [,13]
# 2.5%  0.8466 0.775 0.8001 0.833 0.829 0.9121 0.8948 0.894 0.9392 0.9611 0.9531 0.9583 0.9755
# 97.5% 0.9999 1.000 1.0000 1.000 1.000 1.0000 1.0000 1.000 1.0000 1.0000 1.0000 1.0000 1.0000
#       [,14]  [,15]  [,16]  [,17]  [,18]  [,19]  [,20]  [,21]
# 2.5%  0.983 0.9852 0.9796 0.9907 0.9904 0.9898 0.9898 0.9914
# 97.5% 1.000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000

procr.CA.quants.05 <- apply(procr.CA.05, 2, quantile, c(0.025,0.975), na.rm=TRUE)
round(procr.CA.quants.05[,seq(4,44,2)],4)
#         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]  [,7]   [,8]   [,9]  [,10]  [,11] [,12]  [,13]
# 2.5%  0.9908 0.9883 0.9914 0.9921 0.9916 0.9959 0.997 0.9973 0.9967 0.9984 0.9984 0.999 0.9987
# 97.5% 1.0000 1.0000 0.9999 1.0000 1.0000 1.0000 1.000 1.0000 1.0000 1.0000 1.0000 1.000 1.0000
#        [,14]  [,15]  [,16]  [,17]  [,18]  [,19]  [,20]  [,21]
# 2.5%  0.9995 0.9993 0.9995 0.9997 0.9997 0.9998 0.9997 0.9998
# 97.5% 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000

procr.CA.ones <- rep(0,44)
for(j in seq(4,44,2)) procr.CA.ones[j] <- sum(procr.CA[,j]>0.999)
procr.CA.ones[seq(4,44,2)]
# [1] 31 21 33 38 42 49 50 48 48 60 55 63 67 69 68 75 81 77 81 77 89

procr.CA.ones.05 <- rep(0,44)
for(j in seq(4,44,2)) procr.CA.ones.05[j] <- sum(procr.CA.05[,j]>0.999)
procr.CA.ones.05[seq(4,44,2)]
# [1]  64  54  64  61  70  75  77  80  84  88  93  97  96  99  99 100 100 100 100 100 100

### Figure 11: Levels of coherence for regular chi-square geometry
# pdf(file="Fig_11.pdf", width=5, height=5, useDingbats=FALSE, family="ArialMT")
par(mar=c(5,5,1,1), mgp=c(3.5,0.7,0), font.lab=2, las=1, mfrow=c(1,1))
plot(rep(seq(4,44,2), each=2), as.numeric(procr.CA.quants[,seq(4,44,2)]), xlab="Number of parts in subcomposition",
     ylab="Procrustes correlation", bty="n", xaxt="n", ylim=c(0.75, 1.02), type="n", font.lab=2, xlim=c(4,45))
axis(1, at=seq(4,44,2), labels=seq(4,44,2))
for(j in seq(4,44,2)) segments(j, procr.CA.quants[1,j], j, procr.CA.quants[2,j], col="blue", lwd=2)
eps <- 0.2
for(j in seq(4,44,2)) segments(j-eps, procr.CA.quants[1,j], j+eps, procr.CA.quants[1,j], col="blue", lwd=2, lend=2)
for(j in seq(4,44,2)) segments(j-eps, procr.CA.quants[2,j], j+eps, procr.CA.quants[2,j], col="blue", lwd=2, lend=2)
points(seq(4,44,2), apply(procr.CA[,seq(4,44,2)], 2, median, na.rm=TRUE), pch=21, col="blue", bg="white", cex=0.9)
text(seq(4,44,2), rep(1.01, 21), labels=procr.CA.ones[seq(4,44,2)], font=2, cex=0.6)
# dev.off()

### Figure 12: as before for sqrt profiles (plot window narrower vertically)
# pdf(file="Tellus_coherence_CAsqrt.pdf", width=5, height=2.7, useDingbats=FALSE, family="ArialMT")
par(mar=c(5,5,1,1), mgp=c(3.5,0.7,0), font.lab=2, las=1, mfrow=c(1,1))
plot(rep(seq(4,44,2), each=2), as.numeric(procr.CA.quants.05[,seq(4,44,2)]), xlab="Number of parts in subcomposition",
     ylab="Procrustes correlation", bty="n", xaxt="n", ylim=c(0.90, 1.01), type="n", font.lab=2, xlim=c(4,45), yaxt="n")
axis(1, at=seq(4,44,2), labels=seq(4,44,2))
axis(2, at=c(0.90, 0.95, 1.00), labels=c("0.90", "0.95", "1.00"))
for(j in seq(4,44,2)) segments(j, procr.CA.quants.05[1,j], j, procr.CA.quants.05[2,j], col="blue", lwd=2)
eps <- 0.2
for(j in seq(4,44,2)) segments(j-eps, procr.CA.quants.05[1,j], j+eps, procr.CA.quants.05[1,j], col="blue", lwd=2, lend=2)
for(j in seq(4,44,2)) segments(j-eps, procr.CA.quants.05[2,j], j+eps, procr.CA.quants.05[2,j], col="blue", lwd=2, lend=2)
points(seq(4,44,2), apply(procr.CA.05[,seq(4,44,2)], 2, median, na.rm=TRUE), pch=21, col="blue", bg="white", cex=0.9)
text(seq(4,44,2), rep(1.01, 21), labels=procr.CA.ones.05[seq(4,44,2)], font=2, cex=0.6)
# dev.off()

### for the subcomposition of rare earth minerals
jparts <- rare
D <- as.dist(D.chi[jparts, jparts])
D.rpc <- cmdscale(D, eig=TRUE, k=8)$points
foo <- tellus.pro[,jparts]
# remove parts all zeros
allzero <- which(colSums(foo)==0)
if(length(allzero)>0) {
  jparts <- jparts[-allzero]
  foo <- tellus.pro[,jparts]
}
# remove samples all zeros
allzero <- which(rowSums(foo)==0)
if(length(allzero)>0) {
  foo <- foo[-allzero,]
}
D2 <- chidist(foo, 2)
D2.rpc <- cmdscale(D2, eig=TRUE, k=8)$points
protest(D2.rpc, D.rpc, permutations=0)$t0
# [1] 0.9723521

### same, but with square root transformation
jparts <- rare
D.chi.sqrt  <- as.matrix(chidist(sqrt(tellus.pro), 2))
D.sqrt <- as.dist(D.chi.sqrt[jparts, jparts])
D.sqrt.rpc <- cmdscale(D.sqrt, eig=TRUE, k=8)$points
foo <- tellus.pro[,jparts]
# remove parts all zeros
allzero <- which(colSums(foo)==0)
if(length(allzero)>0) {
  jparts <- jparts[-allzero]
  foo <- tellus.pro[,jparts]
}
# remove samples all zeros
allzero <- which(rowSums(foo)==0)
if(length(allzero)>0) {
  foo <- foo[-allzero,]
}
D2.sqrt <- chidist(sqrt(foo), 2)
D2.sqrt.rpc <- cmdscale(D2.sqrt, eig=TRUE, k=8)$points
protest(D2.sqrt.rpc, D.sqrt.rpc, permutations=0)$t0
# [1] 0.9984683
