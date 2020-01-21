### set working directory to where the data reside
#setwd("C:/Users/Michael Greenacre/Documents/CA/CODA/Martin")

### package easyCODA needs to be installed, from CRAN or latest version from R-Forge
### required packages will be installed automatically
### install.packages("easyCODA")
###   Later, to update to latest version of easyCODA from R-Forge:
### install.packages("easyCODA", repos="http://R-Forge.R-project.org")

library(easyCODA)

###  input FA data and show top left corner of file
FA <- read.csv("amphipods.csv", header=T, check.names=F)
FA[1:3,1:6]
#      Species Stage Season   14:0  15:0   16:0
# 1 C.guilelmi     F Summer 1.1471     0 4.4584 
# 2 C.guilelmi     F Summer 4.6795     0 7.4466
# 3 C.guilelmi     F Summer 0.9734     0 5.7633

### save first three columns and then remove them from FA data frame, then converted to a matrix
sample.names <- FA[,1]
stage        <- FA[,2]
season       <- FA[,3] 
FA <- FA[,-c(1:3)]
FA <- as.matrix(FA)
rownames(FA) <- paste(substr(sample.names,1,3), substr(season,1,1), sep="")
dim(FA)
# [1] 52 27

### average percentages
round(colMeans(FA),1)
# 14:0       15:0       16:0  16:1(n-7)  16:1(n-5)  16:2(n-4)       17:0  16:3(n-4) 
# 5.1        0.1       12.8        7.1        0.3        0.4        0.3        0.5 
# 16:4(n-1)       18:0  18:1(n-9)  18:1(n-7)  18:2(n-6)  18:3(n-3)  18:4(n-3) 20:1(n-11) 
# 0.5        1.2       14.7        2.4        1.5        1.1        4.0        3.3 
# 20:1(n-9)  20:1(n-7)  20:4(n-6)  20:4(n-3)  20:5(n-3) 22:1(n-11)  22:1(n-9)  22:1(n-7) 
# 13.5        1.1        0.4        2.3        8.7        7.0        1.7        0.4 
# 22:5(n-3)  24:1(n-9)  22:6(n-3) 
# 0.5        0.6        8.6 

### replace zeros with NAs in new object FA_NA and find minimum positive FA values
FA_NA <- FA
FA_NA[FA==0] <- NA
FAmin <- apply(FA_NA, 2, min, na.rm=T)

### FA0 has 0s replaced with half the minimum value for each FA
FA0 <- FA
for(j in 1:ncol(FA)) {
  for(i in 1:nrow(FA)) {
    if(FA[i,j]==0) FA0[i,j] <- 0.5 * FAmin[j]
  }
}

### reclose (i.e. renormalize) the data
FA0 <- FA0 / rowSums(FA0)

### function print.ratios() to make output for supplementary material
print.ratios <- function(rationames, R2, procr=NA, N=10) {
# function prints ratios and the corresponding R2, optionally Procrustes correlations
# prints first 10 ratios by default  
# split names of parts in ratio into FA1 and FA2
# notice that the ratios can be reported as FA1/FA2 or FA2/FA1  
  foo    <- strsplit(rationames, "/")
  parts  <- matrix(unlist(foo), ncol=2, byrow=TRUE)
  df   <- as.data.frame(parts)[1:N,]
  if(is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2))
    colnames(df) <- c("FA1", "FA2","R2")
  }
  if(!is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2), round(procr[1:N], 3))
    colnames(df) <- c("FA1", "FA2", "R2","Procrustes")
  }
  print(df[1:N,])
}

### ------
### Step 1
FA.step1 <- STEP(FA0, nsteps = 1, top=20)
print.ratios(FA.step1$names.top, FA.step1$R2.top)
#           FA1        FA2    R2
# 1   18:4(n-3) 22:1(n-11) 42.07
# 2   18:4(n-3)  22:1(n-9) 41.81
# 3   20:5(n-3) 22:1(n-11) 40.95  <- ratio chosen
# 4   20:5(n-3)  22:1(n-9) 38.97
# 5   20:1(n-9)  20:4(n-3) 38.72
# 6  22:1(n-11)  22:6(n-3) 38.15
# 7   18:4(n-3)  20:1(n-9) 38.13
# 8   18:1(n-9)  20:4(n-3) 38.03
# 9   22:1(n-9)  22:6(n-3) 37.96
# 10  18:1(n-9)  18:4(n-3) 37.77

### 3rd ratio chosen by Martin 20:5(n-3)/22:1(n-11)

### logratios1, 2, etc... will gather the logratio values chosen
### numratios1, 2, etc... will gather the numbers of the two parts 
logratios1 <- FA.step1$logratios.top[,3]
numratios1 <- FA.step1$ratios.top[3,]

### ------
### Step 2
FA.step2 <- STEP(FA0, nsteps = 1, top=20, previous=logratios1)
print.ratios(FA.step2$names.top, FA.step2$R2.top)
#          FA1        FA2    R2
# 1       16:0  20:5(n-3) 56.78
# 2       16:0 22:1(n-11) 56.78  <- ratio chosen
# 3       16:0  18:4(n-3) 56.26
# 4       18:0  18:4(n-3) 56.13
# 5       18:0 22:1(n-11) 55.83
# 6       18:0  20:5(n-3) 55.83
# 7       14:0  18:4(n-3) 55.62
# 8       16:0  22:6(n-3) 55.50
# 9  18:1(n-9)  18:4(n-3) 55.44
# 10      18:0  22:6(n-3) 55.01

### 2nd ratio chosen by Martin 16:0/22:1(n-11)
### update logratios1 to logratios2, and numratios1 to numratios2
logratios2 <- cbind(logratios1, FA.step2$logratios.top[,2])
numratios2 <- rbind(numratios1, FA.step2$ratios.top[2,])

### ------
### Step 3
FA.step3 <- STEP(FA0, nsteps = 1, top=20, previous=logratios2)
print.ratios(FA.step3$names.top, FA.step3$R2.top)
#          FA1        FA2    R2
# 1  18:1(n-9)  20:4(n-3) 69.86
# 2       18:0  20:4(n-3) 69.78
# 3  20:4(n-3)  22:6(n-3) 69.66
# 4  18:4(n-3)  22:6(n-3) 69.45
# 5       18:0  18:4(n-3) 69.34  <- ratio chosen
# 6       18:0  20:1(n-9) 69.28
# 7  18:1(n-7)  20:4(n-3) 69.25
# 8  18:1(n-9)  18:4(n-3) 69.08
# 9  20:4(n-3)  20:5(n-3) 68.93
# 10 20:4(n-3) 22:1(n-11) 68.93

### 5th ratio chosen by Martin 18:0/18:4(n-3)
### update logratios2 to logratios3, and numratios2 to numratios3
logratios3 <- cbind(logratios2, FA.step3$logratios.top[,5])
numratios3 <- rbind(numratios2, FA.step3$ratios.top[5,])

### ------
### Step 4
FA.step4 <- STEP(FA0, nsteps = 1, top=20, previous=logratios3)
print.ratios(FA.step4$names.top, FA.step4$R2.top)
#          FA1        FA2    R2
# 1  20:4(n-3)  22:6(n-3) 77.63
# 2       18:0  20:4(n-3) 77.59  <- ratio chosen
# 3  18:4(n-3)  20:4(n-3) 77.59
# 4  20:4(n-3)  20:5(n-3) 77.39
# 5  20:4(n-3) 22:1(n-11) 77.39
# 6       16:0  20:4(n-3) 77.39
# 7  20:1(n-9) 22:1(n-11) 77.37
# 8       16:0  20:1(n-9) 77.37
# 9  20:1(n-9)  20:5(n-3) 77.37
# 10 20:1(n-9)  22:6(n-3) 77.01

### 2nd ratio chosen by Martin 18:0/20:4(n-3)
### update logratios3 to logratios4, and numratios3 to numratios4
logratios4 <- cbind(logratios3, FA.step4$logratios.top[,2])
numratios4 <- rbind(numratios3, FA.step4$ratios.top[2,])

### ------
### Step 5
FA.step5 <- STEP(FA0, nsteps = 1, top=20, previous=logratios4)
print.ratios(FA.step5$names.top, FA.step5$R2.top)
#          FA1        FA2    R2
# 1  18:1(n-9)  20:1(n-9) 83.85   <- ratio chosen
# 2       14:0  18:1(n-9) 83.45
# 3  18:1(n-9)  20:5(n-3) 83.22
# 4       16:0  18:1(n-9) 83.22
# 5  18:1(n-9) 22:1(n-11) 83.22
# 6  20:1(n-9)  22:6(n-3) 82.72
# 7       14:0  22:6(n-3) 82.53
# 8  18:2(n-6)  20:1(n-9) 82.52
# 9  18:1(n-7)  20:1(n-9) 82.39
# 10      14:0  18:2(n-6) 82.39

### 1st ratio chosen by Martin 18:1(n-9)/20:1(n-9)
### update logratios4 to logratios5, and numratios4 to numratios5
logratios5 <- cbind(logratios4, FA.step5$logratios.top[,1])
numratios5 <- rbind(numratios4, FA.step5$ratios.top[1,])

### ------
### Step 6
FA.step6 <- STEP(FA0, nsteps = 1, top=20, previous=logratios5)
print.ratios(FA.step6$names.top, FA.step6$R2.top)
#          FA1        FA2    R2
# 1  16:1(n-7)  22:6(n-3) 87.47   <- ratio chosen
# 2  16:1(n-7)       18:0 87.04
# 3  16:1(n-7)  20:4(n-3) 87.04
# 4  16:1(n-7)  18:4(n-3) 87.04
# 5       14:0  22:6(n-3) 86.98
# 6  18:1(n-9)  22:6(n-3) 86.97
# 7  20:1(n-9)  22:6(n-3) 86.97
# 8       16:0  16:1(n-7) 86.96
# 9  16:1(n-7) 22:1(n-11) 86.96
# 10 16:1(n-7)  20:5(n-3) 86.96

### 1st ratio chosen by Martin 16:1(n-7)/22:6(n-3)
### update logratios5 to logratios6, and numratios5 to numratios6
logratios6 <- cbind(logratios5, FA.step6$logratios.top[,1])
numratios6 <- rbind(numratios5, FA.step6$ratios.top[1,])

### ------
### Step 7
FA.step7 <- STEP(FA0, nsteps = 1, top=20, previous=logratios6)
print.ratios(FA.step7$names.top, FA.step7$R2.top)
#           FA1        FA2    R2
# 1  20:1(n-11)  20:5(n-3) 89.77
# 2        16:0 20:1(n-11) 89.77
# 3  20:1(n-11) 22:1(n-11) 89.77
# 4        14:0 20:1(n-11) 89.76
# 5  20:1(n-11)  22:6(n-3) 89.65
# 6   16:1(n-7) 20:1(n-11) 89.65
# 7   18:1(n-9)  20:5(n-3) 89.63
# 8        16:0  18:1(n-9) 89.63
# 9   20:1(n-9) 22:1(n-11) 89.63   <- ratio chosen
# 10  18:1(n-9) 22:1(n-11) 89.63

### 9th ratio chosen by Martin 20:1(n-9)/22:1(n-11)
### update logratios6 to logratios7, and numratios6 to numratios7
logratios7 <- cbind(logratios6, FA.step7$logratios.top[,9])
numratios7 <- rbind(numratios6, FA.step7$ratios.top[9,])

### ------
### Step 8
FA.step8 <- STEP(FA0, nsteps = 1, top=20, previous=logratios7)
print.ratios(FA.step8$names.top, FA.step8$R2.top)
#           FA1        FA2    R2
# 1   18:1(n-9) 20:1(n-11) 91.57
# 2  20:1(n-11)  20:1(n-9) 91.57
# 3        16:0 20:1(n-11) 91.57
# 4  20:1(n-11) 22:1(n-11) 91.57
# 5  20:1(n-11)  20:5(n-3) 91.57
# 6  20:1(n-11)  22:6(n-3) 91.56
# 7   16:1(n-7) 20:1(n-11) 91.56
# 8        14:0 20:1(n-11) 91.49
# 9   18:1(n-7) 20:1(n-11) 91.45
# 10 20:1(n-11)  22:1(n-9) 91.44

### 7th ratio chosen by Martin 16:1(n-7)/20:1(n-11)
### update logratios7 to logratios8, and numratios7 to numratios8
logratios8 <- cbind(logratios7, FA.step8$logratios.top[,7])
numratios8 <- rbind(numratios7, FA.step8$ratios.top[7,])

### These are the ratios chosen in the 8 steps (numbers first, then names)
rownames(numratios8) <- paste("Step",1:8,sep="")
colnames(numratios8) <- c("FA1","FA2")
finalratios <- as.data.frame(cbind(numratios8, 
            Ratio=paste(colnames(FA0)[numratios8[,1]],"/",colnames(FA0)[numratios8[,2]],sep="")))
finalratios
#       FA1 FA2                Ratio
# Step1  21  22 20:5(n-3)/22:1(n-11)
# Step2   3  22      16:0/22:1(n-11)
# Step3  10  15       18:0/18:4(n-3)
# Step4  10  20       18:0/20:4(n-3)
# Step5  11  17  18:1(n-9)/20:1(n-9)
# Step6   4  27  16:1(n-7)/22:6(n-3)
# Step7  17  22 20:1(n-9)/22:1(n-11)
# Step8   4  16 16:1(n-7)/20:1(n-11)

colnames(logratios8) <- finalratios[,3]

### The 11 parts used in the 8 ratios
partsinratios <- sort(unique(as.numeric(numratios8)))
colnames(FA0)[partsinratios]
# [1] "16:0"       "16:1(n-7)"  "18:0"       "18:1(n-9)"  "18:4(n-3)"  "20:1(n-11)" "20:1(n-9)" 
# [8] "20:4(n-3)"  "20:5(n-3)"  "22:1(n-11)" "22:6(n-3)"

### the LRA of the full data set 
### low-contributing FAs are excluded by eliminating their coordinates
### from the LRA object
FA0.lra <- LRA(FA0)
### .ccc are the column contribution coordinates
### .ctr TRUE for high contributor, otherwise FALSE
FA0.ccc <- FA0.lra$colcoord * sqrt(FA0.lra$colmass)
FA0.ctr <- (FA0.ccc[,1]^2 * FA0.lra$sv[1]^2 + FA0.ccc[,2]^2 * FA0.lra$sv[2]^2) / 
           (FA0.lra$sv[1]^2 + FA0.lra$sv[2]^2) > 1/ncol(FA0)
### only retain parts for which FA0.ctr = TRUE (high contributors)
FA0.lra$colcoord <- FA0.lra$colcoord[FA0.ctr,]
FA0.lra$colmass  <- FA0.lra$colmass[FA0.ctr]
FA0.lra$colnames <- FA0.lra$colnames[FA0.ctr]

### plot the high contributors
# pdf(file="Figure7a.pdf", height=7, width=7.5)   # for saving file
par(mar=c(4.2,4,2,2), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
PLOT.LRA(FA0.lra, axes.inv=c(-1,-1), map="contribution", rescale=2)
# dev.off()

### the PCA of the 8 ratios involving 11 FAs
logratios8.pca <- PCA(logratios8, weight=FALSE)
# pdf(file="Figure7b.pdf", height=6.5, width=8.5) # for saving file
par(mar=c(4.2,4,2,2), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
PLOT.PCA(logratios8.pca, map="contribution", axes.inv=c(1,-1), rescale=2)
# dev.off()

