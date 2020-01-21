### set working directory to where the data reside
#setwd("C:/Users/Michael Greenacre/Documents/CA/CODA/Martin")

### package easyCODA needs to be installed
### required packages will be installed automatically
### install.packages("easyCODA", dependencies=TRUE)

library(easyCODA)

###  input FA data and show top left corner of file
FA <- read.csv("copepods.csv", header=T, check.names=F)
FA[1:3,1:6]
#      Season      14:0 14:1(n-5)    i-15:0    a-15:0
# 1 B5 winter 13.854102 0.2025212 1.1903286 0.4463732
# 2 B6 winter 11.826581 0.1479822 1.2358517 0.4599448
# 3 B7 winter  6.457139 0.0000000 0.7680038 0.2481243

### save first two columns and then remove them from FA data frame, then converted to a matrix
label  <- FA[,1]
season <- FA[,2] 
FA     <- FA[,-c(1,2)]
FA     <- as.matrix(FA)
dim(FA)
# [1] 42 40

### average percentages
round(colMeans(FA),1)
#      14:0  14:1(n-5)     i-15:0     a-15:0       15:0  15:1(n-6)     i-16:0       16:0  16:1(n-9) 
#       8.1        0.1        0.7        0.2        0.5        0.1        0.0        8.9        0.0 
# 16:1(n-7)  16:1(n-5)     i-17:0     a-17:0  16:2(n-4)       17:0  16:3(n-4)  16:4(n-1)       18:0 
#      12.4        0.7        0.2        0.1        0.5        0.1        0.6        0.7        3.0 
# 18:1(n-9)  18:1(n-7)  18:2(n-6)  18:3(n-6)  18:3(n-3)  18:4(n-3)       20:0 20:1(n-11)  20:1(n-9) 
#       5.8        0.9        1.6        0.2        0.9        6.3        0.6        0.1       13.3 
# 20:1(n-7)  20:2(n-6)  20:3(n-6)  20:4(n-6)  20:3(n-3)  20:4(n-3)  20:5(n-3) 22:1(n-11)  22:1(n-9) 
#       1.1        0.2        0.5        0.2        0.1        0.7       10.2        8.5        1.2 
# 22:1(n-7)  22:5(n-3)  22:6(n-3)  24:1(n-9) 
#       0.4        0.6        8.8        0.9 

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
#          FA1       FA2    R2
# 1       16:0 18:4(n-3) 54.28   <-   ratio chosen
# 2  18:1(n-7) 18:4(n-3) 54.27
# 3  18:4(n-3) 22:1(n-7) 53.32
# 4  16:1(n-9) 18:4(n-3) 53.24
# 5  16:2(n-4) 18:4(n-3) 52.98
# 6  16:1(n-7) 18:4(n-3) 52.64
# 7  18:4(n-3)      20:0 52.51
# 8       14:0 18:4(n-3) 52.35
# 9  18:4(n-3) 22:6(n-3) 52.35
# 10 18:4(n-3) 22:1(n-9) 52.29

### 1st ratio chosen by Martin 16:0/18:4(n-3)

### logratios1, 2, etc... will gather the logratio values chosen
### numratios1, 2, etc... will gather the numbers of the two parts 
logratios1 <- FA.step1$logratios.top[,1]
numratios1 <- FA.step1$ratios.top[1,]

### ------
### Step 2
FA.step2 <- STEP(FA0, nsteps = 1, top=20, previous=logratios1)
print.ratios(FA.step2$names.top, FA.step2$R2.top)
#          FA1        FA2    R2
# 1  16:1(n-7) 22:1(n-11) 75.90
# 2  16:1(n-7)  18:3(n-3) 75.82
# 3  16:1(n-7)  20:1(n-9) 75.72
# 4  16:1(n-7)  22:1(n-9) 75.53
# 5  16:1(n-7)  20:4(n-3) 75.46
# 6  16:1(n-7)  18:1(n-9) 75.42
# 7       16:0  16:1(n-7) 75.31
# 8  16:1(n-7)  18:4(n-3) 75.31
# 9  16:1(n-7)  16:1(n-5) 74.94
# 10    i-15:0  16:1(n-7) 74.83

### 7th ratio chosen by Martin 16:0/16:1(n-7)
### update logratios1 to logratios2, and numratios1 to numratios2
logratios2 <- cbind(logratios1, FA.step2$logratios.top[,7])
numratios2 <- rbind(numratios1, FA.step2$ratios.top[7,])

### ------
### Step 3
FA.step3 <- STEP(FA0, nsteps = 1, top=40, previous=logratios2)
print.ratios(FA.step3$names.top, FA.step3$R2.top,N=40)
# FA1       FA2    R2
# 1   20:1(n-9) 24:1(n-9) 82.85
# 2  22:1(n-11) 24:1(n-9) 82.79
# 3   22:1(n-9) 24:1(n-9) 82.76
# 4   16:3(n-4) 24:1(n-9) 82.75
# 5   18:1(n-9) 24:1(n-9) 82.75
# 6   18:3(n-3) 24:1(n-9) 82.75
# # 7   16:1(n-5) 24:1(n-9) 82.74
# 8   20:4(n-3) 24:1(n-9) 82.74
# 9        14:0 24:1(n-9) 82.74
# 10     i-15:0 24:1(n-9) 82.73
# 11       15:0 24:1(n-9) 82.73
# 12  16:2(n-4) 24:1(n-9) 82.68
# 13  16:1(n-7) 24:1(n-9) 82.65
# 14       16:0 24:1(n-9) 82.65
# 15  18:4(n-3) 24:1(n-9) 82.65
# 16       20:0 24:1(n-9) 82.64
# 17  18:1(n-7) 24:1(n-9) 82.61
# 18  16:1(n-9) 24:1(n-9) 82.59
# 19  18:2(n-6) 24:1(n-9) 82.56
# 20  22:1(n-7) 24:1(n-9) 82.56
# 21  20:1(n-7) 24:1(n-9) 82.56
# 22  20:5(n-3) 24:1(n-9) 82.53
# 23  22:6(n-3) 24:1(n-9) 82.53
# 24     i-16:0 24:1(n-9) 82.52
# 25     a-15:0 24:1(n-9) 82.51
# 26  20:3(n-3) 24:1(n-9) 82.42
# 27  15:1(n-6) 24:1(n-9) 82.41
# 28  20:4(n-6) 24:1(n-9) 82.40
# 29       18:0 24:1(n-9) 82.35
# 30  20:2(n-6) 24:1(n-9) 82.35
# 31  16:4(n-1) 24:1(n-9) 82.34
# 32  18:3(n-6) 24:1(n-9) 82.32
# 33  14:1(n-5) 24:1(n-9) 82.27
# 34     a-17:0 24:1(n-9) 82.18
# 35  20:3(n-6) 24:1(n-9) 82.05
# 36  20:1(n-9) 22:6(n-3) 81.94   <- this ratio chosen
# 37  22:5(n-3) 24:1(n-9) 81.87
# 38  20:1(n-9) 22:1(n-7) 81.85
# 39  18:1(n-7) 20:1(n-9) 81.76
# 40 22:1(n-11) 22:6(n-3) 81.72

### 36th ratio chosen by Martin 20:1(n-9)/22:6(n-3)
### update logratios2 to logratios3, and numratios2 to numratios3
logratios3 <- cbind(logratios2, FA.step3$logratios.top[,36])
numratios3 <- rbind(numratios2, FA.step3$ratios.top[36,])

### ------
### Step 4
### Note that 24:1(n-9) is always in the top listed FA ratio to enter, 
### but it is henceforth excluded as it appears only in traces (40th FA) 
### Notice the new form of call to STEP()
FA.step4 <- STEP(data=FA0[,-40],datatarget=FA0, nsteps=1, top=20, previous=logratios3)
print.ratios(FA.step4$names.top, FA.step4$R2.top)
#          FA1       FA2    R2
# 1  18:4(n-3) 22:6(n-3) 85.07
# 2       16:0 20:1(n-9) 85.07   <- this ratio chosen
# 3  16:1(n-7) 22:6(n-3) 85.07
# 4       16:0 22:6(n-3) 85.07
# 5  16:1(n-7) 20:1(n-9) 85.07
# 6  18:4(n-3) 20:1(n-9) 85.07
# 7       14:0 16:1(n-7) 84.89
# 8       14:0      16:0 84.89
# 9       14:0 18:4(n-3) 84.89
# 10 20:1(n-9) 20:5(n-3) 84.87

### 2nd ratio chosen by Martin 16:0/20:1(n-9)
### update logratios3 to logratios4, and numratios3 to numratios4
logratios4 <- cbind(logratios3, FA.step4$logratios.top[,2])
numratios4 <- rbind(numratios3, FA.step4$ratios.top[2,])

### ------
### Step 5
FA.step5 <- STEP(data=FA0[,-40],datatarget=FA0, nsteps=1, top=20, previous=logratios4)
print.ratios(FA.step5$names.top, FA.step5$R2.top)
#       FA1       FA2    R2
# 1    14:0 20:5(n-3) 88.41   <- ratio chosen
# 2  i-15:0 20:5(n-3) 88.36
# 3    14:0 22:1(n-7) 88.32
# 4    14:0 18:4(n-3) 88.11
# 5    14:0      16:0 88.11
# 6    14:0 16:1(n-7) 88.11
# 7    14:0 20:1(n-9) 88.11
# 8    14:0 22:6(n-3) 88.11
# 9    18:0 20:5(n-3) 88.10
# 10   15:0 20:5(n-3) 88.08

### 1st ratio chosen by Martin 14:0/20:5(n-3)
### update logratios4 to logratios5, and numratios4 to numratios5
logratios5 <- cbind(logratios4, FA.step5$logratios.top[,1])
numratios5 <- rbind(numratios4, FA.step5$ratios.top[1,])

### ------
### Step 6
FA.step6 <- STEP(data=FA0[,-40],datatarget=FA0, nsteps=1, top=20, previous=logratios5)
print.ratios(FA.step6$names.top, FA.step6$R2.top)
#          FA1       FA2    R2
# 1       14:0      18:0 91.04
# 2       18:0 20:5(n-3) 91.04   <- ratio chosen (ties for first position)
# 3       15:0      18:0 90.98
# 4     i-15:0      18:0 90.97
# 5       18:0 18:1(n-9) 90.94
# 6       16:0      18:0 90.92
# 7       18:0 20:1(n-9) 90.92
# 8  16:1(n-7)      18:0 90.92
# 9       18:0 22:6(n-3) 90.92
# 10      18:0 18:4(n-3) 90.92

### 2nd ratio chosen by Martin 18:0/20:5(n-3)
### update logratios5 to logratios6, and numratios5 to numratios6
logratios6 <- cbind(logratios5, FA.step6$logratios.top[,2])
numratios6 <- rbind(numratios5, FA.step6$ratios.top[2,])

### These are the ratios chosen in the 6 steps (numbers first, then names)
rownames(numratios6) <- paste("Step",1:6,sep="")
colnames(numratios6) <- c("FA1","FA2")
finalratios <- as.data.frame(cbind(numratios6, 
            Ratio=paste(colnames(FA0)[numratios6[,1]],"/",colnames(FA0)[numratios6[,2]],sep="")))
finalratios
#       FA1 FA2               Ratio
# Step1   8  24      16:0/18:4(n-3)
# Step2   8  10      16:0/16:1(n-7)
# Step3  27  39 20:1(n-9)/22:6(n-3)
# Step4   8  27      16:0/20:1(n-9)
# Step5   1  34      14:0/20:5(n-3)
# Step6  18  34      18:0/20:5(n-3)

colnames(logratios6) <- finalratios[,3]

### The 8 parts used in the 6 ratios
partsinratios <- sort(unique(as.numeric(numratios6)))
colnames(FA0)[partsinratios]
# [1] "14:0"      "16:0"      "16:1(n-7)" "18:0"      "18:4(n-3)" "20:1(n-9)" "20:5(n-3)" "22:6(n-3)"
# 
### the LRA of the full data set 
### low-contributing FAs are de-accentuated by plotting them in light red
### from the LRA object
FA0.lra <- LRA(FA0)
### .ccc are the column contribution coordinates
### .ctr TRUE for high contributor, otherwise FALSE
FA0.ccc <- FA0.lra$colcoord * sqrt(FA0.lra$colmass)
FA0.ctr <- (FA0.ccc[,1]^2 * FA0.lra$sv[1]^2 + FA0.ccc[,2]^2 * FA0.lra$sv[2]^2) / 
           (FA0.lra$sv[1]^2 + FA0.lra$sv[2]^2) > 1/ncol(FA0)
### only show parts for which FA0.ctr = TRUE (high contributors)
FA0.lra$colcoord <- FA0.lra$colcoord[FA0.ctr,]
FA0.lra$colmass  <- FA0.lra$colmass[FA0.ctr]
FA0.lra$colnames <- FA0.lra$colnames[FA0.ctr]

### season numbers and sample coordinates
season.num <- as.numeric(season)
season.col <- c("forestgreen","chocolate","blue")
season.pch <- c(23,22,24)
season.cex <- c(0.9,0.9,0.8)

FA0.rpc <- FA0.lra$rowpcoord

### plot the high contributors
 pdf(file="Figure3a_new.pdf", height=7, width=8.5)   # for saving file
par(mar=c(4.2,4,2,2), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
PLOT.LRA(FA0.lra, map="contribution")
### add the low contributors with small font and in pink
text(FA0.ccc[!FA0.ctr,], labels=colnames(FA0)[!FA0.ctr], col="pink", cex=0.6)
### add samples as coloured symbols for spring (green), summer (brown) and winter (blue)
points(FA0.rpc, pch=season.pch[season.num], col=season.col[season.num], 
       bg=season.col[season.num], cex=season.cex[season.num], font=2)
 dev.off()

### the PCA of the 6 ratios involving 8 FAs
logratios6.pca <- PCA(logratios6, weight=FALSE)
logratios6.rpc <- logratios6.pca$rowpcoord
 pdf(file="Figure3b_new.pdf", height=7, width=8.5) # for saving file
par(mar=c(4.2,4,2,2), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
PLOT.PCA(logratios6.pca, map="contribution")
### add samples as coloured symbols for spring (green), summer (brown) and winter (blue)
points(logratios6.rpc, pch=season.pch[season.num], col=season.col[season.num], 
       bg=season.col[season.num], cex=season.cex[season.num], font=2)
 dev.off()

