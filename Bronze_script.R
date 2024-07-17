### packages required
require(RColorBrewer)
require(easyCODA)
require(boot)
require(ellipse)
require(rpart)
require(rpart.plot)
require(nnet)
require(splitstackshape)
require(randomForest)

### read data from all three dynasties
BRONZE <- read.csv("Shang_West_East.csv", header=TRUE)
dim(BRONZE)
# [1] 362  12
head(BRONZE)
#   Cat.No. Chronological_typology  GROUP GROUP2     Cu     Sn     Pb Zn   Au   Ag    As  Sb
# 1       1                  -1450 ShangC  Shang 851000  86000  53000 28  9.8 2340  1720 490
# 2       2                  -1400 ShangB  Shang 728000  63000 194000 17 10.6  990  2240 219
# 3      80                  -1400 ShangB  Shang 803000 107000  83000 23 16.5 1180  3600 307
# 4       3                  -1325 ShangB  Shang 950000  22100  27300 37  9.4  760  2430 550
# 5       4                  -1325 ShangB  Shang 918000  45000  27600 84 19.5  470  3800 274
# 6      51                  -1325 ShangB  Shang 669000  98000 211000 49 37.0 1090 12700 790


### (we use the classification GROUP2, not the finer classification of Shang included in GROUP)
table(BRONZE$GROUP2)
# Eastern Zhou        Shang Western Zhou 
#          96          109          157 
### the dynasty time order is Shang / Western Zhou / Eastern Zhou
### names, colours, symbols, sizes and numbers for the three dynasty groups
group.names <- c("Shang","WZhou","EZhou")
group.col <- brewer.pal(3, "Dark2")[c(2,1,3)]
group.pch <- c(4,21,22)[c(2,1,3)]
group.cex <- c(0.5,0.5,0.5)[c(2,1,3)]
group <- factor(BRONZE$GROUP2, levels=c("Shang","Western Zhou","Eastern Zhou"))
group.num <- as.numeric(group)
table(group.num)
#   1   2   3 
# 109 157  96 
# Note: the order in group.num is now time order

### element data and compositions
bronze <- BRONZE[,5:12]
bronze.pro <- CLOSE(bronze)

### ---------------------------------------------
### Compositional barplot summary of compositions
### ---------------------------------------------
### Use Dark2 palette in RColorBrewer for the elements
### The 3x2 array group.ends gives the starting and ending sample numbers
### of the time-ordered three dynasties in the data file
group.ends <- matrix(c(1,109, 110,266, 267,362), nrow=3, byrow=TRUE)
colnames(group.ends) <- c("start","end")
rownames(group.ends) <- group.names 
element.order <- order(colMeans(bronze.pro))
element.col <- rev(c(brewer.pal(8,"Dark2")[1:6], "yellow2",brewer.pal(8,"Dark2")[7])) 

### get all cumulative sums of the data in element.order
bronze.pro.cum     <- t(apply(cbind(rep(0,nrow(bronze)), bronze.pro[,element.order]), 1, cumsum)) 

### Figure 2 (pull out window so it is wide and flat, adjusting accordingly)
par(mar=c(4,1,2,1), mgp=c(2,0.7,0), font.lab=2, cex.lab=0.9, las=1)
plot(0,0,xlab="Percentage (%)",ylab="",xlim=c(-3,100),ylim=c(-35,366), yaxt="n", bty="n", main="Bronze data (all elements)", type="n")
nlines <- nrow(bronze)+2
i.extra <- 0
ii <- 1
for(igp in 1:3) {
  for(i in group.ends[igp,1]:group.ends[igp,2]) {
    foo <- 100*bronze.pro.cum[i,]
    for(j in 1:ncol(bronze)) segments(foo[j], ii, foo[j+1], ii, col=element.col[j],lwd=3,lend=1)
    ii <- ii+1
  }
  i.extra <- i.extra+1
  if(igp!=3) segments(0, ii+i.extra, 100, ii+i.extra, col="white",lwd=4, lend=1)
  ii <- ii+i.extra+1
}
foo <- cumsum(c(0,100*colMeans(bronze.pro[,element.order])))
for(j in 1:ncol(bronze)) segments(foo[j], -15, foo[j+1], -15, col=element.col[j],lwd=8, lend=1)
text(c(0,0,0,0,0),c(55,181,318) , labels=group.names, col=group.col[c(2,1,3)],pos=2, font=2, cex=0.8)
text(0, -15, labels="Average", cex=0.8, font=3, pos=2)
### just label top three elements
labels.mid <- foo[6:8] + (foo[7:9] - foo[6:8])/2
text(labels.mid, -30, labels=colnames(bronze)[element.order][6:8], col=element.col[6:8], cex=0.8, font=4)

### plot the 5 small percentages separately
### get subcomposition (first 5 in element.order) and cumulative sums
bronze.sub.pro <- CLOSE(bronze.pro[,element.order][,1:5])
element.sub.order <- order(colMeans(bronze.sub.pro))
bronze.sub.pro.cum     <- t(apply(cbind(rep(0,nrow(bronze)), bronze.sub.pro[,element.sub.order]), 1, cumsum)) 

### Figure 3 (pull out window so it is wide and flat, adjusting accordingly)
par(mar=c(4,1,2,1), mgp=c(2,0.7,0), font.lab=2, cex.lab=0.9, las=1)
plot(0,0,xlab="Percentage (%)",ylab="",xlim=c(-3,100),ylim=c(-35,366), yaxt="n", bty="n", main="Bronze data (5 trace elements as subcomposition)", type="n")
nlines <- nrow(bronze)+2
i.extra <- 0
ii <- 1
for(igp in 1:3) {
  for(i in group.ends[igp,1]:group.ends[igp,2]) {
    foo <- 100*bronze.sub.pro.cum[i,]
    for(j in 1:5) segments(foo[j], ii, foo[j+1], ii, col=element.col[j],lwd=3,lend=1)
    ii <- ii+1
  }
  i.extra <- i.extra+1
  if(igp!=3) segments(0, ii+i.extra, 100, ii+i.extra, col="white",lwd=5, lend=1)
  ii <- ii+i.extra+1
}

foo <- cumsum(c(0,100*colMeans(bronze.sub.pro[,element.sub.order])))
for(j in 1:1:6) segments(foo[j], -15, foo[j+1], -15, col=element.col[j],lwd=8, lend=1)
text(c(0,0,0,0,0),c(55,181,318) , labels=group.names, col=group.col[c(2,1,3)],pos=2, font=2, cex=0.8)
text(0, -15, labels="Average", cex=0.8, font=3, pos=2)
### just label bottom 5 elements
labels.mid <- foo[1:6] + (foo[2:7] - foo[1:6])/2
text(labels.mid, -30, labels=colnames(bronze.sub.pro)[element.sub.order], col=element.col[1:5], cex=0.8, font=4)


### ---------------------------------
### Univariate summaries of logratios
### ---------------------------------
### Total variance computed using pairwise logratios and CLRs
### Both ways gives identical result
### All J*(J-1)/2 pairwise logratios
### all pairwise logratios (PLRs)-- function is LR(), LRVAR() computes total and individual variances
bronze.PLR <- LR(bronze.pro, weight=FALSE)
bronze.PLR.var <- LR.VAR(bronze.PLR, vars=TRUE)
bronze.PLR.var$LRtotvar
# [1] 0.8020681
sum(bronze.PLR.var$LRvars)
# [1] 0.8020681
### The J centered logratios (CLRs)
bronze.CLR <- CLR(bronze.pro, weight=FALSE)
bronze.CLR.var <- LR.VAR(bronze.CLR, vars=TRUE)
bronze.CLR.var$LRtotvar
# [1] 0.8020681
sum(bronze.CLR.var$LRvars)
# [1] 0.8020681

### Percentage contributions of LRs and CLRs to total variance
bronze.PLR.ctr <- bronze.PLR.var$LRvars / sum(bronze.PLR.var$LRvars) 
bronze.PLR.ctr.order <- order(bronze.PLR.ctr, decreasing=TRUE)
round(100*bronze.PLR.ctr[bronze.PLR.ctr.order],2)
# Table 1, column (i) -- top 10 only
# Pb/Zn Pb/Au Sn/Pb Zn/Sb Cu/Pb Zn/Au Pb/As Pb/Sb Sn/Sb Zn/As Sn/Zn Zn/Ag Cu/Sb Au/Sb 
#  7.73  6.41  5.64  5.58  5.23  4.94  4.28  4.26  4.24  4.22  4.08  3.76  3.71  3.67 
# Pb/Ag Sn/Au Sn/As Au/As Cu/Zn Cu/Au Cu/As Ag/Sb Au/Ag Sn/Ag Ag/As As/Sb Cu/Ag Cu/Sn 
#  3.55  3.48  3.43  2.98  2.95  2.72  2.61  2.05  1.83  1.74  1.70  1.35  1.08  0.75 
bronze.CLR.ctr <- bronze.CLR.var$LRvars / sum(bronze.CLR.var$LRvars) 
bronze.CLR.ctr.order <- order(bronze.CLR.ctr, decreasing=TRUE)
# These not listed in article, but will be important in PCA of CLRs (i.e., LRA)
#    Pb    Zn    Au    Sb    Sn    As    Cu    Ag 
# 24.61 20.76 13.53 12.36 10.87  8.09  6.56  3.22 

### Percentages of logratio variance explained by each pairwise logratio
bronze.PLR.expvar <- rep(0, ncol(bronze.PLR$LR))
for(j in 1:ncol(bronze.PLR$LR)) {
  foo <- rda(bronze.PLR$LR ~ bronze.PLR$LR[,j]) 
  bronze.PLR.expvar[j] <- foo$CCA$tot.chi / foo$tot.chi
}
names(bronze.PLR.expvar) <- colnames(bronze.PLR$LR)
bronze.PLR.expvar.order <- order(bronze.PLR.expvar, decreasing=TRUE)
round(100*bronze.PLR.expvar[bronze.PLR.expvar.order],2)
# Table 1, column (ii) -- top 10 only
# Pb/Zn Cu/Pb Sn/Pb Zn/Sb Pb/Au Pb/Ag Zn/As Cu/Sb Sn/Sb Zn/Ag Pb/As Cu/As Zn/Au Sn/As 
# 33.44 31.96 30.49 28.27 28.26 27.45 25.54 25.23 25.02 23.77 23.46 22.28 21.61 21.28 
# Pb/Sb Au/Sb Sn/Zn Cu/Zn Cu/Ag Ag/Sb Sn/Au Au/As Cu/Au Sn/Ag Au/Ag Ag/As As/Sb Cu/Sn 
# 21.04 19.65 19.06 18.77 18.14 17.98 16.71 16.45 16.34 15.89 14.11 12.48  8.37  5.92  

### Order by F statistics in an ANOVA
bronze.PLR.aov <- rep(0, ncol(bronze.PLR$LR)) 
for(jk in 1:ncol(bronze.PLR$LR)) {
  foo <- aov(bronze.PLR$LR[,jk] ~ factor(group))
  bronze.PLR.aov[jk] <- summary(foo)[[1]]$F[1]
} 
names(bronze.PLR.aov) <- colnames(bronze.PLR$LR)
bronze.PLR.aov.order <- order(bronze.PLR.aov, decreasing=TRUE)
round(bronze.PLR.aov[bronze.PLR.aov.order],2)
# Table 1, column (iii) -- top 10 only
# Pb/Zn Cu/Sb Pb/Sb Zn/Sb Cu/Pb Pb/Au Cu/Ag Sn/Pb Sn/Sb Pb/As Zn/As Pb/Ag Cu/As Zn/Ag 
# 45.56 42.97 42.38 42.17 41.50 40.16 39.89 37.06 34.42 30.50 30.05 27.21 25.82 25.41 
# Cu/Au Ag/Sb Sn/Ag Au/Ag Zn/Au Sn/As Sn/Au As/Sb Au/Sb Au/As Cu/Zn Ag/As Sn/Zn Cu/Sn 
# 	24.88 21.39 20.73 19.49 18.93 17.53 17.44 11.19  8.32  5.53  4.25  4.15  3.33  0.11

### ----------------------------------------------------------------------------------
### Confidence plots comparing logratios of Pb/Zn, Cu/Sb and Pb/Sb across three groups
### ----------------------------------------------------------------------------------
### (also show their means and ANOVAs here)

PbZn <- log(bronze$Pb/bronze$Zn)
aggregate(PbZn ~ factor(group), FUN=mean)
#   factor(group)     PbZn
# 1         Shang 6.602638
# 2  Western Zhou 5.659250
# 3  Eastern Zhou 7.865428

summary(aov(PbZn ~ factor(group)))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# factor(group)   2  290.8  145.41   45.56 <2e-16 ***
# Residuals     359 1145.7    3.19                            
                    
CuSb <- log(bronze$Cu/bronze$Sb)
aggregate(CuSb ~ factor(group), FUN=mean)
#   factor(group)     CuSb
# 1         Shang 8.286858
# 2  Western Zhou 7.558402
# 3  Eastern Zhou 6.671532

summary(aov(CuSb ~ factor(group)))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# factor(group)   2  133.3   66.63   42.97 <2e-16 ***
# Residuals     359  556.6    1.55                   

PbSb <- log(bronze$Pb/bronze$Sb)
aggregate(PbSb ~ factor(group), FUN=mean)
#   factor(group)     PbSb
# 1         Shang 5.553224
# 2  Western Zhou 4.040664
# 3  Eastern Zhou 4.899196

summary(aov(PbSb ~ factor(group)))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# factor(group)   2  151.2   75.60   42.38 <2e-16 ***
# Residuals     359  640.4    1.78      

### Figure 4, CIplot_uni function is on GitHub or www.multivariatestatistics.org
par(mar=c(4.2,4,2,1), mgp=c(2,0.7,0), font.lab=2, mfrow=c(1,3), cex.lab=1.2)
CIplot_uni(PbZn, group=group.num, cols=group.col, names=c("Shang","WZhou","EZhou"), ylab="log(Pb/Zn)",  
           shade=TRUE, ylim=c(5,8.5), thin=1.4, CONF=0.99)
CIplot_uni(CuSb, group=group.num, cols=group.col, names=c("Shang","WZhou","EZhou"), ylab="log(Cu/Sb)", 
           shade=TRUE, ylim=c(6.5,9.0), thin=1.4, CONF=0.99)
CIplot_uni(PbSb, group=group.num, cols=group.col, names=c("Shang","WZhou","EZhou"), ylab="log(Pb/Sb)", 
           shade=TRUE, ylim=c(3.5,6), thin=1.4, CONF=0.99)

### -----------------------
### Logratio Analysis (LRA)
### -----------------------
bronze.pro <- as.matrix(CLOSE(bronze)) 
rownames(bronze.pro) <- group.num
bronze.LRA <- LRA(bronze.pro, weight=FALSE)
### Row principal coordinates
bronze.rpc <- bronze.LRA$rowpcoord
### Column standard coordinates 
bronze.csc <- bronze.LRA$colcoord
### Reverse the second dimension
bronze.rpc[,2] <- -bronze.rpc[,2]
bronze.csc[,2] <- -bronze.csc[,2]
percs   <- 100 * bronze.LRA$sv^2 / sum(bronze.LRA$sv^2)
sum(bronze.LRA$sv^2)
# [1] 0.8020681       # agrees with total variance computed before in different ways
### Plot Figure 5 "by hand"
par(mar=c(3.8,3.6,2.5,1), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.7)
plot(1.05 * rbind(bronze.rpc, bronze.csc), type = "n", asp = 1,
     xlab = paste("LRA dimension ", 1, " (", round(percs[1], 1), "%)", sep = ""), 
     ylab = paste("LRA dimension ", 2, " (", round(percs[2], 1), "%)", sep = ""),
     main = "Logratio Analysis")
abline(h = 0, v = 0, col = "gray", lty = 2)
points(bronze.rpc, pch = group.pch[group.num], col = group.col[group.num], font = 1, 
       cex = group.cex[group.num])
text(bronze.csc, labels = colnames(bronze), col = "red", cex = 0.9, font = 4)    
legend("bottomleft", legend=c("S:Shang","WZ:WZou","EZ:EZhou"), bty="n",
       pch=group.pch, col=group.col, pt.cex=group.cex, 
       text.col=group.col, cex=0.9, text.font=2)
for (j in 1:3){
 points     <- bronze.rpc[group.num==j,1:2]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,], lty = 3, col=group.col[j])
}

set.seed(123)
CIplot_biv(bronze.rpc[,1], bronze.rpc[,2], group=group.num, groupcols=group.col, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           groupnames=c("S","WZ","EZ"))
set.seed(123)
CIplot_biv(bronze.rpc[,1], bronze.rpc[,2], group=group.num, groupcols=group.col, 
           alpha=0.99, add=TRUE, shade=FALSE, groupnames=c("S","WZ","EZ"))


### ------------------------------------------------------------
### Non-hierarchical k-means cluster analysis of all the samples
### ------------------------------------------------------------
### Looping on k-means algorithm to decide how many clusters
set.seed(1234)
bronze.BW <- rep(0, 10)
for(nc in 2:10) {
  bronze.kmeans <- kmeans(bronze.CLR$LR, centers=nc, nstart=20, iter.max=200)
  bronze.BW[nc] <- bronze.kmeans$betweenss/bronze.kmeans$totss
}
bronze.BW
#   [1] 0.0000000 0.2362804 0.3680456 0.4422665 0.4924735 0.5326702 0.5591886 0.5872021 0.6136536 0.6345105
### Plot the proportion of between-cluster variance and the increments in between-cluster variance
### Figure 6 barplots for deciding number of clusters
par(mar=c(4.2,4,1,2), mgp=c(2.5,0.7,0), las=1, font.lab=2, cex.lab=0.9, cex.axis=0.8, mfrow=c(1,2))
bronze.BW[1] <- 0.005
plot(bronze.BW, xlab="Number of clusters", ylab="BSS/TSS", type="n", bty="n", xaxt="n", xlim=c(1.5,10.5), ylim=c(0,0.7))
axis(1, at=2:10, labels=2:10, tick=FALSE)
for(g in 2:10) segments(g, 0, g, bronze.BW[g], lwd=10, col="gray", lend=1)
bronze.BWinc <- bronze.BW[2:10]-bronze.BW[1:9]  
plot(2:10, bronze.BWinc[1:9], type="n", xlab="Number of clusters",  bty="n", xaxt="n", xlim=c(2.5,12), 
     ylim=c(0,0.13), ylab="Improvement in BSS/TSS")
axis(1, at=3:10, labels=3:10, tick=FALSE)
for(g in 3:6) segments(g, 0, g, bronze.BWinc[g-1], lwd=10, col=rep("pink",5), lend=1)
for(g in 7:10) segments(g, 0, g, bronze.BWinc[g-1], lwd=10, col=rep("lightblue",3), lend=1)
### (looks like 6-cluster solution is a good choice)
set.seed(1649)  # set this seed if you want to get same solution as paper
                # otherwise clusters come out in different order
bronze.kmeans6 <- kmeans(bronze.CLR$LR, centers=6, nstart=20, iter.max=200)
bronze.kmeans6$betweenss/bronze.kmeans6$totss
# 0.5326932
### cluster sizes (these can be in a different order)
bronze.kmeans6$size
# [1]  38  14  61  61 120  68

### Table 2
table(bronze.kmeans6$cluster, BRONZE$GROUP2)
#     Eastern Zhou Shang Western Zhou
#   1            2     6           30
#   2            3     9            2
#   3           16    18           27
#   4            3    11           47
#   5           67    16           37
#   6            5    49           14

colSums(table(bronze.kmeans6$cluster, BRONZE$GROUP2))
# Eastern Zhou        Shang Western Zhou 
#           96          109          157
# correct assignments in the clusters (majority wins in each)
# 67/96 = 69.8% (EZhou)  58/109 = 53.2% (Shang)  104/157 = 66.2% (WZhou)	    

### Sum of maximum cluster assignments
(67+58+104)/362 
# [1] 0.6325967

### In order to see which samples are in which clusters, use the cluster number
### in $cluster of the kmeans object. For example, to see the 14 samples in clsuter 2:
# Samples in each cluster can be extracted, e.g. for cluster 2
BRONZE[bronze.kmeans6$cluster==2,]
#     Cat.No. Chronological_typology    GROUP       GROUP2     Cu     Sn     Pb     Zn    Au   Ag    As     Sb
# 8        27                  -1250   ShangA        Shang 787000 167000  26000 1670.0 0.325   81  1020   79.0
# 12       42                  -1250   ShangA        Shang 819000 188000   9500   23.3 0.270   61   130   80.0
# 13       43                  -1250   ShangA        Shang 806000  52000 120000    8.3 0.840  147   277   79.0
# etc...    :                      :        :            :      :      :      :     :   :       :     :     :
# 350      21                   -400 EastZhou Eastern Zhou 675000  99000 193000   94.0 0.340 2820  4100 1330.0
# 351      21                   -400 EastZhou Eastern Zhou 678000 100000 187000   80.0 0.380 2780  4200 1340.0


### ----------------------------------
### LRA as before, adding the clusters
### following code repeats earlier LRA
### Figure 7 plot "by hand" again: up to the adding of the clusters, the code is the same as before
par(mar=c(3.8,3.6,2.5,1), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.7, mfrow=c(1,1))
plot(1.05 * rbind(bronze.rpc, bronze.csc), type = "n", asp = 1,
     xlab = paste("LRA dimension ", 1, " (", round(percs[1], 1), "%)", sep = ""), 
     ylab = paste("LRA dimension ", 2, " (", round(percs[2], 1), "%)", sep = ""),
     main = "Logratio Analysis (with clusters)")
abline(h = 0, v = 0, col = "gray", lty = 2)
points(bronze.rpc, pch = group.pch[group.num], col = group.col[group.num], font = 1, 
       cex = group.cex[group.num])
text(bronze.csc, labels = colnames(bronze), col = "red", cex = 0.9, font = 4)    
legend("bottomleft", legend=c("S:Shang","WZ:WZou","EZ:EZhou"), bty="n",
       pch=group.pch, col=group.col, pt.cex=group.cex, 
       text.col=group.col, cex=0.9, text.font=2)
for (j in 1:3){
 points     <- bronze.rpc[group.num==j,1:2]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,], lty = 3, col=group.col[j])
}
set.seed(123)
CIplot_biv(bronze.rpc[,1], bronze.rpc[,2], group=group.num, groupcols=group.col, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           groupnames=c("S","WZ","EZ"))
set.seed(123)
CIplot_biv(bronze.rpc[,1], bronze.rpc[,2], group=group.num, groupcols=group.col, 
           alpha=0.99, add=TRUE, shade=FALSE, groupnames=c("S","WZ","EZ"))

### now add the clusters
set.seed(123)
CIplot_biv(bronze.rpc[,1], bronze.rpc[,2], group=bronze.kmeans6$cluster, groupcols=rep("gray5",6), 
           add=TRUE, shade=TRUE, alpha=0.99, alpha.f=0.1,
           groupnames=1:6)
set.seed(123)
CIplot_biv(bronze.rpc[,1], bronze.rpc[,2],group=bronze.kmeans6$cluster, groupcols=rep("gray40",6),
           alpha=0.99, add=TRUE, shade=FALSE, shownames=FALSE)
 
### check log(Pb/Zn), log(Zn/Sb) and log(Pb/Sb) means for the three clusters 1,3,4
PbZn <- log(bronze[,"Pb"]/bronze[,"Zn"])
aggregate(PbZn ~factor(bronze.kmeans6$cluster), FUN=mean)
#   factor(bronze.kmeans6$cluster)     PbZn
# 1                              1 2.808672
# 2                              2 7.203565
# 3                              3 5.655373
# 4                              4 5.829333
# 5                              5 8.190461
# 6                              6 6.945133
CuSb <- log(bronze[,"Cu"]/bronze[,"Sb"])
aggregate(CuSb ~factor(bronze.kmeans6$cluster), FUN=mean)
  factor(bronze.kmeans6$cluster)     CuSb
# 1                              1 9.251981
# 2                              2 8.290710
# 3                              3 7.040530
# 4                              4 7.299986
# 5                              5 6.470486
# 6                              6 8.993067
PbSb <- log(bronze[,"Pb"]/bronze[,"Sb"])
aggregate(PbSb ~factor(bronze.kmeans6$cluster), FUN=mean)
#   factor(bronze.kmeans6$cluster)     PbSb
# 1                              1 2.977139
# 2                              2 5.634079
# 3                              3 5.048424
# 4                              4 3.444391
# 5                              5 4.747746
# 6                              6 6.326600


### -------------------------------
### clustering the parts (Figure 8)
### -------------------------------
### A. Ward clustering of the columns of bronze.pro
bronze.ward <- WARD(CLR(CLOSE(t(bronze.pro)))$LR, weight=FALSE)
par(mar=c(1,3,3,1), mgp=c(2,0.7,0), font.lab=2, mfrow=c(1,2))
plot(bronze.ward, hang=-1, main="Ward clustering", labels=colnames(bronze.pro))
### B. amalgamation clustering
bronze.aclust <- ACLUST(bronze.pro, weight=FALSE, close=FALSE)
plot(bronze.aclust, hang=-1, main="Amalgamation clustering")
# (Total variance = 0.8021)

### ----------------------------
### find best ALR transformation
### ----------------------------
colnames(bronze.pro)
# [1] "Cu" "Sn" "Pb" "Zn" "Au" "Ag" "As" "Sb"
bronze.FINDALR <- FINDALR(bronze.pro, weight=FALSE)
bronze.FINDALR$procrust.cor
# [1] 0.9528238 0.9201740 0.8874291 0.8871596 0.8960438 0.9584928 0.9268800 0.9133485
# (Ag gives highest Procrustes, followed closely by Cu)
bronze.FINDALR$var.log
# [1] 0.008942823 0.358539951 2.460116665 1.516660452 1.370384677 0.485934010 1.252131422 1.778686302
# (Cu gives very small variance)

### plot the PCA of the selected ALRs
bronze.ALR <- ALR(bronze.pro, denom=1, weight=FALSE)$LR
bronze.ALR.PCA <- PCA(bronze.ALR, weight=FALSE)
### Row principal coordinates
bronze.ALR.rpc <- bronze.ALR.PCA$rowpcoord
### Column standard coordinates 
bronze.ALR.csc <- bronze.ALR.PCA$colcoord
bronze.ALR.rpc[,1] <- -bronze.ALR.rpc[,1]     # reverse first axis
bronze.ALR.csc[,1] <- -bronze.ALR.csc[,1]     # reverse first axis
percs.ALR.PCA  <- 100 * bronze.ALR.PCA$sv^2 / sum(bronze.ALR.PCA$sv^2)
sum(bronze.ALR.PCA$sv^2)
# [1] 1.397757

### Figure 9 plot "by hand"
par(mar=c(4.2,4,4.5,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
rescale <- 1.5
plot(1.05 * rbind(bronze.ALR.rpc, rescale * bronze.ALR.csc), type = "n", asp = 1,
     xlab = paste("PCA dimension ", 1, " (", round(percs.ALR.PCA[1], 1), "%)", sep = ""), 
     ylab = paste("PCA dimension ", 2, " (", round(percs.ALR.PCA[2], 1), "%)", sep = ""),
     xaxt = "n", yaxt = "n", main="PCA of ALRs with reference Cu")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = "red", col.axis = "red")
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
     col = "black", col.ticks = "red", col.axis = "red")
arrows(0,0,0.95*rescale*bronze.ALR.csc[,1],0.95*rescale*bronze.ALR.csc[,2],angle=10, length=0.1,col="pink", lwd=2)
points(bronze.ALR.rpc, pch = group.pch[group.num], col = group.col[group.num], font = 1, cex = group.cex[group.num])
text(rescale * bronze.ALR.csc, labels = colnames(bronze.ALR), col = "red", cex = 0.9, font = 4)    
legend("bottomleft", legend=c("S:Shang","WZ:WZou","EZ:EZhou"), bty="n",
       pch=group.pch, col=group.col, pt.cex=group.cex, 
       text.col=group.col, cex=0.9, text.font=2)
for (j in 1:3){
 points     <- bronze.ALR.rpc[group.num==j,1:2]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,], lty = 3, col=group.col[j])
}
require(ellipse)
set.seed(123)
CIplot_biv(bronze.ALR.rpc[,1], bronze.ALR.rpc[,2], group=group.num, groupcols=group.col, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           groupnames=c("S","WZ","EZ"), shownames=TRUE)
set.seed(123)
CIplot_biv(bronze.ALR.rpc[,1], bronze.ALR.rpc[,2], group=group.num, groupcols=group.col, 
           alpha=0.99, add=TRUE, shade=FALSE, groupnames=c("S","WZ","EZ"), shownames=TRUE)


### -------------------
### SUPERVISED LEARNING
### -------------------
### predicting chronology
chrono <- BRONZE[,2]
chrono.STEPR <- STEPR(bronze.pro, chrono, method=1)
# [1] "Criterion increases when 3-th ratio enters"
chrono.STEPR$ratios
      row col
Cu/Sb   1   8
Zn/Au   4   5

### linear model with two predictors
chrono.lm <- lm(chrono ~ ., data=as.data.frame(chrono.STEPR$logratios))
summary(chrono.lm)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -335.504     64.519  -5.200 3.35e-07 ***
# `Cu/Sb`      -71.258      8.569  -8.316 1.91e-15 ***
# `Zn/Au`      -35.687      7.427  -4.805 2.27e-06 ***
# ---
# Residual standard error: 220.7 on 359 degrees of freedom
# Multiple R-squared:  0.2381,    Adjusted R-squared:  0.2339 
# F-statistic: 56.11 on 2 and 359 DF,  p-value: < 2.2e-16

### predictions of dynasty groups from predicted chronologies
chrono.lm.pred <- predict(chrono.lm)
# recall:
#   Shang        (c. 1500-1046 BCE) 
#   Western Zhou (c. 1046–771 BCE)
#   Eastern Zhou (c. 771–256 BCE)
chrono.lm.pred.dyn <- rep("WZhou", 362)
chrono.lm.pred.dyn[chrono.lm.pred < (-1046)] <- "Shang" 
chrono.lm.pred.dyn[chrono.lm.pred > (-771)]  <- "EZhou" 
table(chrono.lm.pred.dyn, BRONZE$GROUP2)[c(2,3,1),c(2,3,1)]
# chrono.lm.pred.dyn Shang Western Zhou Eastern Zhou
#              Shang    38           22            1
#              WZhou    66          129           75
#              EZhou     5            6           20
(38+129+20)/362
# [1] 0.5165746

### regression tree
bronze.ratios <- as.data.frame(exp(LR(bronze.pro, weight=FALSE)$LR))
bronze.rtree  <- rpart(chrono ~ ., data=bronze.ratios)
prp(bronze.rtree, extra=1, font=2)
bronze.rtree.pred <- predict(bronze.rtree)
# recall:
#   Shang        (c. 1500-1046 BCE) 
#   Western Zhou (c. 1046–771 BCE)
#   Eastern Zhou (c. 771–256 BCE)
bronze.rtree.pred.dyn <- rep("WZhou", 362)
bronze.rtree.pred.dyn[bronze.rtree.pred < (-1046)] <- "Shang" 
bronze.rtree.pred.dyn[bronze.rtree.pred > (-771)]  <- "EZhou" 
table(bronze.rtree.pred.dyn, BRONZE$GROUP2)[c(2,3,1),c(2,3,1)]
#  bronze.rtree.pred.dyn Shang Western Zhou Eastern Zhou
#                 Shang    69           29            0
#                 WZhou    38          124           28
#                 EZhou     2            4           68
(69+124+68)/362
# [1] 0.7209945


### ----------------------------------------------------------
### Multinomial logistic regression predicting all five phases
### Uses package 'nnet'. Stepwise procedure followed step-by-step
### Uses Bonferroni as stopping criterion with nmumber of 
### parameters m = 2 x the estimated parameters in one equation
### Bonferroni = (AIC - 2m = -2 log Likelihood) + 9.23 * m
bronze.LR <- LR(bronze.pro, weight=FALSE)
bronze.groups <- group
multlogit <- rep(999, 28)
names(multlogit) <- colnames(bronze.LR$LR) 
for(jk in 1:28){
  test <- multinom(bronze.groups ~ bronze.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Cu/Pb 
#     2
step1 <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"])
summary(step1)
# AIC: 694.3577  
# Bonferroni
694.3577 - 2*4 + 9.23*4
# [1] 723.2777
# Success of prediction
sum(diag(table(predict(step1), bronze.groups)))/362
# [1] 0.5745856

# step2
multlogit <- rep(999, 28)
names(multlogit) <- colnames(bronze.LR$LR) 
for(jk in (1:28)[-2]){
  test <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"] + bronze.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Cu/Sb 
#    7
step2 <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"]+ bronze.LR$LR[,"Cu/Sb"])
summary(step2)
# AIC: 626.3514 
# Bonferroni
626.3514 - 2*6 + 9.23*6
# [1] 669.7314
# Success of prediction
sum(diag(table(predict(step2), bronze.groups)))/362
# [1] 0.6685083

# step3
multlogit <- rep(999, 28)
names(multlogit) <- colnames(bronze.LR$LR) 
for(jk in (1:28)[-c(2,7)]){
  test <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"] + bronze.LR$LR[,"Cu/Sb"] + bronze.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Zn/Ag
#    20 
step3 <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"]+ bronze.LR$LR[,"Cu/Sb"]+ bronze.LR$LR[,"Zn/Ag"])
summary(step3)
# AIC: 602.6249
# Bonferroni
602.6249 - 2*8 + 9.23*8
# [1] 660.4649
# Success of prediction
sum(diag(table(predict(step3), bronze.groups)))/362
# [1] 0.6740331

# step4
multlogit <- rep(999, 28)
names(multlogit) <- colnames(bronze.LR$LR) 
for(jk in (1:28)[-c(2,7,20)]) {
  test <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"]+ bronze.LR$LR[,"Cu/Sb"]+ bronze.LR$LR[,"Zn/Ag"]+ bronze.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Au/Ag
#    23 
step4 <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"]+ bronze.LR$LR[,"Cu/Sb"]+ bronze.LR$LR[,"Zn/Ag"]+ bronze.LR$LR[,"Au/Ag"])
summary(step4)
# AIC: 584.0358 
584.0358 - 2*10 + 9.23*10
# Bonferroni
# [1] 656.3358
# Success of prediction
sum(diag(table(predict(step4), bronze.groups)))/362
# [1] 0.698895

# step5
multlogit <- rep(999, 28)
names(multlogit) <- colnames(bronze.LR$LR) 
for(jk in (1:28)[-c(2,7,20,23)]){
  test <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"]+ bronze.LR$LR[,"Cu/Sb"]+ bronze.LR$LR[,"Zn/Ag"]+ bronze.LR$LR[,"Au/Ag"] + bronze.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Zn/Au
#    19 
step5 <- multinom(bronze.groups ~ bronze.LR$LR[,"Cu/Pb"]+ bronze.LR$LR[,"Cu/Sb"]+ bronze.LR$LR[,"Zn/Ag"]+ bronze.LR$LR[,"Au/Ag"]+ bronze.LR$LR[,"Zn/Au"])
summary(step5)
# AIC: 584.0358    ### almost same as before!
584.0358 - 2*12 + 9.23*12
# Bonferroni
# [1] 670.7958
# Success of prediction
sum(diag(table(predict(step5), bronze.groups)))/362
# [1] 0.698895   ### no change, must be a cycle!

### hence, stop with four terms
sum(diag(table(predict(step4), bronze.groups)))/362
# [1] 0.698895
table(predict(step4), bronze.groups)
#                Shang Western Zhou Eastern Zhou
#   Shang           68           22           10
#   Western Zhou    26          116           17
#   Eastern Zhou    15           19           69
#             68/109=62%   116/157=74%    69/96=72%
# (68+116+69)/362 = 69.9%

summary(step4)
# Coefficients:
#              (Intercept) bronze.LR$LR[, "Cu/Pb"] bronze.LR$LR[, "Cu/Sb"] bronze.LR$LR[, "Zn/Ag"] bronze.LR$LR[, "Au/Ag"]
# Western Zhou    8.460987               0.6864055              -0.9005493             -0.01582178               0.8566738
# Eastern Zhou    7.340674              -0.2572281              -0.9726675             -0.61740113               0.3670683

# Std. Errors:
#              (Intercept) bronze.LR$LR[, "Cu/Pb"] bronze.LR$LR[, "Cu/Sb"] bronze.LR$LR[, "Zn/Ag"] bronze.LR$LR[, "Au/Ag"]
# Western Zhou    1.312340               0.1255817               0.1484844               0.1144972               0.1933090
# Eastern Zhou    1.456138               0.1941253               0.1758611               0.1516167               0.1736836

# Residual Deviance: 564.0358 
# AIC: 584.0358 

### coefficients of log-contrast
beta.Cu <-  summary(step4)$coefficients[,2] + summary(step4)$coefficients[,3]
beta.Pb <- -summary(step4)$coefficients[,2]
beta.Sb <- -summary(step4)$coefficients[,3]
beta.Zn <-  summary(step4)$coefficients[,4]
beta.Ag <- -summary(step4)$coefficients[,4] - summary(step4)$coefficients[,5] 
beta.Au <-  summary(step4)$coefficients[,5]
betas <- cbind(beta.Cu, beta.Pb, beta.Sb, beta.Zn, beta.Ag, beta.Au)
betas
#                 beta.Cu    beta.Pb   beta.Sb     beta.Zn    beta.Ag   beta.Au
# Western Zhou -0.2141438 -0.6864055 0.9005493 -0.01582178 -0.8408521 0.8566738
# Eastern Zhou -1.2298956  0.2572281 0.9726675 -0.61740113  0.2503328 0.3670683

# check that coefficients of log-contrast sum to 0
rowSums(betas)
# Western Zhou Eastern Zhou 
# -3.122502e-17  1.110223e-16 

### perform bootstrapping on the log-contrast coefficients (output for each bootstrap)
nboot <- 1000
coefs.boot <- array(0, c(2,6,nboot))
set.seed(123)
for(iboot in 1:nboot) {
  sample.boot <- sample(1:nrow(bronze.pro), replace=TRUE)
  bronze.boot <- bronze.LR$LR[sample.boot,]
  groups.boot <- bronze.groups[sample.boot] 
  step.boot <- multinom(groups.boot ~ bronze.boot[,"Cu/Pb"]+ bronze.boot[,"Cu/Sb"]+ bronze.boot[,"Zn/Ag"]+ bronze.boot[,"Au/Ag"])
  coefs.boot[,1,iboot] <-  summary(step.boot)$coefficients[,2] + summary(step.boot)$coefficients[,3]
  coefs.boot[,2,iboot] <- -summary(step.boot)$coefficients[,2]
  coefs.boot[,3,iboot] <- -summary(step.boot)$coefficients[,3]
  coefs.boot[,4,iboot] <-  summary(step.boot)$coefficients[,4]
  coefs.boot[,5,iboot] <- -summary(step.boot)$coefficients[,4] - summary(step.boot)$coefficients[,5]
  coefs.boot[,6,iboot] <-  summary(step.boot)$coefficients[,5]
}

### bootstrap confidence intervals, for WZhou(W) vs. Shang and EZhou (E) vs. Shang
intervalsW <- apply(coefs.boot[1,,], 1, quantile, c(0.01,0.99))
colnames(intervalsW) <- c("Cu","Pb","Sb","Zn","Ag","Au")
intervalsW
#              Cu        Pb        Sb         Zn         Ag        Au
# 1%  -0.61763891 -1.040585 0.5497325 -0.2384277 -1.6511015 0.4013546
# 99%  0.09880202 -0.451761 1.4078563  0.2699125 -0.3434184 1.6159717
intervalsE <- apply(coefs.boot[2,,], 1, quantile, c(0.01,0.99))
colnames(intervalsE) <- c("Cu","Pb","Sb","Zn","Ag","Au")
intervalsE
#             Cu         Pb        Sb         Zn         Ag        Au
# 1%  -2.0872587 -0.2462838 0.5602343 -1.2084723 -0.6845684 0.0583307
# 99% -0.7339686  1.0892344 1.5942542 -0.2203947  0.8634290 1.3296408

### plot bootstrap confidence intervals
par(mar=c(4.5,3,2.5,2), font.lab=2, mgp=c(2,0.7,0), cex.lab=1, cex.axis=0.7, font.axis=2)
plot(0,0,type="n",xlab="Effect size on log-odds", ylab="",bty="n", main="Log-contrast effect sizes for WZhou & EZhou vs. Shang",
     xlim=c(-3.5,3),ylim=c(0.5,12.5), yaxt="n")
axis(2, at=13-seq(1.6,11.6,2), labels=colnames(intervalsW), tick=FALSE, las=1, cex.axis=1)
for(j in seq(1.5,11.5,2)) abline(h=j-0.2, lwd=45, col="gray95")
text(rep(-3.35,12), seq(0.8,11.8,1), labels=rep(c("WZhou","EZhou"), 6), col=rep(c("blue","red"),6), cex=0.8)
segments(0,0.5,0,12.3, lwd=2, col="gray", lty=2)

for(j in 1:6) {
  segments(intervalsW[1,7-j],2*j-1.2,intervalsW[2,7-j],2*j-1.2, lwd=2, col="blue", lend=2)
}
eps <- 0.1
for(j in 1:6) {
  segments(intervalsW[1,7-j],2*j-1.2-eps,intervalsW[1,7-j],2*j-1.2+eps, lwd=2, col="blue", lend=2)
  segments(intervalsW[2,7-j],2*j-1.2-eps,intervalsW[2,7-j],2*j-1.2+eps, lwd=2, col="blue", lend=2)
}
points(betas[1,6:1],seq(1,11,2)-0.2, pch=21, col="blue", bg="white", cex=1)

for(j in 1:6) {
  segments(intervalsE[1,7-j],2*j-0.2,intervalsE[2,7-j],2*j-0.2, lwd=2, col="red", lend=2)
}
eps <- 0.1
for(j in 1:6) {
  segments(intervalsE[1,7-j],2*j-0.2-eps,intervalsE[1,7-j],2*j-0.2+eps, lwd=2, col="red", lend=2)
  segments(intervalsE[2,7-j],2*j-0.2-eps,intervalsE[2,7-j],2*j-0.2+eps, lwd=2, col="red", lend=2)
}
points(betas[2,6:1],seq(2,12,2)-0.2, pch=21, col="red", bg="white", cex=1)

### ----------------------------------------
### classification tree predicting dynasties
bronze.short <- rep("S", 362)
bronze.short[BRONZE[,4]=="Western Zhou"] <- "WZ"
bronze.short[BRONZE[,4]=="Eastern Zhou"] <- "EZ"
bronze.ratios <- as.data.frame(exp(LR(bronze.pro, weight=FALSE)$LR))
bronze.ctree <- rpart(factor(bronze.short) ~ ., data=bronze.ratios)
par(mar=c(.1,.1,.1,.1))
require(rpart.plot)
rpart.plot(bronze.ctree, cex=0.7, font=2, extra=104)#, box.palette=gp2.col[c(2,3,1)])
prp(bronze.ctree, extra=1) 

### Predictions of the groups
foo.pred <- predict(bronze.ctree)
foo.pred.max <- apply(foo.pred, 1, max)
bronze.ctree.pred <- rep(0, nrow(bronze))
for(i in 1:nrow(bronze)) bronze.ctree.pred[i] <- which(foo.pred[i,]==foo.pred.max[i])
table(bronze.ctree.pred, bronze.short)[c(2,3,1),c(2,3,1)]
#                  bronze.short
# bronze.ctree.pred   S  WZ  EZ
#                 2  79  13   5
#                 3  25 133  22
#                 1   5  11  69
sum(diag(table(bronze.ctree.pred, BRONZE[,4])))/nrow(bronze)
# [1] 0.7762431

# correct predictions in each group
apply(table(bronze.ctree.pred, bronze.short),2,max)/colSums(table(bronze.ctree.pred, bronze.short))
#        EZ         S        WZ 
# 0.7187500 0.7247706 0.8471338 
# (e.g., for EZ: 69/96 = 0.71875)

### ---------------------------------------------------------------------------
### cross-validation of regression tree followed by classification into periods
### define stratified training and test sets
foo  <- data.frame(list(ratios=bronze.ratios, periods=bronze.short))
rownames(foo) <- 1:nrow(foo)
set.seed(1539)    # use his seed in the cross-validaions to agree with article
                  # different results can be obtained depending on R version
test <- stratified(foo, "periods", 0.3333333, keep.rownames=TRUE)
# indices of test set in first column $rn
test.index <- as.numeric(test$rn)
train.index <- (1:nrow(foo))[-test.index]
test  <- foo[test.index,]
train <- foo[train.index,]
# remove last columns of train and test, containing the periods
train <- train[,-ncol(train)]
test  <- test[,-ncol(test)]
train.y <- chrono[train.index]
test.y  <- chrono[test.index]

### cross-validation of regression + classification
set.seed(1539)
# we need 242 random allocations
groups10 <- c(rep(1:10, each=24),1:2) 
groups10 <- groups10[sample(1:242)]
cvpreds <- rep(0, length(train.y))
for(cv in 1:10) {
  cvset <- (1:nrow(train))[-which(groups10==cv)]
  foo <- rpart(train.y[cvset]~., data=train[cvset,])
  cvset <- which(groups10==cv)
  foo.pred <- predict(foo, newdata=train[cvset,])
  cvpreds[cvset] <- foo.pred
}

# recall
# Shang (c. 1500-1046 BCE), 
# Western Zhou (c. 1046–771 BCE)
# Eastern Zhou (c. 771–256 BCE)
cvpred.dyn <- rep("WZhou", length(train.y))
cvpred.dyn[cvpreds < (-1046)] <- "Shang" 
cvpred.dyn[cvpreds > (-771)]  <- "EZhou" 
table(cvpred.dyn, bronze.short[train.index])[c(2,3,1),c(2,3,1)]
# cvpred.dyn  S WZ EZ
#      Shang 33 18  2
#      WZhou 32 74 21
#      EZhou  8 13 41
(33+74+41)/242
# [1] 0.6115702

# applied to test set
foo <- rpart(train.y ~., data = train)
foo.pred <- predict(foo, newdata=test)
cvpred.dyn <- rep("WZhou", nrow(test))
cvpred.dyn[foo.pred < (-1046)] <- "Shang" 
cvpred.dyn[foo.pred > (-771)]  <- "EZhou" 
table(cvpred.dyn, bronze.short[test.index])[c(2,3,1),c(2,3,1)]
# cvpred.dyn  S WZ EZ
#      Shang 13  7  3
#      WZhou 15 42 16
#      EZhou  8  3 13
(13+42+13)/120
# [1] 0.5666667

### ---------------------------------------------------------------------------
### cross-validation on multinomial logit
### define stratified training and test sets
foo  <- data.frame(list(logratios=bronze.LR$LR, periods=bronze.short))
rownames(foo) <- 1:nrow(foo)
set.seed(1539)
test <- stratified(foo, "periods", 0.3333333, keep.rownames=TRUE)
# indices of test set in first column $rn
test.index <- as.numeric(test$rn)
train.index <- (1:nrow(foo))[-test.index]
test  <- foo[test.index,]
train <- foo[train.index,]
# remove last columns of train and test, containing the periods
train <- train[,-ncol(train)]
test  <- test[,-ncol(test)]
train.y <- bronze.short[train.index]
test.y  <- bronze.short[test.index]

### cross-validation on multinomial logit
# we need 242 random allocations
set.seed(1539)
groups10 <- c(rep(1:10, each=24),1:2) 
groups10 <- groups10[sample(1:242)]
cvpreds <- matrix(0, length(train.y), 3)
cvpreds <- rep(0, length(train.y))
for(cv in 1:10) {
  cvset <- (1:246)[-which(groups10==cv)]
  foo <- multinom(train.y[cvset] ~ logratios.Cu.Pb+logratios.Cu.Sb+logratios.Zn.Ag+logratios.Au.Ag, data=train[cvset,])
  cvset <- which(groups10==cv)
  foo.pred <- predict(foo, newdata=train[cvset,]) 
  cvpreds[cvset] <- foo.pred
}
table(cvpreds, train.y)[c(2,3,1), c(2,3,1)]
#         train.y
# cvpreds  S WZ EZ
#       2 45 19  6
#       3 19 71 13
#       1  9 15 45
(45+71+45)/242
# [1] 0.6652893

# apply to test data
foo <- multinom(train.y ~ logratios.Cu.Pb+logratios.Cu.Sb+logratios.Zn.Ag+logratios.Au.Ag, data=train)
test.pred <- predict(foo, newdata=as.data.frame(test))
table(test.pred, bronze.short[test.index])[c(2,3,1), c(2,3,1)]
# test.pred  S WZ EZ
#        S  19  5  4
#        WZ  8 40  8
#        EZ  9  7 20
(19+40+20)/120
# [1] 0.6583333


### --------------------------------------------------------------------
### cross-validation on classification tree
### define stratified training and test sets (this is a repeat of above)
foo  <- data.frame(list(ratios=bronze.ratios, periods=bronze.short))
rownames(foo) <- 1:nrow(foo)
set.seed(1539)
test <- stratified(foo, "periods", 0.3333333, keep.rownames=TRUE)
# indices of test set in first column $rn
test.index <- as.numeric(test$rn)
train.index <- (1:nrow(foo))[-test.index]
test  <- foo[test.index,]
train <- foo[train.index,]
# remove last columns of train and test, containing the periods
train <- train[,-ncol(train)]
test  <- test[,-ncol(test)]
train.y <- bronze.short[train.index]
test.y  <- bronze.short[test.index]

# we need 242 random allocations
set.seed(1539)
groups10 <- c(rep(1:10, each=24),1:2) 
groups10 <- groups10[sample(1:242)]
cvpreds <- matrix(0, length(train.y), 3)
for(cv in 1:10) {
  cvset <- (1:nrow(train))[-which(groups10==cv)]
  foo <- rpart(factor(train.y[cvset]) ~ ., data=train[cvset,])
  cvset <- which(groups10==cv)
  foo.pred <- predict(foo, newdata=train[cvset,])
  cvpreds[cvset,] <- foo.pred
}
cv.pred.period <- rep(0, length(train.y))
for(i in 1:length(train.y)) {
  tmp <- which(cvpreds[i,] == max(cvpreds[i,]))
  if(length(tmp)>1) cv.pred.period[i] <- tmp[sample(tmp,1)]
  if(length(tmp) == 1) cv.pred.period[i] <- tmp
}
table(cv.pred.period, train.y)[c(2,3,1),c(2,3,1)]
#               train.y
# cv.pred.period  S WZ EZ
#              2 49 15  5
#              3 17 74 13
#              1  7 16 46
(49+74+46)/242
# [1] 0.6983471

# applied to test set
foo <- rpart(factor(train.y) ~ ., data=train)
foo.pred <- predict(foo, newdata=test)
test.pred <- rep(0, length(test.index))
for(i in 1:length(test.index)) test.pred[i] <- which(foo.pred[i,] == max(foo.pred[i,]))
table(test.pred, bronze.short[test.index])[c(2,3,1),c(2,3,1)]
# test.pred  S WZ EZ
#         2 13  4  2
#         3 14 41  9
#         1  9  7 21
(13+41+21)/120
# [1] 0.625

### Random forest classification
foo  <- data.frame(list(ratios=bronze.ratios, periods=bronze.short))
rownames(foo) <- 1:nrow(foo)
# use same training and test sets as before
set.seed(1549)
test <- stratified(foo, "periods", 0.3333333, keep.rownames=TRUE)
# indices of test set in first column $rn
test.index <- as.numeric(test$rn)
train.index <- (1:nrow(foo))[-test.index]
test  <- foo[test.index,]
train <- foo[train.index,]
# remove last columns of train and test, containing the periods
train <- train[,-ncol(train)]
test  <- test[,-ncol(test)]
train.y <- bronze.short[train.index]
test.y  <- bronze.short[test.index]
test.x  <- data.frame(test)
train.x <- data.frame(train)
set.seed(123)
bronze.rf <- randomForest(x=train.x,
                       y=factor(train.y),
                       xtest=test.x,
                       ytest=factor(test.y),     
                       ntree=1000, do.trace=FALSE, importance=TRUE)
# (the do.trace=T gives results for each of the 1500 trees).
bronze.rf
# No. of variables tried at each split: 5
# 
#         OOB estimate of  error rate: 25.21%
# Confusion matrix:
#    EZ  S WZ class.error
# EZ 47  4 13   0.2656250
# S   6 47 20   0.3561644
# WZ  8 10 87   0.1714286
#                 Test set error rate: 21.67%
# Confusion matrix:
#    EZ  S WZ class.error
# EZ 24  5  3   0.2500000
# S   2 26  8   0.2777778
# WZ  3  5 44   0.1538462

# in the above "true" dynasties are the rows, predictions in columns
# also note that the values have to be re-ordered to obtain Table 8
# the tranposed and re-ordered tables can be obtained as follows:

# our OOB confusion matrix, transposed version of above result in bronze$rf
table(bronze.rf$predicted, train.y)[c(2,3,1), c(2,3,1)]
#      train.y
#       S WZ EZ
#   S  47 10  4
#   WZ 20 87 13
#   EZ  6  8 47

# our test set confusion matrix, transposed version of abve result in bronze.rf
table(bronze.rf$test$predicted, test.y)[c(2,3,1), c(2,3,1)]
#     test.y
#       S WZ EZ
#   S  26  5  5
#   WZ  8 44  3
#   EZ  2  3 24




