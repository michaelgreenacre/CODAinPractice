### ------------------- KIMBERLITE DATA ------------------------
### First set the working directory to where the data file is
### using setwd() function or pull-down menu
	
### Read the data
kim270 <- read.table("kimberlite.cation.txt", header=T)
dim(kim270)
# [1] 270  24
colnames(kim270)
# [1] "SampleID"  "StratUnit" "Si"        "Ti"        "Al"        "Fe"       
# [7] "Mg"        "Ca"        "Na"        "K"         "P"         "Rb"       
# [13] "Nb"        "Zr"        "Th"        "V"         "Cr"        "Co"       
# [19] "Ni"        "La"        "Er"        "Yb"        "Y"         "Ga" 

### The compositional data are in the 22 columns from 3 to 24
### and the rows sum to 100%
### The column "StratUnit" contains the names of the 5 phases

### The elements in decreasing order of average percetnage
order.kim <- order(colMeans(kim270[,3:24]), decreasing=TRUE)
round(colMeans(kim270[,3:24])[order.kim], 4)
#      Si      Mg      Fe      Ca      Al      Ti       P      Na       K      Ni      Cr 
# 42.3344 36.9554 10.0453  5.1805  2.6756  1.5662  0.3837  0.2858  0.2245  0.1660  0.1309 
#      Co      Nb       V      La      Zr      Th      Rb       Y      Ga      Er      Yb 
#  0.0112  0.0110  0.0087  0.0084  0.0081  0.0015  0.0014  0.0009  0.0005  0.0001  0.0001 

### Counts of the five phases and convert them to numeric
table(kim270$StratUnit)
# Cantuar     eJF     lJF     mJF   Pense 
#      21     154      28      40      27 
kim.su <- as.numeric(as.factor(kim270$StratUnit))
table(kim.su)    
#  1   2   3   4   5 
# 21 154  28  40  27 

### Present order of phases in file is eJF lJF mJF Pense Cantuar 
### ending at sample numbers 154 182 222 249 270 respectively
### The time order we want is Cantuar Pense eJF mJF lJF
### The 5x2 array su.ends gives the starting and ending sample numbers
### of the time-ordered five phases in the data file
phase.names <- c("Cantuar","Pense","eJF","mJF","lJF")
su.ends <- matrix(c(155,182, 183,222, 1,154, 223,249, 250,270), nrow=5, byrow=TRUE)
colnames(su.ends) <- c("start","end")
rownames(su.ends) <- rev(phase.names)  
### (for bar plot, most recent (lJF) at the top to earliest (Cantuar) at the bottom
su.ends
#         start end
# lJF       155 182
# mJF       183 222
# eJF         1 154
# Pense     223 249
# Cantuar   250 270

### Extract the compositional data and express percentages as proportions
kim <- kim270[,3:24] / 100
I <- nrow(kim)
J <- ncol(kim)
### (You can check that the rows all sum to 1)

### plot means and variances on log scales
kim.mean <- apply(kim, 2, mean)
kim.var  <- apply(kim, 2, var)
plot(kim.mean, kim.var, type="n", log="xy", xlab="Means (log scale)", ylab="Variances (log scale)")
text(kim.mean, kim.var, labels=colnames(kim))

### Centred logratios, and their variances similarly plotted against the means
### Needs package easyCODA
require(easyCODA)
kim.CLR <- CLR(kim, weight=FALSE)$LR
kim.CLR.var <- apply(kim.CLR, 2, var)
plot(kim.mean, kim.CLR.var, type="n", log="xy", xlab="Means (log scale)", ylab="CLR variances (log scale)")
text(kim.mean, kim.CLR.var, labels=colnames(kim))
### (This plot shows that Na, K and Rb will be amongst the most important elements in 
###  the unsupervised multivariate analyses, since they contribute high logratio variance)
  
### -----------------------------------------------------------------
### Plotting raw values after monotonic transformation of proportions
### Use Dark2 palette in RColorBrewer for the elements -- see Figure 2
require(RColorBrewer)
elements.col <- c(brewer.pal(8,"Dark2"),brewer.pal(8,"Dark2"),brewer.pal(6,"Dark2")) 
elements.order <- rev(order.kim) 
### Perform monotonic transformation from kim to kim.up.log (closed)
kim.up <- kim*(1/min(kim))
kim.up.log <- CLOSE(log(kim.up))
### Plot (pull out window so it is wide and flat, adjusting accordingly)
lab.start <- 1
par(mar=c(4,1,2,1), mgp=c(2,0.7,0), font.lab=2, cex.lab=0.9, las=1)
plot(0,0,xlab="Percentages",ylab="",xlim=c(-0.02,1),ylim=c(-25,274), xaxt="n", yaxt="n", bty="n", main="Kimberlite data", type="n")
nlines <- I+4
i.extra <- 0
ii <- 1
for(iphase in 1:5) {
  for(i in su.ends[iphase,1]:su.ends[iphase,2]) {
    foo <- cumsum(c(0,kim.up.log[i,elements.order]))
    for(j in 1:J) segments(foo[j], nlines-ii-i.extra, foo[j+1], nlines-ii-i.extra, col=elements.col[j],lwd=3,lend=1)
    ii <- ii+1
  }
  i.extra <- i.extra+1
  segments(0, nlines-ii-i.extra, 1, nlines-ii-i.extra, col="white",lwd=4)
}
foo <- cumsum(c(0,colMeans(kim.up.log[,elements.order])))
for(j in 1:J) segments(foo[j], -12, foo[j+1], -12, col=elements.col[j],lwd=10, lend=1)
text(c(0,0,0,0,0),c(10,35,121,224,259) , labels=phase.names, pos=2, font=2, cex=0.8)
text(0, -12, labels="Average", cex=0.8, font=4, pos=2)
labels.mid <- rep(0,J)
kim.mean.cum <- cumsum(colMeans(kim.up.log[,elements.order]))
kim.mean.cum.actual <- cumsum(colMeans(kim[,elements.order])) 
labels.mid[1] <- kim.mean.cum[1]/2
for(j in 2:J) labels.mid[j] <- kim.mean.cum[j-1]+(kim.mean.cum[j]-kim.mean.cum[j-1])/2
text(labels.mid[lab.start:22],rep(-23,6),labels=colnames(kim)[elements.order][lab.start:22],col=elements.col[lab.start:22], cex=0.8, font=4)
ticks <- c(0,0.00001,0.0001,0.001,0.01,0.1,0.2,0.5,1)
ticks <- ticks[2:length(ticks)]-ticks[1:(length(ticks)-1)]
sum(ticks)
cumsum(ticks)
ticks.up <- ticks*(1/min(kim))
ticks.up.log <- log(ticks.up)/sum(log(ticks.up))
ticks.up.cum <- cumsum(ticks.up.log)
axis(1, at=c(0,ticks.up.cum), labels=c("0","0.001","0.01","0.1","1","10","20","50","100"), cex.axis=0.7, font=2)

### ------------------------------------------------------------
### Two selected bivariate scatterplots in new window (Figure 3)
kim.su.col <- c("blue", "cyan","red","yellowgreen","black")
kim.su.pch  <- c(21,20,22,24,4)
kim.su.bg <- c("white","cyan","white","white","white")
kim.su.cex <- c(0.8, 1, 0.7, 0.6, 0.7)  
kim.su.order <- c(1,5,2,4,3)
kim.pair <- 100 * cbind(kim[,"Yb"], kim[,"P"]) ; colnames(kim.pair) <- c("Yb","P")
par(mar=c(4.2,4,2,1), mgp=c(2,0.7,0), font.lab=2, las=1, cex.axis=0.8)
plot(kim.pair[,1], kim.pair[,2], xlim=c(0,1.05*max(kim.pair[,1])), ylim=c(0, 1.05*max(kim.pair[,2])), 
     xlab=colnames(kim.pair)[1], ylab=colnames(kim.pair)[2], type="n")
points(kim.pair[,1], kim.pair[,2], col=kim.su.col[kim.su], pch=kim.su.pch[kim.su], bg=kim.su.bg[kim.su],
       cex=kim.su.cex[kim.su])
legend("topleft", legend=phase.names, text.col=c("blue", "cyan3","red","forestgreen","black")[kim.su.order], 
       pch=kim.su.pch[kim.su.order], col=kim.su.col[kim.su.order], pt.bg=kim.su.bg[kim.su.order],
       pt.cex=kim.su.cex[kim.su.order], cex=0.9, text.font=2, bty="n")

kim.pair <- 100 * cbind(kim[,"Zr"], kim[,"La"]) ; colnames(kim.pair) <- c("Zr","La")
par(mar=c(4.2,4,2,1), mgp=c(2,0.7,0), font.lab=2, las=1, cex.axis=0.8)
plot(kim.pair[,1], kim.pair[,2], xlim=c(0,1.05*max(kim.pair[,1])), ylim=c(0, 1.05*max(kim.pair[,2])), 
     xlab=colnames(kim.pair)[1], ylab=colnames(kim.pair)[2], type="n")
points(kim.pair[,1], kim.pair[,2], col=kim.su.col[kim.su], pch=kim.su.pch[kim.su], bg=kim.su.bg[kim.su],
       cex=kim.su.cex[kim.su])
legend("topleft", legend=phase.names, text.col=c("blue", "cyan3","red","forestgreen","black")[kim.su.order], 
       pch=kim.su.pch[kim.su.order], col=kim.su.col[kim.su.order], pt.bg=kim.su.bg[kim.su.order],
       pt.cex=kim.su.cex[kim.su.order], cex=0.9, text.font=2, bty="n")

### -------------------------------------------
### Compute contributed and explained variances
### Total variance computed using pairwise logratios and CLRs
### Both ways givs identical result
### All J*(J-1)/2 pairwise logratios
require(easyCODA)
kim.LR <- LR(kim, weight=FALSE)
kim.LR.var <- LR.VAR(kim.LR, vars=TRUE)
kim.LR.var$LRtotvar
# [1] 0.1132472
sum(kim.LR.var$LRvars)
# [1] 0.1132472
### The J centered logratios (CLRs)
kim.CLR <- CLR(kim, weight=FALSE)
kim.CLR.var <- LR.VAR(kim.CLR, vars=TRUE)
kim.CLR.var$LRtotvar
# [1] 0.1132472
sum(kim.CLR.var$LRvars)
# [1] 0.1132472
### Note that the sum of squared singular values in the LRA is also the total variance
sum(LRA(kim, weight=FALSE)$sv^2)
# [1] 0.1132472

### Percentage contributions of LRs and CLRs to total variance
kim.LR.ctr <- kim.LR.var$LRvars / sum(kim.LR.var$LRvars) 
kim.LR.ctr.order <- order(kim.LR.ctr, decreasing=TRUE)
head(round(100*kim.LR.ctr[kim.LR.ctr.order],2), 10)
# Table 1, column (i)
# Na/Ni Mg/Na Na/Co Ca/Na  Na/P Na/La Na/Nb Fe/Na Ti/Na Na/Th 
#  1.78  1.76  1.74  1.69  1.66  1.66  1.61  1.61  1.56  1.51 
kim.CLR.ctr <- kim.CLR.var$LRvars / sum(kim.CLR.var$LRvars) 
kim.CLR.ctr.order <- order(kim.CLR.ctr, decreasing=TRUE)
round(100*kim.CLR.ctr[kim.CLR.ctr.order], 2)
#    Na     K    Rb    Ca    Ni    Co     P    La    Mg    Th    Nb    Fe    Al 
# 25.50 17.11 12.98  6.34  4.73  3.44  3.24  3.23  3.10  2.59  2.42  2.11  2.04 
#    Si    Ti    Yb    Cr    Zr     V    Er    Ga     Y 
#  2.02  1.73  1.36  1.24  1.02  1.00  0.99  0.96  0.85

### Percentages of logratio variance explained by each pairwise logratio
kim.LR.expvar <- rep(0, ncol(kim.LR$LR))
for(j in 1:ncol(kim.LR$LR)) {
  foo <- rda(kim.LR$LR ~ kim.LR$LR[,j]) 
  kim.LR.expvar[j] <- foo$CCA$tot.chi / foo$tot.chi
}
names(kim.LR.expvar) <- colnames(kim.LR$LR)
kim.LR.expvar.order <- order(kim.LR.expvar, decreasing=TRUE)
head(round(100*kim.LR.expvar[kim.LR.expvar.order],2), 10)
# Table 1, column (ii)
#  Mg/K  Fe/K  K/Co   K/V  K/Cr  Rb/V  Ti/K Mg/Rb Fe/Rb  Si/K 
# 47.28 47.11 46.77 46.53 46.00 45.97 45.66 45.55 45.55 45.52 

### To get percentages of between-group logratio variance, the easiest is
### to set up a response matrix where the individual logratios are 
### replaced by the mean logratios of the five phases

### First, compute the mean logratios of the phases
kim.LR.phase <- aggregate(kim.LR$LR ~ factor(kim270$StratUnit), FUN=mean)
rownames(kim.LR.phase) <- kim.LR.phase[,1]
kim.LR.phase <- kim.LR.phase[,-1]
kim.LR.phase[, 1:8]   # (to see some of the logratio means)

### Second, replace the individual logratios by their phase means
kim.LR.phase.all <- matrix(0, nrow=nrow(kim), ncol=ncol(kim.LR$LR))
colnames(kim.LR.phase.all) <- colnames(kim.LR.phase)
cumsum(table(factor(kim270$StratUnit)))
# Cantuar     eJF     lJF     mJF   Pense 
#      21     175     203     243     270 
for(i in 1:21)    kim.LR.phase.all[i,] <- as.numeric(kim.LR.phase[1,])
for(i in 22:175)  kim.LR.phase.all[i,] <- as.numeric(kim.LR.phase[2,])
for(i in 176:203) kim.LR.phase.all[i,] <- as.numeric(kim.LR.phase[3,])
for(i in 204:243) kim.LR.phase.all[i,] <- as.numeric(kim.LR.phase[4,])
for(i in 244:270) kim.LR.phase.all[i,] <- as.numeric(kim.LR.phase[5,])

### Third, compute the explained variance between groups by each logratio
kim.LR.expvar.phase <- rep(0, ncol(kim.LR.phase.all))
names(kim.LR.expvar.phase) <- colnames(kim.LR$LR)
for(j in 1:ncol(kim.LR.phase.all)) {
  foo <- rda(kim.LR.phase.all ~ kim.LR.phase.all[,j]) 
  kim.LR.expvar.phase[j] <- foo$CCA$tot.chi / foo$tot.chi
}
kim.LR.expvar.phase.order <- order(kim.LR.expvar.phase, decreasing=TRUE)
head(round(100*kim.LR.expvar.phase[kim.LR.expvar.phase.order],2), 10)
# Table 1, column (iii)
#  K/Co  Fe/K  K/Ni  V/Yb  Mg/K Rb/Co Cr/Ga Al/Mg Fe/Rb Rb/Ni 
# 56.44 56.16 56.00 55.91 55.67 55.62 55.56 55.44 55.37 55.03 

### Table 1 in summary
Table1 <- 
cbind(colnames(kim.LR$LR)[kim.LR.ctr.order],round(100*kim.LR.ctr[kim.LR.ctr.order],2),
      colnames(kim.LR$LR)[kim.LR.expvar.order],round(100*kim.LR.expvar[kim.LR.expvar.order], 2),
      colnames(kim.LR$LR)[kim.LR.expvar.phase.order],round(100*kim.LR.expvar.phase[kim.LR.expvar.phase.order], 2))[1:10,]
rownames(Table1) <- 1:10
Table1  # shows the values in Table 1


### ----------
### Star plots
star.elts <- kim[,c("Si","Mg","Fe","Cr","Co","Ni","Ti",      # Mantle
                    "Al","Rb","Na","K","Ga",                 # Crust
                    "Nb","La","Th","Zr","P","Er","Yb",       # Kimberlite
                    "Ca","Y","V")]                           # Unused
star.ends <- matrix(c(1,7, 8,12, 13,19, 20,22), nrow=4, ncol=2, byrow=TRUE)    
star.cols <- c(rep("blue",7), rep("red",5), rep("green", 7), rep("gray",3))
ng <-5
nv <- 22
temp<-matrix(nrow=ng,ncol=nv)
geolunits <- c("Cantuar", "Pense", "eJF", "mJF", "lJF")
temp<-as.data.frame(temp)
names(temp)<-colnames(star.elts)
for (i in 1:ng) {	
  x <- star.elts[kim270$StratUnit==geolunits[i],]
  temp[i,] <- colMeans(as.matrix(x))
}
### Figure 4 star plots
palette(star.cols)
stars(temp, labels=geolunits, draw.segments=TRUE, key.loc=c(5,2), axes=FALSE,
      xlim=c(1,6.5),ylim=c(0.5,8),
      frame.plot=FALSE,scale=TRUE,cex=0.9)
# (the arcs and labelling of the amalgamations was done in Powerpoint)

### --------------------------------------------------
### Ternary plots, logratios and logs of amalgamations

### Amalgamations
colnames(kim)
 [1] "Si" "Ti" "Al" "Fe" "Mg" "Ca" "Na" "K"  "P"  "Rb" "Nb" "Zr" "Th" "V"  "Cr" "Co" "Ni" "La" "Er" "Yb" "Y"  "Ga"
# Mantle Contamination -- Si+Mg+Fe+Cr+Co+Ni+Ti
# Crustal Contamination -- Al+Rb+Na+K+Ga
# Kimberlite Fractionation - Nb+La+Th+Zr+P+Er+Yb
mantle <- rowSums(kim[,c(1,5,4,15,16,17,2)]) 
crustal <- rowSums(kim[,c(3,10,7,8,22)])
kimberlite <- rowSums(kim[,c(11,18,13,12,9,19,20)])
rest <- 1-rowSums(cbind(mantle, crustal, kimberlite))
amalgs <- CLOSE(cbind(mantle, crustal, kimberlite))
amalg.means <- aggregate(as.matrix(amalgs) ~ kim270$StratUnit, FUN=mean)
rownames(amalg.means) <- amalg.means[,1]
amalg.means <- amalg.means[kim.su.order,-1]
round(100 * amalg.means, 2)     
#         mantle crustal kimberlite     # Table 2
# Cantuar  95.69    3.13       1.18
# Pense    97.00    2.45       0.55
# eJF      96.70    2.97       0.32
# mJF      95.21    4.40       0.39
# lJF      94.36    5.14       0.50

### Ternary plot of top corner
require(Ternary)
# Data object
dat <- data.frame(mantle = 100*amalgs[,1],
                  crustal = 100*amalgs[,2],
                  kimberlite = 100*amalgs[,3])
### Compute coordinates of top corner 
TernaryCoords(c(100,0,0)) # 0.0000000 0.8660254
TernaryCoords(c(80,20,0)) # 0.1000000 0.6928203
TernaryCoords(c(80,0,20)) #-0.1000000  0.6928203
# Figure 5A (axis labels and arrows added in Powerpoint)
par(mar = rep(0.3, 4))
TernaryPlot("Mantle", "Crustal", "Kimberlite", 
            xlim = c(-0.03, 0.03), ylim = c(0.8, 0.8666), padding = 0.01,
            grid.lines = 20, grid.lty = "dotted",
            grid.minor.lines = 5, grid.minor.lty = "dotted")
pointSize <- kim.su.cex[kim.su]
pointCol <- kim.su.col[kim.su]
TernaryPoints(dat[, c("mantle", "crustal", "kimberlite")],
              cex = kim.su.cex[kim.su],      # Point size
              col = kim.su.col[kim.su],      # Point colour
              bg  = kim.su.bg[kim.su],       # Point background colour
              pch = kim.su.pch[kim.su],      # Point character
              )

### Centered ternary plot
dat2 <- scale(log(dat), scale=FALSE)
dat2 <- CLOSE(exp(dat2))
# Figure 5B (axis labels and arrows added in Powerpoint)
par(mar = rep(0.3, 4))
# TernaryPlot(alab = "Mantle \u2192", blab = "Crustal \u2192", clab = "\u2190 Kimberlite",
TernaryPlot("","","",
            grid.lines = 0, #grid.lty = "dotted",
            grid.minor.lines = 0)  #, grid.minor.lty = "dotted")
TernaryPoints(100*dat2[, c("mantle", "crustal", "kimberlite")],
              cex = kim.su.cex[kim.su],   
              col = kim.su.col[kim.su],   
              pch = kim.su.pch[kim.su],  
              bg  = kim.su.bg[kim.su]
              )

### Plotting two logratios of amalgamations, with convex hulls and confidence ellipses
crustal.mantle <- log(crustal/mantle)
kimberlite.mantle <- log(kimberlite/mantle)
# Figure 5C (some labels enhanced afterwards in graphics editor)
par(mar=c(4,4,2,1), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, las=1)
plot(crustal.mantle, kimberlite.mantle, type="n",xlab="Crustal/Mantle (log)", ylab="Kimberlite/Mantle (log)")
points(crustal.mantle, kimberlite.mantle, col=kim.su.col[kim.su], pch=kim.su.pch[kim.su], bg="white",# bg=kim.su.bg[kim.su],
       cex=kim.su.cex[kim.su])
for (j in 1:5){
 points     <- cbind(crustal.mantle, kimberlite.mantle)[kim.su==j,]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,], lty = 3, col=kim.su.col[j])
}
require(ellipse)
set.seed(123)
CIplot_biv(crustal.mantle, kimberlite.mantle, group=kim.su, groupnames=c("Cantaur","eJF","lJF","mJF","Pense"),
           groupcols=kim.su.col, shade=TRUE, add=TRUE, cex=0.8, alpha=0.99)
set.seed(123)
CIplot_biv(crustal.mantle, kimberlite.mantle, group=kim.su, groupnames=c("Cantaur","eJF","lJF","mJF","Pense"),
           groupcols=kim.su.col, add=TRUE, cex=0.8, alpha=0.99)
legend("topleft", legend=c("Cantuar","eJF","lJF","mJF","Pense")[kim.su.order], 
       pch=kim.su.pch[kim.su.order], col=kim.su.col[kim.su.order], 
       text.col=kim.su.col[kim.su.order], cex=0.8, bty="n")

### Plotting two log-transformed amalgamations, with convex hulls and confidence ellipses
log.crustal <- log(crustal)
log.kimberlite <- log(kimberlite)
# Figure 5D (some labels enhanced afterwards in graphics editor)
par(mar=c(4,4,2,1), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, las=1)
plot(log.crustal, log.kimberlite, type="n",xlab="Crustal (log)", ylab="Kimberlite (log)")
points(log.crustal, log.kimberlite, col=kim.su.col[kim.su], pch=kim.su.pch[kim.su], bg="white",# bg=kim.su.bg[kim.su],
       cex=kim.su.cex[kim.su])
for (j in 1:5){
 points     <- cbind(log.crustal, log.kimberlite)[kim.su==j,]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,], lty = 3, col=kim.su.col[j])
}
require(ellipse)
set.seed(123)
CIplot_biv(log.crustal, log.kimberlite, group=kim.su, groupnames=c("Cantaur","eJF","lJF","mJF","Pense"),
           groupcols=kim.su.col, shade=TRUE, add=TRUE, cex=0.8, alpha=0.99)
set.seed(123)
CIplot_biv(log.crustal, log.kimberlite, group=kim.su, groupnames=c("Cantaur","eJF","lJF","mJF","Pense"),
           groupcols=kim.su.col, add=TRUE, cex=0.8, alpha=0.99)
legend("topleft", legend=c("Cantuar","eJF","lJF","mJF","Pense")[kim.su.order], 
       pch=kim.su.pch[kim.su.order], col=kim.su.col[kim.su.order], 
       text.col=kim.su.col[kim.su.order], cex=0.8, bty="n")


### ---------------------------------------------------------------------------
### Hierarchical cluster analysis of all the samples,  with horizontal plotting
### Ward clustering the CLRs (CLRs already computed above, repeated here)
### The horizontal plotting needs package 'ape'
kim.CLR <- CLR(kim, weight=FALSE)
kim.CLR.clust <- hclust(dist(kim.CLR$LR), method="ward.D2")
# recall point sizes and order of phases
kim.su.cex <- c(0.8, 1, 0.7, 0.6, 0.7)  
kim.su.order <- c(1,5,2,4,3)
require(ape)
# only use the following statement if you want to save the PDF of the clustering
# pdf("AHC_Ward.pdf", width=7, height=12)
### Figure 6
par(xpd=TRUE, mar=c(0.5,1,0.5,0.1))
plot(as.phylo(kim.CLR.clust), font=2, tip.color = kim.su.col[kim.su], label.offset = 0.1,
     cex=0.3, pch=20, bg=kim.su.bg[kim.su], show.tip.label=FALSE)   
tiplabels(text="", pch=20, cex=0.5, col=kim.su.col[kim.su], frame="none", bg=kim.su.bg[kim.su], adj=c(0.6, 0.5))
segments(7.6,0, 7.6, 270, lty=2, lwd=2, col="gray")
legend("topleft", legend=c("Cantaur","Pense","eJF","mJF","lJF"), text.col=kim.su.col[kim.su.order], 
       pch=20, col=kim.su.col[kim.su.order], pt.cex=0.9, text.font=2, bty="n")  
# only use the following to close the PDF file if saved
# dev.off()
### Cut tree at 5 clusters
kim.clust5 <- cutree(kim.CLR.clust, 5)
### Table 3
table(kim.clust5, kim270$StratUnit)[,kim.su.order]
# kim.clust5 Cantuar Pense eJF mJF lJF             
#          1       0     0  32  18   1
#          2      13     9  10   0   0
#          3       0    12  37   0   0
#          4       0     1  71   0   0
#          5       8     5   4  22  27
### Sum of maximum cluster assignments
(32+13+37+71+27)/270   # Pense and mJF not assigned
# [1] 0.6666667

### ------------------------------------------------------------
### Non-hierarchical k-means cluster analysis of all the samples
### Looping on k-means algorithm to decide how many clusters
set.seed(1234)
kim.BW <- rep(0, 10)
for(nc in 2:10) {
  kim.kmeans <- kmeans(kim.CLR$LR, centers=nc, nstart=20, iter.max=200)
  kim.BW[nc] <- kim.kmeans$betweenss/kim.kmeans$totss
}
kim.BW
# [1] 0.0000000 0.3987519 0.5087816 0.5756866 0.6227841 0.6503426 0.6742384 0.6988051 0.7189742 0.7324721
### Plot the proportion of between-cluster variance and the increments in between-cluster variance
### Plot with bars
### Figure 7
par(mar=c(4.2,4,1,2), mgp=c(2.5,0.7,0), las=1, font.lab=2, cex.lab=0.9, cex.axis=0.8, mfrow=c(1,2))
kim.BW[1] <- 0.005
plot(kim.BW, xlab="Number of clusters", ylab="BSS/TSS", type="n", bty="n", xaxt="n", xlim=c(1.5,10.5), ylim=c(0,0.7))
axis(1, at=2:10, labels=2:10, tick=FALSE)
for(g in 2:10) segments(g, 0, g, kim.BW[g], lwd=10, col="gray", lend=1)
kim.BWinc <- kim.BW[2:10]-kim.BW[1:9]  
plot(2:10, kim.BWinc[1:9], type="n", xlab="Number of clusters",  bty="n", xaxt="n", xlim=c(2.5,12), 
     ylim=c(0,0.12), ylab="Improvement in BSS/TSS")
axis(1, at=3:10, labels=3:10, tick=FALSE)
for(g in 3:5) segments(g, 0, g, kim.BWinc[g-1], lwd=10, col=c(rep("pink",4),rep("lightblue",4)), lend=1)
for(g in 6:10) segments(g, 0, g, kim.BWinc[g-1], lwd=10, col=c(rep("lightblue",4),rep("lightblue",4)), lend=1)
### (looks like 5-cluster solution is a good choice)
set.seed(1234)
kim.kmeans5 <- kmeans(kim.CLR$LR, centers=5, nstart=20, iter.max=200)
kim.kmeans5$betweenss/kim.kmeans5$totss
# [1] 0.6227841
### cluster sizes
kim.kmeans5$size
# [1] 71 32 70 44 53
### Table 4
table(kim.kmeans5$cluster, kim270$StratUnit)[,kim.su.order]
#     Cantuar Pense eJF mJF lJF
#   1       0     0  35  32   4
#   2      13    11   8   0   0
#   3       0     0  70   0   0
#   4       8     4   0   8  24
#   5       0    12  41   0   0
### Sum of maximum cluster assignments
(35+13+70+24+41)/270
# [1] 0.6777778

### ------------------------------------
### Amalgamation clustering of the parts
kim.aclust <- ACLUST(kim, weight=FALSE)
### Figure 8 (label of top level of 0.1132 added in Powerpoint)
plot(kim.aclust, hang=-1, ylab="Height = Logratio variance lost", main="Amalgamation clustering")

### ------------------------------------------
### Logratio analysis (LRA) of kimberlite data
### (LRA = PCA of the CLRs, equivalently PCA of all LRs)
kim.LRA <- LRA(kim, weight=FALSE)
### Row principal coordinates
kim.rpc <- kim.LRA$rowpcoord
### Column coordinates 
kim.crd <- kim.LRA$colcoord 
percs   <- 100 * kim.LRA$sv^2 / sum(kim.LRA$sv^2)
### Figure 9 (some labels enhanced in graphics editor)
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
rescale <- 0.25
plot(1.05 * rbind(kim.rpc, rescale * kim.crd), type = "n", asp = 1,
     xlab = paste("LRA dimension ", 1, " (", round(percs[1], 1), "%)", sep = ""), 
     ylab = paste("LRA dimension ", 2, " (", round(percs[2], 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = "red", col.axis = "red")
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
     col = "black", col.ticks = "red", col.axis = "red")
points(kim.rpc, pch = kim.su.pch[kim.su], col = kim.su.col[kim.su], font = 1, cex = kim.su.cex[kim.su])
text(rescale * kim.crd, labels = colnames(kim), col = "red", cex = 0.9, font = 4)    
legend("bottomright", legend=c("Cantaur","Pense","eJF","mJF","lJF"), bty="n",
       pch=kim.su.pch[kim.su.order], col=kim.su.col[kim.su.order], pt.cex=kim.su.cex[kim.su.order], 
       text.col=kim.su.col[kim.su.order], cex=0.9, text.font=2)
for (j in 1:5){
 points     <- kim.rpc[kim.su==j,1:2]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,], lty = 3, col=kim.su.col[j])
}
require(ellipse)
set.seed(123)
CIplot_biv(kim.rpc[,1], kim.rpc[,2], group=kim.su, groupcols=kim.su.col, 
           add=TRUE, shade=TRUE, alpha=0.995, 
           groupnames=c("Cantuar","eJF","lJF","mJF","Pense"))
CIplot_biv(kim.rpc[,1], kim.rpc[,2], group=kim.su, groupcols=kim.su.col, 
           alpha=0.995, add=TRUE, shade=FALSE, shownames=FALSE)

### ----------------------------------------------------
### Principal component analysis (PCA) of ALRs w.r.t. Zr
### First identify Zr as being the best reference element
### (most isometric) for an additive logratio transform
FINDALR(kim, weight=FALSE)
# $totvar
# [1] 0.1132472
# $procrust.cor
#  [1] 0.9657580 0.9586441 0.9373568 0.9678512 0.9551558 0.8232680 0.8346451 0.8769878 0.9194245 0.8892915
# [11] 0.9623413 0.9772571 0.9435395 0.9721373 0.9732015 0.9451002 0.9182467 0.9412665 0.9609117 0.9383257
# [21] 0.9769708 0.95452000
# $procrust.max
# [1] 0.9772571
$procrust.ref
[1] 12
# $var.log
#  [1] 0.003046547 0.061977127 0.079009722 0.004858624 0.006784848 0.227910235 0.852322165 0.606497843
#  [9] 0.149071830 0.494010137 0.119399769 0.095607494 0.144020770 0.042618539 0.016720746 0.017030990
# [17] 0.028567142 0.160182671 0.095157478 0.102457123 0.096669090 0.061117178
# $var.min
# [1] 0.003046547
# var.ref
# [1] 1
# (the best reference is part number 12) 
kim.ALR <- ALR(kim, denom=12, weight=FALSE)
kim.PCA <- PCA(kim.ALR$LR, weight=FALSE)
### Row principal coordinates
kim.alr.rpc <- kim.PCA$rowpcoord
### Column standard coordinates 
kim.alr.csc <- kim.PCA$colcoord 
percs.PCA  <- 100 * kim.PCA$sv^2 / sum(kim.PCA$sv^2)
### Figure 10 (some labels enhanced in graphics editor)
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
rescale <- 0.25
# Note that both dimensions are reversed to facilitate comparison with the LRA
plot(1.05 * rbind(-kim.alr.rpc, rescale * (-kim.alr.csc)), type = "n", asp = 1,
     xlab = paste("PCA dimension ", 1, " (", round(percs.PCA[1], 1), "%)", sep = ""), 
     ylab = paste("PCA dimension ", 2, " (", round(percs.PCA[2], 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = "red", col.axis = "red")
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
     col = "black", col.ticks = "red", col.axis = "red")
points(-kim.alr.rpc, pch = kim.su.pch[kim.su], col = kim.su.col[kim.su], font = 1, cex = kim.su.cex[kim.su])
text(rescale * (-kim.alr.csc), labels = colnames(kim.ALR$LR), col = "red", cex = 0.9, font = 4)    
legend("bottomright", legend=c("Cantaur","Pense","eJF","mJF","lJF"), bty="n",
       pch=kim.su.pch[kim.su.order], col=kim.su.col[kim.su.order], pt.cex=kim.su.cex[kim.su.order], 
       text.col=kim.su.col[kim.su.order], cex=0.9, text.font=2)
for (j in 1:5){
 points     <- -kim.alr.rpc[kim.su==j,1:2]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,], lty = 3, col=kim.su.col[j])
}
require(ellipse)
set.seed(123)
CIplot_biv(-kim.alr.rpc[,1], -kim.alr.rpc[,2], group=kim.su, groupcols=kim.su.col, 
           add=TRUE, shade=TRUE, alpha=0.995, 
           groupnames=c("Cantuar","eJF","lJF","mJF","Pense"))
CIplot_biv(-kim.alr.rpc[,1], -kim.alr.rpc[,2], group=kim.su, groupcols=kim.su.col, 
           alpha=0.995, add=TRUE, shade=FALSE, shownames=FALSE)	

### ------------------------------------------------------------
### Supervised analysis predicting eJF using logistic regression
### with default Bonferroni stopping criterion
eJF <- rep(0, nrow(kim))
eJF[kim270$StratUnit == "eJF"] <- 1
table(eJF)
#   0   1 
# 116 154 
set.seed(1234567)
require(easyCODA)
### method 1, any logratio can enter
eJF.STEPR <- STEPR(kim, as.factor(eJF), method=1, family="binomial")
eJF.STEPR$names
[1] "Mg/La" "Mg/V" 
# (note: if you get "Mg/La" "V/La", this is equivalent to the above,
#  there is a tie between "Mg/V" and "V/La" for the second ratio and
#  one is chosen randomly) 
eJF.STEPR$Bonferroni
# [1] 117.3022 100.0730 102.2857
eJF.MgLa.MgV <- glm(factor(eJF.Rest) ~ kim.LR$LR[,"Mg/La"] + kim.LR$LR[,"Mg/V"], family="binomial")
summary(eJF.MgLa.MgV)
# Coefficients:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -228.807     40.901  -5.594 2.22e-08 ***
# kim.LR$LR[, "Mg/La"]   14.016      2.579   5.436 5.46e-08 ***
# kim.LR$LR[, "Mg/V"]    13.146      3.271   4.019 5.85e-05 ***
eJF.MgLa.MgV.pred <- predict(eJF.MgLa.MgV)
table(eJF.MgLa.MgV.pred>0, factor(eJF.Rest))
#           0   1
#   FALSE 107   9
#   TRUE    9 145
(107+145)/270
# [1] 0.9333333

eJF.STEPR2 <- STEPR(kim, as.factor(eJF), method=2, family="binomial")
eJF.STEPR2$names
# [1] "Mg/La" "V/Ni" 
eJF.STEPR2$Bonferroni
# [1] 117.3022 100.5510 102.0931
eJF.MgLa.NiV <- glm(factor(eJF) ~ kim.LR$LR[,"Mg/La"] + kim.LR$LR[,"V/Ni"], family="binomial")
summary(eJF.MgLa.NiV)
# Coefficients:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -119.103     19.171  -6.213 5.21e-10 ***
# kim.LR$LR[, "Mg/La"]   10.202      2.083   4.899 9.65e-07 ***
# kim.LR$LR[, "V/Ni"]   -11.033      2.734  -4.036 5.45e-05 ***
# (we will use Ni/V rather than V/Ni for the plot to ensure both coefficients are positive)
eJF.MgLa.NiV.pred <- predict(eJF.MgLa.NiV)
table(eJF.MgLa.NiV.pred>0, factor(eJF))
        0   1
  FALSE 108   9
  TRUE    8 145
(108+145)/270
# [1] 0.937037

### Two plots of two predicting logratios
### Figure 11
par(mar=c(4.2,4,1,1), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, mfrow=c(1,2))
plot(kim.LR$LR[,"Mg/La"],kim.LR$LR[,"Mg/V"], type="n", asp=1, xlab="log(Mg/La)", ylab="log(Mg/V)")
points(kim.LR$LR[,"Mg/La"],kim.LR$LR[,"Mg/V"], pch=21, col=c("hotpink1","cyan2")[eJF+1], 
       bg=c("hotpink1","cyan2")[eJF+1], cex=0.7)
for (j in 1:2){
 points <- cbind(kim.LR$LR[,"Mg/La"],kim.LR$LR[,"Mg/V"])[eJF+1==j,]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,],lty = 3, col=c("hotpink1","cyan3")[j])
}
legend("bottomright", legend=c("eJF phase","other phases"), bty="n",
       pch=21, col=c("cyan2","hotpink1"), pt.bg=c("cyan2","hotpink1"), text.col=c("cyan3","hotpink1"), 
       text.font=2, cex=0.9, pt.cex=0.7)

plot(kim.LR$LR[,"Mg/La"],-kim.LR$LR[,"V/Ni"], type="n", asp=1, xlab="log(Mg/La)", ylab="log(Ni/V)")
points(kim.LR$LR[,"Mg/La"],-kim.LR$LR[,"V/Ni"], pch=21, col=c("hotpink1","cyan2")[eJF+1], 
       bg=c("hotpink1","cyan2")[eJF+1], cex=0.7)
for (j in 1:2){
 points <- cbind(kim.LR$LR[,"Mg/La"],-kim.LR$LR[,"V/Ni"])[eJF+1==j,]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,],lty = 3, col=c("hotpink1","cyan3")[j])
}
legend("bottomright", legend=c("eJF phase","other phases"), bty="n",
       pch=21, col=c("cyan2","hotpink1"), pt.bg=c("cyan2","hotpink1"), text.col=c("cyan3","hotpink1"), 
       text.font=2, , cex=0.9, pt.cex=0.7)

### -----------------------------------------------------------------------------------------
### Plot of second pair of logratios, with contours of prediction function (all untransformed)
### (notice that the inverse of V/Ni is used to obtain Ni/V that has a positive coefficient --
### see the logistic regression resujlts above)
### Figure 12
par(mar=c(4.2,4,1,1), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, mfrow=c(1,1))
plot(exp(kim.LR$LR[,"Mg/La"]),exp(-kim.LR$LR[,"V/Ni"]), type="n", xlab="Mg/La", ylab="Ni/V")
points(exp(kim.LR$LR[,"Mg/La"]),exp(-kim.LR$LR[,"V/Ni"]), pch=21, col=c("hotpink1","cyan2")[eJF+1], 
       bg=c("hotpink1","cyan2")[eJF+1], cex=0.7)
for (j in 1:2){
 points     <- cbind(exp(kim.LR$LR[,"Mg/La"]),exp(-kim.LR$LR[,"V/Ni"]))[eJF+1==j,]
 hpts       <- chull(points)
 hpts       <- c(hpts, hpts[1])
 lines(points[hpts,],lty = 3, col=c("hotpink1","cyan3")[j])
}
X <- seq(1.05*min(exp(kim.LR$LR[,"Mg/La"])), 0.95*max(exp(kim.LR$LR[,"Mg/La"])), length=50)
Y <- seq(1.05*min(exp(-kim.LR$LR[,"V/Ni"])), 0.95*max(exp(-kim.LR$LR[,"V/Ni"])), length=50)
XY <- cbind(rep(X, each=50), rep(Y,50))
XY.pred <- exp(-119.103+10.202*log(XY[,1])+11.033*log(XY[,2])) / (1+exp(-119.103+10.202*log(XY[,1])+11.033*log(XY[,2]))) 
ordisurf(XY, XY.pred, col="grey20", add=TRUE, levels=c(0.01,seq(0.1,0.9,0.1),0.95, 0.99), labcex=0.7, lty=2)
text(2500, 12.5, "Rest", font=2, cex=1.2, col="hotpink1")
text(7700, 28.5, "eJF", font=2, cex=1.2, col="cyan2")

### ---------------------------------------------------------------
### Classification tree to predict eJF phase, using package 'rpart‘
require(rpart)
### Figure 13 classification tree, using pairwise logratios
### The labelling has been enhanced afterwards
eJF.tree <- rpart(factor(eJF) ~ ., data=as.data.frame(kim.LR$LR))
plot(eJF.tree, mar=0.1)
text(eJF.tree, use.n=TRUE)
eJF.tree <- rpart(factor(eJF) ~ ., data=as.data.frame(kim.ratios))
plot(eJF.tree, mar=0.1)
text(eJF.tree, use.n=TRUE)
### Predictions of eJF 
eJF.pred   <- predict(eJF.tree)
eJF.pred01 <- eJF.pred[,2]>0.5
table(eJF.pred01, eJF)
#             eJF
# eJF.pred   0   1
#    FALSE 111   2
#    TRUE    5 152
(111+152)/270
# [1] 0.9740741
### 10-fold cross-validation; notice that luckily 270 divides exactly by 10
set.seed(123)
groups10 <- rep(1:10, each=27)
groups10 <- groups10[sample(1:270)]
cvpreds <- rep(0, nrow(kim))
for(cv in 1:10) {
  cvset <- (1:nrow(kim))[-which(groups10==cv)]
  foo <- rpart(factor(eJF[cvset]) ~ ., data=as.data.frame(kim.LR$LR[cvset,]) )
  cvset <- which(groups10==cv)
  foo.pred <- predict(foo, newdata=as.data.frame(kim.LR$LR[cvset,]))
  cvpreds[cvset] <- foo.pred[,2]
}
table(cvpreds>0.5, factor(eJF))
#           0   1
#   FALSE 106   9
#   TRUE   10 145
(106+145)/270
# [1] 0.9296296

### ----------------------------------------------------------
### Multinomial logistic regression predicting all five phases
### Uses package 'nnet'. Stepwise procedure followed step-by-step
require(nnet)
### Uses Bonferroni as stopping criterion with nmumber of 
### parameters m = 4 x the estimated parameters in one equation
### Bonferroni = (AIC - 2m = -2 log Likelihood) + 9.23 * m
multlogit <- rep(0, 231)
names(multlogit) <- colnames(kim.LR$LR) 
for(jk in 1:231){
  test <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Si/Zr 
#    11
step1 <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,"Si/Zr"])
summary(step1)
# AIC: 332.1833 
332.1833 - 2*8 + 9.23*8
# [1] 390.0233
sum(diag(table(predict(step1), kim270$StratUnit)))/270
# [1] 0.7740741

# step2
multlogit <- rep(999, 231)
names(multlogit) <- colnames(kim.LR$LR) 
for(jk in (1:231)[-11]){
  test <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,"Si/Zr"] + kim.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# P/Rb 
#  141 
step2 <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,"Si/Zr"]+ kim.LR$LR[,"P/Rb"])
summary(step2)
# AIC: 174.704 
174.704 - 2*12 + 9.23*12
# [1] 261.464
sum(diag(table(predict(step2), kim270$StratUnit)))/270
table(predict(step2), kim270$StratUnit)[kim.su.order,kim.su.order]
# [1] 0.8962963

# step3
multlogit <- rep(999, 231)
names(multlogit) <- colnames(kim.LR$LR) 
for(jk in (1:231)[-c(11,141)]){
  test <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,"Si/Zr"] + kim.LR$LR[,"P/Rb"] + kim.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Si/Ni
#  16 
step3 <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,"Si/Zr"]+kim.LR$LR[,"P/Rb"]+kim.LR$LR[,"Si/Ni"])
summary(step3)
# AIC: 148.5424
148-54.704 - 2*16 + 9.23*16
# [1] 208.976
sum(diag(table(predict(step3), kim270$StratUnit)))/270
# [1] 0.9296296

# step4
multlogit <- rep(999, 231)
names(multlogit) <- colnames(kim.LR$LR) 
for(jk in (1:231)[-c(11,16,141)]){
  test <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,"Si/Zr"] + kim.LR$LR[,"P/Rb"] + kim.LR$LR[,"Si/Ni"] + kim.LR$LR[,jk])
  multlogit[jk] <- AIC(test)
}
which(multlogit == min(multlogit))
# Mg/Cr
#  88 
step4 <- multinom(factor(kim270$StratUnit) ~ kim.LR$LR[,"Si/Zr"]+kim.LR$LR[,"P/Rb"]+kim.LR$LR[,"Si/Ni"]+kim.LR$LR[,"Mg/Cr"])
summary(step4)
Coefficients:
      (Intercept) kim.LR$LR[, "Si/Zr"] kim.LR$LR[, "P/Rb"] kim.LR$LR[, "Si/Ni"] kim.LR$LR[, "Mg/Cr"]
eJF    -348.20474            40.166364          -0.4104366            -10.51273            13.006976
lJF     -85.16834             8.851836          -3.5485629             16.26679           -10.740482
mJF    -204.86144            17.358544          -4.4333070             10.28428             5.229471
Pense  -146.04069            25.380807           5.3624595            -22.98457             5.091146

Std. Errors:
      (Intercept) kim.LR$LR[, "Si/Zr"] kim.LR$LR[, "P/Rb"] kim.LR$LR[, "Si/Ni"] kim.LR$LR[, "Mg/Cr"]
eJF      71.53365             8.139722            1.847785             8.250520             6.652180
lJF      52.49849             5.960789            1.793604             8.280532             7.673767
mJF      60.82092             6.494835            1.762123             7.712160             7.679369
Pense    61.84604             7.927122            1.813085             7.938223             5.923046

Residual Deviance: 90.83202   # 90.832++8.166*16 = 221.488
AIC: 130.832 

130.83 - 2*20 + 9.23*20
# [1] 275.43

### hence, stop after three terms
summary(step3)
Coefficients:
      (Intercept) kim.LR$LR[, "Si/Zr"] kim.LR$LR[, "P/Rb"] kim.LR$LR[, "Si/Ni"]
eJF    -270.95534            38.865838          -0.2607031            -9.539437
lJF     -55.29997             3.733136          -4.2146014             8.738966
mJF    -132.33121            14.988885          -4.3254225             5.991670
Pense  -108.36150            23.223598           4.9639302           -21.049631

Std. Errors:
      (Intercept) kim.LR$LR[, "Si/Zr"] kim.LR$LR[, "P/Rb"] kim.LR$LR[, "Si/Ni"]
eJF      60.46482             6.969529            1.681542             6.889565
lJF      46.60317             3.742451            1.671015             5.965183
mJF      50.38389             4.401147            1.617929             5.837064
Pense    51.68934             6.562230            1.706403             7.256361

Residual Deviance: 116.5424 
AIC: 148.5424 
	

### -----------------------------------------------------------------
### Classification tree on pairwise ratios (not logratios) predicting
### all five phases
require(rpart)
kim.ratios <- as.data.frame(exp(LR(kim)$LR))
kim.tree <- rpart(factor(kim270$StratUnit) ~ ., data=kim.ratios)
### Figure 14 classification tree using pairwise ratios
### The labelling has been enhanced afterwards
plot(kim.tree, mar=0.1)
text(kim.tree, use.n=TRUE)
### Predictions of the phases
kim.phase.pred <- rep(0, nrow(kim))
kim.pred <- predict(kim.tree)
kim.pred.max <- apply(predict(kim.tree), 1, max)
for(i in 1:nrow(kim)) kim.phase.pred[i] <- which(kim.pred[i,]==kim.pred.max[i])
table(kim.phase.pred, kim270$StratUnit)[kim.su.order,kim.su.order]
# Table 6
# kim.phase.pred Cantuar Pense eJF mJF lJF
#              1      20     0   0   0   0
#              5       1    20   1   0   0
#              2       0     7 153   2   1
#              4       0     0   0  36   0
#              3       0     0   0   2  27
sum(diag(table(kim.phase.pred, kim270$StratUnit)))/nrow(kim)
# [1] 0.9481481

### cross-validation
set.seed(123)
groups10 <- rep(1:10, each=27)
groups10 <- groups10[sample(1:270)]
cvpreds <- matrix(0, nrow(kim), 5)
for(cv in 1:10) {
  cvset <- (1:270)[-which(groups10==cv)]
  foo <- rpart(factor(kim270$StratUnit[cvset]) ~ ., data=kim.ratios[cvset,])
  cvset <- which(groups10==cv)
  foo.pred <- predict(foo, newdata=kim.ratios[cvset,])
  cvpreds[cvset,] <- foo.pred
}
cv.pred.phase <-rep(0, nrow(kim))
for(i in 1:nrow(kim)) cv.pred.phase[i] <- which(cvpreds[i,] == max(cvpreds[i,]))
table(cv.pred.phase, kim270$StratUnit)[kim.su.order,kim.su.order]
# cv.pred.phase Cantuar Pense eJF mJF lJF
#             1      19     2   0   1   0
#             5       2    15   0   1   1
#             2       0     6 146   3   2
#             4       0     2   8  31   3
#             3       0     2   0   4  22
sum(diag(table(cv.pred.phase, kim270$StratUnit)))/nrow(kim)
# 37 mispredictions, 233 correct ones 
# [1] 0.862963



