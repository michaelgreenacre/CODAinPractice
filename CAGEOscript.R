### Install easyCODA package from CRAN in usual way or from R-Forge:

install.packages("easyCODA", repos="http://R-Forge.R-project.org")

### (version 0.29 of easyCODA was used here)

### 3-part data set of wine, beer and spirits assumed in data frame 'wbs'

### Model and plot of amalgamation logratio vs. logratio, and ILR vs. logratio
##  transformations
#       alog = amalgamation logratio
#       plog = single logratio
#       ilog = isometric logratio

attach(wbs)
alog <- log(spirits/(beer+wine))
plog <- log(beer/wine) 
# ILR "by hand" classical definition
ilog <- sqrt(2/3) * log(spirits / sqrt(beer*wine))  # using counts
# using function ILR in easyCODA weights are used
# so the result for option weight=FALSE (equal weighting) is 
# the classical definition divided by sqrt(number of parts) = sqrt(3)
ILR(wbs, numer=3, denom=c(1,2), weight=FALSE)$LR
# compare above result with
ilog/sqrt(3)

##  Model with amalgamation logratio and plotting of Fig. 1a

mod1 <- lm(alog ~ plog)
summary(mod1)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.85503    0.03796 -48.872  < 2e-16 ***
# plog         0.18915    0.04413   4.286 8.72e-05 ***

mod1.pred <- predict(mod1, type="response", se.fit=T, 
                     newdata=data.frame(plog=seq(range(plog)[1], range(plog)[2],length=100)))

par(mar=c(4.2,4,3,1), font.lab=2, cex.lab=1.2, mgp=c(2.7,0.7,0), las=1, mfrow=c(1,2))
plot(plog, alog, type="n", bty="n", xaxt="n", yaxt="n", xlab="log(beer/wine)",ylab="log[spirits/(beer+wine)]", 
     xlim=log(c(0.25,1)), ylim=c(-2.2,-1.8), las=1, main="Amalgamation logratio response")
axis(1)
axis(2)
segments(range(plog)[1], mod1.pred$fit[1], range(plog)[2], mod1.pred$fit[100], lwd = 3, col = "blue")
lines(seq(range(plog)[1], range(plog)[2],length=100), mod1.pred$fit+1.96*mod1.pred$se, lty = 3, lwd = 2,
      col = "blue")
lines(seq(range(plog)[1], range(plog)[2],length=100), mod1.pred$fit-1.96*mod1.pred$se, lty = 3, lwd = 2, 
      col = "blue")
symbols(plog, alog, fg="black", bg="white", circles=rep(0.02, length(alog)), inches=F, add=T, lwd=2)

##  Model with isometric logratio and plotting of Fig. 1b
#   Notice that plog is now divided by sqrt(2) to be the ILR definition

ilog2 <- plog/sqrt(2)
mod2 <- lm(ilog ~ ilog2)
summary(mod2)
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.99574    0.03066 -32.473   <2e-16 ***
# ilog2        0.01366    0.05042   0.271    0.788    

mod2.pred <- predict(mod2, type="response", se.fit=T, 
                     newdata=data.frame(ilog2=seq(range(ilog2)[1], range(ilog2)[2],length=100)))
plot(ilog2, ilog, type="n", bty="n", xaxt="n", yaxt="n", xlab="ILR(beer:wine)", ylab="ILR(spirits:beer,wine)", 
     xlim=c(-1.0, 0), ylim=c(-1.15,-0.85), las=1, main="ILR response")
axis(1) 
axis(2)
segments(range(ilog2)[1], mod2.pred$fit[1], range(ilog2)[2], mod2.pred$fit[100], lwd = 3, col = "blue")
lines(seq(range(ilog2)[1], range(ilog2)[2],length=100), mod2.pred$fit+1.96*mod2.pred$se, lty = 3, lwd = 2, 
      col = "blue")
lines(seq(range(ilog2)[1], range(ilog2)[2],length=100), mod2.pred$fit-1.96*mod2.pred$se, lty = 3, lwd = 2, 
      col = "blue")
symbols(ilog2, ilog, fg="black", bg="white", circles=rep(0.014, length(ilog)), inches=F, add=T, lwd=2)


### Ternary plot (Figure 2)
require(Ternary)

#   full plot (Fig.2a)
par(mfrow=c(1, 1), mar=rep(0.3, 4))
TernaryPlot(alab="Spirits \u2192", blab="Beer \u2192", clab="\u2190 Wine",
            point='up', lab.cex=1.5, grid.minor.lines = 0,
            grid.lty='solid', col=rgb(0.9, 0.9, 0.9), grid.col='white', 
            axis.col=rgb(0.6, 0.6, 0.6), ticks.col=rgb(0.6, 0.6, 0.6),
            padding=0.08)
AddToTernary(points, wbs[,c(3,2,1)], pch=21, cex=0.9, col="blue", bg="lightblue")

#   partial plot (Fig. 2b) with regression model in Fig. 1a back-transformed
for(i in 1:100) {
  beerwine.seq[i] <- beer.seq[i]/wine.seq[i]
  spirits.seq[i]  <- -1.855 + 0.1892 * log(beerwine.seq[i])
  spirits.seq[i]  <- exp(spirits.seq[i])/(1+exp(spirits.seq[i]))
}
wbs.add <- cbind(wine.seq, beer.seq, spirits.seq)
TernaryPlot(xlim=c(-0.3,-0.08), ylim=c(0,0.17), alab="Spirits \u2192", blab="Beer \u2192", clab="\u2190 Wine",
            point='up', lab.cex=1.5, grid.minor.lines = 0, grid.lty='solid', col=rgb(0.9, 0.9, 0.9),
            grid.col='white', axis.col=rgb(0.6, 0.6, 0.6), ticks.col=rgb(0.6, 0.6, 0.6), padding=0.22)
AddToTernary(points, wbs[,c(3,2,1)], pch=21, cex=0.9, col="blue", bg="lightblue")
AddToTernary(lines, wbs.add[,c(3,2,1)], lwd=2, col="gray30")


### For a fixed value of beer+wine=0.88, how does the geometric mean vary? (Figure 3)

beer.sim <- seq(range(beer)[1], range(beer)[2], length=100)
wine.sim <- 0.88-beer.sim
beerwine.gm.sim <- sqrt(beer.sim*wine.sim)

par(mar=c(4.2,4,3,1), font.lab=2, cex.lab=1.2, mgp=c(2.7,0.7,0), mfrow=c(1,1), las=1)
plot(beer.sim/wine.sim, beerwine.gm.sim, type="n", xlab="beer/wine ratio", ylim=c(0.36,0.44),
     ylab="geometric mean of beer and wine", main="For fixed sum beer+wine=0.88")
lines(beer.sim/wine.sim, beerwine.gm.sim, lwd=2, col="blue")

### Aar Massif data assumed in data.frame 'aar'

### Comparison of ILRs and SLRs (Figure 4) for aar and aar.sub (without MnO)

aar.sub <- aar[,-4]
aar.sub <- aar.sub / apply(aar.sub, 1, sum)

# SiO2, Na2O, MnO, P2O5 are numbers 1, 7, 4, 9 in aar
# SiO2, Na2O, P2O5 are numbers 1, 6, 8 in aar.sub

# remember that ILR in easyCODA divides the classic definition by the square root
# of the number of parts in the (sub)composition

ilr1 <- ILR(aar, numer=c(1,7,4), denom=9, weight=FALSE)$LR
ilr1.sub <- ILR(aar.sub, numer=c(1,6), denom=8, weight=FALSE)$LR
par(mar=c(4.2,4,1,2), font.lab=2, cex.lab=1.3, mgp=c(3,0.7,0), mfrow=c(1,3), las=1)
plot(ilr1.sub, ilr1, xlab="ILR of {SiO2, Na2O}:P2O5", ylab="ILR of {SiO2, Na2O, MnO}:P2O5")

gm1 <- (aar[,1]*aar[,7]*aar[,4])^(1/3)
gm1.sub <- (aar.sub[,1]*aar.sub[,6])^(1/2)
plot(gm1.sub, gm1, xlab="GM of {SiO2, Na2O}", ylab="GM of {SiO2, Na2O, MnO}")

slr1 <- SLR(aar, numer=c(1,7,4), denom=9, weight=FALSE)$LR
slr1.sub <- SLR(aar.sub, numer=c(1,6), denom=8, weight=FALSE)$LR
plot(slr1.sub, slr1, xlab="SLR of {SiO2, Na2O}:P2O5", ylab="SLR of {SiO2, Na2O, MnO}:P2O5")

### Stepwise selection of ratios, including the three amalgamations

#   Define the amalgamations and add them to the set of 10 parts
mafic      <- apply(aar[,c(5,10,4)], 1, sum)
felsic     <- apply(aar[,c(7,1,3,8)], 1, sum)
carbonate  <- apply(aar[,c(6,9)], 1, sum)
aar.amalg  <- cbind(aar, mafic, felsic, carbonate)

#   Perform the stepwise analysis
aar.step   <- STEP(aar.amalg, aar, weight=FALSE)

# Table 1
cbind(aar.step$ratios, 
      round(100*aar.step$R2max,1), 
      round(aar.step$pro.cor,3))

### LRA of original data and PCA of the 9 selected logratios

rownames(aar) <- 1:nrow(aar)
par(mar=c(4.2,4,3,3), font.lab=2, cex.lab=1.2, mgp=c(2.7,0.7,0), las=1, mfrow=c(1,2))

### LRA of original data and PCA of the 9 selected logratios

rownames(aar) <- 1:nrow(aar)
par(mar=c(4.2,4,3,3), font.lab=2, cex.lab=1.2, mgp=c(2.7,0.7,0), las=1, mfrow=c(1,2), cex.axis=0.8)

#   LRA (logratio analysis, = PCA of the CLRs)
aar.lra <- LRA(aar, weight=FALSE)
PLOT.LRA(aar.lra, map="contribution")

#   PCA of the selected logratios
rownames(aar.step$logratios) <- 1:nrow(aar)
# invert K2O/P2O5
aar.step$logratios[,2] <- -aar.step$logratios[,2]
colnames(aar.step$logratios)[2] <- "P2O5/K2O"
aar.ratios.pca <- PCA(aar.step$logratios, weight=FALSE)
PLOT.PCA(aar.ratios.pca, map="contribution", axes.inv=c(1,-1), rescale=2)

#   Procrustes correlation of two configurations of samples in two dimensions
protest(aar.ratios.pca$rowpcoord[,1:2], aar.lra$rowpcoord[,1:2])$t0
# [1] 0.9971076

#   Procrustes correlation of full-space geometry of samples
protest(aar.ratios.pca$rowpcoord, aar.lra$rowpcoord)$t0
# [1] 0.993197
