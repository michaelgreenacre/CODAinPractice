### Script for analysis of single-cell gene data
### Set working directory for loading the data objects
load("SingleCell.RData")
ls()
 [1] "allcols"  "cellcols" "cells"    "ctypes"   "M"        "Mi"       "Mj"       "myind"   
 [9] "subcomp"  "supind"   "types"  

### Description of each object
# M		raw count matrix of all pooled samples
dim(M)  # [1]  724 6147
# Mi		like M but containing 0-imputed compositions
dim(Mi) # [1]  724 6147
# Mj		joint matrix of active samples from M and raw counts of all single cells  45 active + 12611 passive
dim(Mj) # [1] 12656  6147
# allcols	"colors" (numbers 1 to 5) for all 724 pooled samples
length(allcols)  # [1] 724
# cellcols	"colors" (numbers 1 to 6) for all single cells
length(cellcols)  # [1] 12611
cells		cell indices (of rows of matrix Mj)
length(cells) # [1] 12611
# ctypes      cell type labels for all cells, including NA
table(ctypes)
#   Antigen presenting              Stromal    Thymic epithelial           Thymocytes Vascular endothelial 
#                  465                  112                  259                 7852                   91 
sum(is.na(ctypes))
# [1] 3832
sum(table(ctypes))+ sum(is.na(ctypes))
# [1] 12611

### Some preliminary data re-organization
ctypes[is.na(ctypes)] <- "Unknown"
table(ctypes)
#   Antigen presenting              Stromal    Thymic epithelial           Thymocytes              Unknown 
#                  465                  112                  259                 7852                 3832 
# Vascular endothelial 
#                   91 

ctypes.ind <- as.numeric(as.factor(ctypes))
table(ctypes.ind)
#    1    2    3    4    5    6 
#  465  112  259 7852 3832   91 
ctypes.ind[ctypes.ind==5] <- 7
ctypes.ind[ctypes.ind==6] <- 5
ctypes.ind[ctypes.ind==7] <- 6
table(ctypes.ind)
#    1    2    3    4    5    6 
#  465  112  259 7852   91 3832 

length(ctypes)  # [1] 12611
# myind		active sample indices of rows of M and Mi
myind
#  [1]   1   2   3  72  73 542 543 546 547 651 652 653 654 655 658 659 660 681 682 683 689 690 691 692 
# [25] 693 694 695 696 708 709 710 711 712 713 714 715 716 717 718 719 720 721 722 723 724
# subcomp		gene labels after variable selection
length(subcomp)
# [1] 1402
# supind		supplementary sample indices of rows of M and Mi
length(supind)  # [1] 679
# types		cell type labels for pooled samples
length(types)   # [1] 724

### Percentage of missing data in whole matrix and subcomposition
100*sum(M==0)/(nrow(M)*ncol(M))
# [1] 64.44088
100*sum(M[,subcomp]==0)/(nrow(M[,subcomp])*ncol(M[,subcomp]))
# [1] 76.09601

### LRA and %s of variance
require(easyCODA)
lra <- LRA(CLOSE(Mi), suprow=supind)
round(100*lra$sv^2/sum(lra$sv^2),2)[1:10]
[1] 25.31  9.45  4.35  4.19  4.13  3.84  3.79  3.61  3.47  3.04
rownames(lra$colcoord) <- colnames(Mi)
require(RColorBrewer)
genecols <- brewer.pal(7, "Dark2")[c(1,6,4,3,7)]
plot(lra$colcoord[,1]*lra$sv[1],lra$colcoord[,2]*lra$sv[2],asp=1,pch=20,col="lightgrey",cex=0.2, 
     xlab="LRA dimension 1 (25.3%)", ylab="LRA dimension 2 (9.5%)",
     main="Gene-principal LRA biplot", xlim=c(-1.2,2.4), ylim=c(-1.5,2))
points(lra$colcoord[subcomp,1]*lra$sv[1],lra$colcoord[subcomp,2]*lra$sv[2],asp=1,pch=20,col="skyblue1",cex=0.7)
points(lra$rowcoord[myind,1],lra$rowcoord[myind,2],col=genecols[allcols[myind]],asp=1,pch=20)
legend("bottomright", legend=c("Thymocytes","Antigen presenting","Thymic epithelial","Stromal","Vascular endothelial"), 
       pch=19, col=genecols, text.col=genecols, pt.cex=0.6, cex=0.7, text.font=2)

plot(lra$colcoord[,3]*lra$sv[3],lra$colcoord[,4]*lra$sv[4],asp=1,pch=20,col="lightgrey",cex=0.2, 
     xlab="LRA dimension 3 (4.4%)", ylab="LRA dimension 2 (4.2)",
     main="Gene-principal LRA biplot", xlim=c(-1.1,2), ylim=c(-1.2,1.3))
points(lra$colcoord[subcomp,3]*lra$sv[3],lra$colcoord[subcomp,4]*lra$sv[4],asp=1,pch=20,col="skyblue1",cex=0.7)
points(lra$rowcoord[myind,3],lra$rowcoord[myind,4],col=genecols[allcols[myind]],asp=1,pch=20)
legend("right", legend=c("Thymocytes","Antigen presenting","Thymic epithelial","Stromal","Vascular endothelial"), 
       pch=19, col=genecols, text.col=genecols, pt.cex=0.6, cex=0.7, text.font=2)


### CA and %s of inertia
rca=ca(CLOSE(M),suprow=supind)
round(100*rca$sv^2/sum(rca$sv^2),2)[1:10]
#  [1] 12.01  9.37  5.86  5.02  2.50  2.30  2.16  2.11  2.04  2.03

genecols <- brewer.pal(7, "Dark2")[c(1,6,4,3,7)]
plot(rca$colcoord[,1]*rca$sv[1],-rca$colcoord[,2]*rca$sv[2],asp=1,pch=20,col="lightgrey",cex=0.2, 
     xlab="CA dimension 1 (12.0%)", ylab="CA dimension 2 (9.4%)",
     main="Gene-principal CA biplot", xlim=c(-1.2,2.4), ylim=c(-1.5,1.8))
points(rca$colcoord[subcomp,1]*rca$sv[1],-rca$colcoord[subcomp,2]*rca$sv[2],asp=1,pch=20,col="skyblue1",cex=0.7)
points(rca$rowcoord[myind,1],-rca$rowcoord[myind,2],col=genecols[allcols[myind]],asp=1,pch=20)
legend("bottomright", legend=c("Thymocytes","Antigen presenting","Thymic epithelial","Stromal","Vascular endothelial"), 
       pch=19, col=genecols, text.col=genecols, pt.cex=0.6, cex=0.7, text.font=2)

plot(rca$colcoord[,3]*rca$sv[3],rca$colcoord[,4]*rca$sv[4],asp=1,pch=20,col="lightgrey",cex=0.2, 
     xlab="CA dimension 3 (5.9%)", ylab="CA dimension 4 (5.0%)",
     main="Gene-principal CA biplot", xlim=c(-2.6,1.5), ylim=c(-1.9,2.1))
abline(h = 0, v = 0, col = "gray", lty = 2)
points(rca$colcoord[subcomp,3]*rca$sv[3],rca$colcoord[subcomp,4]*rca$sv[4],asp=1,pch=20,col="skyblue1",cex=0.7)
points(rca$rowcoord[myind,3],rca$rowcoord[myind,4],col=genecols[allcols[myind]],asp=1,pch=20)
legend("topleft", legend=c("Thymocytes","Antigen presenting","Thymic epithelial","Stromal","Vascular endothelial"), 
       pch=19, col=genecols, text.col=genecols, pt.cex=0.6, cex=0.7, text.font=2)


### CA coherence exercise, looking at subcompositions of 6147 genes 
### (subcompositions in same proportions as Tellus study)
set.seed(1234567)
genes.pro <- CLOSE(M[myind,])
genes.CA.coherence <- matrix(0, 100, 11) 
genes.CA.cpc <- CA(genes.pro)$colpcoord
k <- 1
for(j in seq(4,44,4)) {
  nparts <- round((j/52)*6147)
  for(i in 1:100) {
# find the subcompositional parts  
    jparts <- sample(1:6147, nparts)
    foo.pro <- CLOSE(genes.pro[,jparts])
# remove samples all zeros
    allzero <- which(rowSums(foo.pro)==0)
    if(length(allzero)>0) foo.pro <- foo.pro[-allzero,]
    genes.foo.cpc <- CA(foo.pro)$colpcoord
    genes.CA.coherence [i,k] <- protest(genes.CA.cpc[jparts,], genes.foo.cpc, permutations=0)$t0
  }
  k <- k+1
}

genes.CA.quants <- apply(genes.CA.coherence, 2, quantile, c(0.025,0.5,0.975), na.rm=TRUE)
round(genes.CA.quants,4)
        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]  [,10]  [,11]
2.5%  0.9867 0.9943 0.9957 0.9963 0.9975 0.9978 0.9983 0.9990 0.9987 0.9992 0.9997
50%   0.9954 0.9975 0.9982 0.9989 0.9993 0.9995 0.9996 0.9997 0.9998 0.9998 0.9999
97.5% 0.9980 0.9991 0.9996 0.9996 0.9997 0.9998 0.9999 0.9999 0.9999 1.0000 1.0000
genes.CA.ones <- rep(0,11)
for(j in 1:11) genes.CA.ones[j] <- sum(genes.CA.coherence[,j]>0.999)
genes.CA.ones
[1]  0  9 26 47 64 78 85 97 95 99 99

### Figure 17
# pdf(file="SingleCell_CA_coherence.pdf", width=7.5, height=4, useDingbats=FALSE, family="ArialMT")
par(mar=c(5,5,1,1), mgp=c(3.5,0.7,0), font.lab=2, las=1, mfrow=c(1,1))
plot(rep(1:11, each=3), as.numeric(genes.CA.quants), xlab="Number of parts in subcomposition",
     ylab="Procrustes correlation", bty="n", xaxt="n", ylim=c(0.98, 1.001), type="n", font.lab=2, xlim=c(1,11))
axis(1, at=1:11, labels=round((seq(4,44,4)/52)*6147))
for(j in 1:11) segments(j, genes.CA.quants[1,j], j, genes.CA.quants[3,j], col="blue", lwd=2)
eps <- 0.06
for(j in 1:11) segments(j-eps, genes.CA.quants[1,j], j+eps, genes.CA.quants[1,j], col="blue", lwd=2, lend=2)
for(j in 1:11) segments(j-eps, genes.CA.quants[3,j], j+eps, genes.CA.quants[3,j], col="blue", lwd=2, lend=2)
points(1:11,genes.CA.quants[2,], pch=21, col="blue", bg="white", cex=0.9)
text(1:11, rep(1.001, 11), labels=genes.CA.ones, font=2, cex=0.8)
# dev.off()


### CA coherence study for compositions of different sizes
genes.CA.comp <- matrix(0,100,9)
set.seed(1234567)
for(j in 1:9) {
  nparts <- round(6147*j/10)
  for(k in 1:100) {
    foo <- M[myind,sample(1:6147, nparts)]
# remove samples all zeros
    allzero <- which(rowSums(foo)==0)
    if(length(allzero)>0) foo <- foo[-allzero,]
### 20% sample in subcomposition
    subsample <- sample(1:nparts, round(nparts/5))
    foo.sub <- foo[,subsample]
# remove samples all zeros
    allzero <- which(rowSums(foo.sub)==0)
    if(length(allzero)>0) foo.sub <- foo.sub[-allzero,] 
    foo.pro <- CLOSE(foo)
    foo.sub <- CLOSE(foo.sub)
    foo.cpc <- CA(foo)$colpcoord[subsample,]
    foo.sub.cpc <- CA(foo.sub)$colpcoord
    genes.CA.comp[k,j] <- protest(foo.cpc, foo.sub.cpc, permutations=0)$t0  
  } 
}
genes.CA.comp.quant <- apply(genes.CA.comp, 2, quantile, c(0.025,0.5,1))  
round(genes.CA.comp.quant, 3)
      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9]
2.5% 0.979 0.983 0.985 0.988 0.990 0.988 0.992 0.992 0.993
50%  0.993 0.995 0.995 0.996 0.997 0.997 0.998 0.998 0.998
100% 0.997 0.998 0.999 0.999 0.999 0.999 0.999 0.999 1.000

nparts <- round(6147*(1:9)/10)
genes.CA.comp.ones <- rep(0,9)
for(j in 1:9) genes.CA.comp.ones[j] <- sum(genes.CA.comp[,j]>0.999)
genes.CA.comp.ones
[1]  0  0  0  0  0  0  6  9 23

### for square-rooted data
genes.CA.comp05 <- matrix(0,100,9)
set.seed(1234567)
for(j in 1:9) {
  nparts <- round(6147*j/10)
  for(k in 1:100) {
    foo <- M[myind,sample(1:6147, nparts)]
### 20% sample in subcomposition
    subsample <- sample(1:nparts, round(nparts/5))
    foo.sub <- foo[,subsample]
    foo <- CLOSE(foo^0.5)
    foo.sub <- CLOSE(foo.sub^0.5)
    foo.cpc <- CA(foo)$colpcoord[subsample,]
    foo.sub.cpc <- CA(foo.sub)$colpcoord
    genes.CA.comp05[k,j] <- protest(foo.cpc, foo.sub.cpc, permutations=0)$t0  
  } 
}
genes.CA.comp.quant05 <- apply(genes.CA.comp05, 2, quantile, c(0.025,0.5,1))  
round(genes.CA.comp.quant05, 3)
      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9]
2.5% 0.994 0.995 0.997 0.998 0.998 0.999 0.999 0.999 0.999
50%  0.996 0.998 0.998 0.999 0.999 0.999 0.999 0.999 0.999
100% 0.998 0.999 0.999 0.999 0.999 1.000 1.000 1.000 1.000

nparts <- round(6147*(1:9)/10)
genes.CA.comp.ones05 <- rep(0,9)
for(j in 1:9) genes.CA.comp.ones05[j] <- sum(genes.CA.comp05[,j]>0.999)
genes.CA.comp.ones05
[1]  0  0  4 25 58 72 89 96 94

### Figure 18
# pdf(file="SingleCell_CA_comp.pdf", width=7.5, height=8, useDingbats=FALSE, family="ArialMT")
par(mar=c(5,5,1,1), mgp=c(2.5,0.7,0), font.lab=2, las=1, mfrow=c(2,1), cex.axis=0.9)
plot(rep(seq(10,90,10), each=3), as.numeric(genes.CA.comp.quant), xlab="Numbers of parts in composition",
     ylab="Procrustes correlation (coherence)", bty="n", xaxt="n", ylim=c(0.975, 1.002), type="n", font.lab=2,
     main="20% subcompositions of compositions of increasing sizes")
axis(1, at=seq(10,90,10), labels=round(6147*(1:9)/10))
for(j in 1:9) segments(10*j, genes.CA.comp.quant[1,j], 10*j, genes.CA.comp.quant[3,j], col="blue", lwd=2)
eps <- 0.5
for(j in 1:9) segments(10*j-eps, genes.CA.comp.quant[1,j], 10*j+eps, genes.CA.comp.quant[1,j], col="blue", lwd=2, lend=2)
for(j in 1:9) segments(10*j-eps, genes.CA.comp.quant[3,j], 10*j+eps, genes.CA.comp.quant[3,j], col="blue", lwd=2, lend=2)
points(seq(10,90,10), genes.CA.comp.quant[2,], pch=21, col="blue", bg="white")
text(seq(10,90,10), rep(1.001, 11), labels=genes.CA.comp.ones, font=2, cex=0.8)

plot(rep(seq(10,90,10), each=3), as.numeric(genes.CA.comp.quant05), xlab="Numbers of parts in composition",
     ylab="Procrustes correlation (coherence)", bty="n", xaxt="n", ylim=c(0.975, 1.002), type="n", font.lab=2,
     main="20% subcompositions of compositions^0.5 of increasing sizes")
axis(1, at=seq(10,90,10), labels=round(6147*(1:9)/10))
for(j in 1:9) segments(10*j, genes.CA.comp.quant05[1,j], 10*j, genes.CA.comp.quant05[3,j], col="blue", lwd=2)
eps <- 0.5
for(j in 1:9) segments(10*j-eps, genes.CA.comp.quant05[1,j], 10*j+eps, genes.CA.comp.quant05[1,j], col="blue", lwd=2, lend=2)
for(j in 1:9) segments(10*j-eps, genes.CA.comp.quant05[3,j], 10*j+eps, genes.CA.comp.quant05[3,j], col="blue", lwd=2, lend=2)
points(seq(10,90,10), genes.CA.comp.quant05[2,], pch=21, col="blue", bg="white")
text(seq(10,90,10), rep(1.001, 11), labels=genes.CA.comp.ones05, font=2, cex=0.8)
# dev.off()
  

### now the 12000+ single cells
sca <- ca(CLOSE(Mj), suprow=cells)
sca.rsc <- sca$rowcoord
dim(sca.rsc)
[1] 12656    44
round(100*sca$sv^2/sum(sca$sv^2),2)
[1]  12.01  9.37  5.86  5.02  2.50  2.30  2.16  2.11  2.04  2.03  1.99  1.98  1.95  1.94  1.92  1.87  1.85  1.82
[19]  1.80  1.79  1.78  1.76  1.73  1.71  1.69  1.68  1.66  1.64  1.61  1.59  1.57  1.50  1.47  1.46  1.43  1.41
[37]  1.38  1.33  1.32  1.29  1.24  1.22  1.18  1.04

### Figure S3
par(mar=c(4.2,4,4,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.6, mfrow=c(1,2))
plot(sca.rsc[cells,1],-sca.rsc[cells,2],asp=1,pch=20,col=c(genecols,"yellow")[cellcols],cex=0.2, 
     xlab="CA dimension 1 (12.0%)", ylab="CA dimension 2 (9.4%)",
     main="12611 single cells as supplementary points")
abline(h = 0, v = 0, col = "gray", lty = 2)
legend("bottomright", legend=c("Thymocytes","Antigen presenting","Thymic epithelial","Stromal","Vascular endothelial","Unknown"), 
       pch=21, col=c(genecols,"yellow"), text.col=c(genecols,"gray60"), pt.bg=c(genecols,"yellow"), pt.cex=0.4, cex=0.6, text.font=2)

plot(sca.rsc[cells,3],sca.rsc[cells,4],asp=1,pch=20,col=c(genecols,"yellow")[cellcols],cex=0.2, 
     xlab="CA dimension 3 (5.9%)", ylab="CA dimension 4 (5.0%)",
     main="12611 single cells as supplementary points")
abline(h = 0, v = 0, col = "gray", lty = 2)
legend("bottomleft", legend=c("Thymocytes","Antigen presenting","Thymic epithelial","Stromal","Vascular endothelial","Unknown"), 
       pch=21, col=c(genecols,"yellow"), text.col=c(genecols,"gray60"), pt.bg=c(genecols,"yellow"), pt.cex=0.4, cex=0.6, text.font=2)



















points(rca$rowcoord[supind,1],rca$rowcoord[supind,2],col=allcols[supind])

points(rca$rowcoord[myind,1]*rca$sv[1],-rca$rowcoord[myind,2]*rca$sv[2],col=allcols[myind],pch=20)
points(-rca$rowcoord[myind,3]*rca$sv[3],-rca$rowcoord[myind,4]*rca$sv[4],col=allcols[myind],pch=20)
abline(h=0, v=0, col="grey", lty=2)

#subcomposition & reference gene MALAT:
points(rca$colcoord[,1]*rca$sv[1],-rca$colcoord[,2]*rca$sv[2],asp=1,pch=20,col="lightgrey",cex=0.3)#, xlab="CA dim1", ylab="CA dim2",
     main="Row-principal biplot")
plot(-rca$colcoord[,3],-rca$colcoord[,4],asp=1,pch=20,col="lightgrey",cex=0.3, xlab="CA dim3", ylab="CA dim4",
     main="Row-principal biplot")
plot(5*rca$colcoord[,1]*sqrt(rca$colmass),-5*rca$colcoord[,2]*sqrt(rca$colmass),asp=1,pch=20,col="lightgrey",cex=0.5, xlab="CA dim1", ylab="CA dim2",
     main="Row-principal contribution biplot")

which(5*rca$colcoord[,1]*sqrt(rca$colmass)>1.5)
ENSG00000164692.13  ENSG00000108821.9  ENSG00000168542.8 
                94                183                318 

CLOSE(M[myind,subcomp])[,c(94,183,318)]

points(rca$colcoord[subcomp,1],-rca$colcoord[subcomp,2],pch=20,col="slateblue",cex=0.3)
points(-rca$colcoord["ENSG00000251562.3",1],rca$colcoord["ENSG00000251562.3",2],col="orange",pch=20)


round(rca$sv^2/ sum(rca$sv^2),4)[1:10]
 [1] 0.1959 0.1294 0.0839 0.0600 0.0257 0.0233 0.0224 0.0218 0.0211 0.0206

# CA Ionas style
rescale <- 5
dim <- c(1,2)
col <- c("blue","red")
perc.hor <- 19.6; perc.ver <- 12.9
par(mar=c(4.2,4,4,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.6, mfrow=c(1,1))
plot(rescale*rca$colcoord[,1]*sqrt(rca$colmass),-rescale*rca$colcoord[,2]*sqrt(rca$colmass), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "LRA of 1402 genes, 45 active samples")
# abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*rca$colcoord[,1]*sqrt(rca$colmass),-rescale*rca$colcoord[,2]*sqrt(rca$colmass), pch = 19, col = gray(genes3.ca.ctr), cex = 0.5)
#text(lras.rpc[myind,], labels = substr(lras$rownames[myind],1,1), col =allcols[myind], cex = 0.6, font = 2) 
points(rca$rowcoord[myind,1]*rca$sv[1],-rca$rowcoord[myind,2]*rca$sv[2], pch=20, col =allcols[myind],  font = 2) 


points(lras.rpc[supind,], col =allcols[supind], cex=0.9) 

legend("topright", legend=c("Antigen presenting","Stromal","Thymic epithelial","Thymocytes","Vascular endothelial"), 
       pch=19, col=1:5, text.col=1:5, pt.cex=0.6, cex=0.7, text.font=2)




genes3.ca <- CA(CLOSE(M[,subcomp]), suprow=supind)

genes3.ca.rpc <- genes3.ca$rowpcoord 
genes3.ca.csc <- genes3.ca$colcoord
genes3.ca.ccc <- genes3.ca$colcoord * sqrt(genes3.ca$colmass)
genes3.ca.ctr <- genes3.ca.ccc[,1]^2 + genes3.ca.ccc[,2]^2 
genes3.ca.ctr <- 1-genes3.ca.ctr^0.5/max(genes2.ca.ctr^0.5) 
ca.grey <- round(genes3.ca.ctr*80)+1
round(100*genes3.ca$sv^2/sum(genes3.ca$sv^2),2)[1:20]
 [1] 12.65  9.48  4.05  3.40  0.84  0.53  0.49  0.47  0.45  0.43  0.42  0.38  0.35  0.34  0.33  0.32  0.32  0.32
[19]  0.31  0.30



### invert 2nd axis
genes3.ca.rpc[,2] <- -genes3.ca.rpc[,2] 
genes3.ca.ccc[,2] <- -genes3.ca.ccc[,2]

rescale <- 0.5
dim <- c(1,2)
col <- c("blue","red")
perc.hor <- 6.4; perc.ver <- 5.4
par(mar=c(4.2,4,4,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, mfrow=c(1,1))
plot(1.05 * rbind(genes3.ca.rpc, rescale*genes3.ca.csc), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "CA of 1402 genes, with zeros, 45 active samples")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*genes3.ca.csc, pch = 19, col = gray(genes3.ca.ctr), cex = 0.5)
text(genes3.ca.rpc[myind,], labels = substr(genes3.ca$rownames[myind],1,1), col =allcols[myind], cex = 0.7, font = 2) 

legend("bottomright", legend=ABnames, 
       pch=tellus.pch, col=tellus.col, text.col=tellus.col, pt.cex=0.6, cex=0.8)


### Power study
### Now with power transformation and with original data zeros
### target.lra.rpc <- lras.rpc[myind,] is the active set target
target.lra.rpc <- lras.rpc[myind,]
genes3.pro <- CLOSE(M[myind,subcomp])
BoxCox3.genes <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- genes3.pro^alpha
  foo.ca <- CA(CLOSE(foo))
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox3.genes[k] <- protest(target.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

### What is maximum Procrustes correlation
max(BoxCox3.genes)
# [1] 0.96288

### For which power?
c(seq(1,0.01,-0.01),0.0001)[which(BoxCox3.genes == max(BoxCox3.genes))]
# [1] 0.5

#CA figure for all cells:

bca=ca(Mj,suprow=cells)

plot(bca$rowcoord[cells,1],bca$rowcoord[cells,2],asp=1,col=cellcols,pch=20,cex=0.3)
legend("topright",c(unique(types),"unassigned"),text.col=c(1:6))



# some example command lines:



colors=as.numeric(substr(rownames(Mt),1,1))

The unimputed matrix "Mt" I use for CA, the imputed one "Mi" I use for LRA. 
(I multiply the coordinates of the 2nd CA dimension by -1 to make the display more similar to the LRA.)

By the way, I quickly checked and your intuition is right: 
the subcomposition hardly makes a difference to the CA! 
This makes the variable selection appear a bit silly though... 
It looks like choosing the highly expressed genes already separates the groups extremely well, 
and the ANOVA is not really necessary for visualisation. 
Perhaps I should start off with more than 3596 genes. 
(Before I tried with much more genes, but I used genes with counts that were too small, 
and the ANOVA results didn't make any sense.)

I use the column-principal map because I like the samples in the corners of the simplex. 
The genes are then sitting between the corners according to their expression, 
which I find intuitive. 
I think the LRA also looked prettier because the variance of columns and rows is similar, 
and you don't need different axes for columns and rows. 
I also checked the row-principal CA, and I don't think the weighting made a big difference to the plot. (I don't insist on this, of course, I understand that row-principal seems more logical here.)

Oh, and I didn't explore the optimal power transformation for the CA. 
I did try the square root though, and it rather made the plot look worse, 
so I didn't want to complicate things.

Now that Thom seems to be available again, 
I wonder whether we should also consider adding an alternative variable selection using CODACORE? 
(I haven't discussed this with him yet.)

genes0 <- Mt
dim(genes0)
[1]   24 3569
genes0[,1:10]
types <- as.numeric(substr(rownames(Mt),1,1))
sum(genes0==0) / (nrow(genes0)*ncol(genes0))
[1] 0.2825721

genes0.size <- rowSums(genes0)
  1_1   1_2   1_3   1_4   1_5   1_6   2_1   2_2   2_3   2_4   2_5   2_6   3_1   3_2   3_3   3_4   3_5   3_6   4_1 
 8975 13120  9368  9269  9444 11320 11760 11058  9893  8033  6733  7555  9727  9133  8627  9775  9328 11054  8532 
  4_2   4_3   5_1   5_2   5_3 
 7498 13382  9291 12066  9727 

model0 <- Mi
Mi[,1:10]
rowSums(Mi)
1_1 1_2 1_3 1_4 1_5 1_6 2_1 2_2 2_3 2_4 2_5 2_6 3_1 3_2 3_3 3_4 3_5 3_6 4_1 4_2 4_3 5_1 5_2 5_3 
  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
sum(Mi==0)
[1] 0

require(easyCODA)
genes0.ca <- CA(CLOSE(genes0))
model0.ca <- CA(model0)

genes0.rpc <- genes0.ca$rowpcoord
model0.rpc <- model0.ca$rowpcoord
genes0.cpc <- genes0.ca$colpcoord
model0.cpc <- model0.ca$colpcoord

protest(genes0.rpc, model0.rpc, permutations=0)$t0
[1] 0.9982627
protest(genes0.cpc, model0.cpc, permutations=0)$t0
[1] 0.9963113

model0.lra <- LRA(model0, weight=FALSE)	
model0.lra.rpc <- model0.lra$rowpcoord
model0.lra.cpc <- model0.lra$colpcoord

protest(model0.rpc, model0.lra.rpc, permutations=0)$t0
[1] 0.9581373

protest(model0.cpc, model0.lra.cpc, permutations=0)$t0
[1] 0.8275609

### ----------------------------------------------------------------------
### ----------------------------------------------------------------------
### study of the Box-Cox transformation in CA on the geometry of the samples
### 
model.lra <- LRA(model0, weight=FALSE)
model.lra.rpc <- model.lra$rowpcoord 

### for original CA/chi-square, no zero replacement
genes0.ca <- CA(CLOSE(genes0))
genes0.ca.rpc <- genes0.ca$rowpcoord 
protest(model.lra.rpc, genes0.ca.rpc, permutations=0)$t0
# [1] 0.9557203

### Now CA/chi-square with power transformation on data with zeros replaced
### Sequence of powers down to the smallest 0.0001 
BoxCox.genes <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- model0^alpha
  foo.ca <- CA(foo)
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox.genes[k] <- protest(model.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

### Now with power transformation and with original data zeros
genes0.pro <- CLOSE(genes0)
BoxCox0.genes <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- genes0.pro^alpha
  foo.ca <- CA(CLOSE(foo))
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox0.genes[k] <- protest(model.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

### What is maximum Procrustes correlation
max(BoxCox0.genes)
# [1] 0.9705153


### For which power?
c(seq(1,0.01,-0.01),0.0001)[which(BoxCox0.genes == max(BoxCox0.genes))]
# [1] 0.41

### -----------------------------------------------------------------------
### plots of Procrustes correlations for Box-Cox transformation 
###            CA for data with zeros replaced and data with original zeros
# pdf(file="Fig_X.pdf", width=5, height=5, useDingbats=FALSE, family="ArialMT")

### use square window
par(mar=c(4.2,4,1,1), font.lab=2, las=1)
plot(c(seq(1,0.01,-0.01),0.0001), BoxCox.genes, xlab="Power of Box-Cox transformation", 
     ylab="Procrustes correlation", type="l", lwd=2, col="blue", ylim=c(0.95,1),
     bty="n", xaxt="n", yaxt="n")
axis(1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
axis(2)
lines(c(seq(1,0.01,-0.01),0.0001), BoxCox0.genes, lwd=2, col="red", lty=3)
segments(0.41,0,0.41,BoxCox0.genes[60], col="pink", lwd=2, lty=2)
legend("bottomright", legend=c("zeros replaced","with zeros"), 
       bty="n",
       col=c("blue","red"), 
       lwd=c(2,2), lty=c(1,3), cex=0.8)

# dev.off()

### ---------------------------------------------------------------------------------
### should we weight the parts?
genes.clr.unw <- CLR(model0, weight=FALSE)$LR
genes.clr.unw.var <- apply(genes.clr.unw, 2, var)
genes.pro.cm <- colMeans(model0)
# pdf(file="Fig_1.pdf", width=5, height=5, useDingbats=FALSE, family="ArialMT")
### use square window for this plot
par(mar=c(4.2,4,1,1), font.lab=2, las=1, mfrow=c(1,1))
plot(genes.pro.cm, genes.clr.unw.var, # log="xy", 
     xlab="Average compositional values (log-scale)", ylab="Variance of CLR (log-scale)")
### seems like unweighted OK: no tendency for rarer elements to have much higher variance
### in fact Al, the second highest component, has the lowest CLR variance
### ---------------------------------------------------------------------------------


(model.findalr <- FINDALR(model0))

$tot.var
[1] 1.39058
...
$procrust.max
[1] 0.9957712

$procrust.ref
[1] 2138
...

$var.min
[1] 0.08112942

$var.ref
[1] 2138

colnames(model0)[2138]


which(genes.pro.cm == max(genes.pro.cm))
ENSG00000251562.3 
             2138

apply(log(model0), 2, sd)[2138]
ENSG00000251562.3 
        0.2848323 
apply(log(model0), 2, range)[,2138]
[1] -5.088666 -3.844784


#####################################################################################
### some plots

require(RColorBrewer)
type.cols <- brewer.pal(5, "Dark2")


### LRA
model.lra.rpc <- -model.lra$rowpcoord 
model.lra.ccc <- -model.lra$colcoord * sqrt(model.lra$colmass)
model.lra.ctr <- model.lra.ccc[,1]^2 + model.lra.ccc[,2]^2 
model.lra.ctr <- 1-model.lra.ctr/max(model.lra.ctr) 

model.grey <- round(model.lra.ctr*80)+1
100*model.lra$sv^2/sum(model.lra$sv^2)

rescale <- 20
dim <- c(1,2)
col <- c("blue","red")
perc.hor <- 20.5; perc.ver <- 15.9
par(mar=c(4.2,4,4,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, mfrow=c(1,3))
plot(1.05 * rbind(model.lra.rpc, rescale*model.lra.ccc), type = "n", asp = 1, 
     xlab = paste("LRA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("LRA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "LRA of 3569 genes")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*model.lra.ccc, pch = 19, col = gray(model.lra.ctr), cex = 0.5)
text(model.lra.rpc, labels = types, col = type.cols[types], cex = 0.7, font = 2) 

   
legend("bottomleft", legend=ABnames, 
       pch=tellus.pch, col=tellus.col, text.col=tellus.col, pt.cex=0.6, cex=0.8)


### ALR
model.alr <- PCA(ALR(model0, denom=2138, weight=FALSE)$LR, weight=FALSE)
model.alr.rpc <- -model.alr$rowpcoord 
model.alr.ccc <- -model.alr$colcoord * sqrt(model.alr$colmass)
100*model.alr$sv^2/sum(model.alr$sv^2)
model.alr.ctr <- model.alr.ccc[,1]^2 + model.alr.ccc[,2]^2 
model.alr.ctr <- 1-model.alr.ctr/max(model.alr.ctr) 

rescale <- 20
dim <- c(1,2)
col <- c("blue","red")
perc.hor <- 20.4; perc.ver <- 18.1

plot(1.05 * rbind(model.alr.rpc, rescale*model.alr.ccc), type = "n", asp = 1, 
     xlab = paste("PCA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("PCA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "PCA of 3568 ALRs w.r.t. reference #2138")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*model.alr.ccc, pch = 19, col = gray(model.alr.ctr), cex = 0.5)
text(model.alr.rpc, labels = types, col = type.cols[types], cex = 0.7, font = 2) 

   
### CA with 0.5 power transformation
genes05.ca <- CA(CLOSE(genes0^0.5))
genes05.rpc <- -genes05.ca$rowpcoord/0.5
genes05.ccc <- -genes05.ca$colcoord * sqrt(genes05.ca$colmass)
genes05.ctr <- genes05.ccc[,1]^2 + genes05.ccc[,2]^2 
genes05.ctr <- 1-sqrt(genes05.ctr)/max(sqrt(genes05.ctr)) 
100*genes05.ca$sv^2/sum(genes05.ca$sv^2)
[1] 22.093012 14.848020  8.329700  8.191394  3.483925  3.137939  2.981124  2.877830  2.695056  2.642760  2.586955
[12]  2.528250  2.479204  2.422827  2.363432  2.286030  2.233785  2.189060  2.134450  2.032936  1.997756  1.805325
[23]  1.659232

rescale <- 15
dim <- c(1,2)
col <- c("blue","red")
perc.hor <- 22.1; perc.ver <- 14.8

plot(1.05 * rbind(genes05.rpc, rescale*genes05.ccc), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "CA of square-root compositions (with zeros)",
     ylim=c(-1,3), xlim=c(-3,1))
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*genes05.ccc, pch = 19, col = gray(genes05.ctr), cex = 0.5)
text(genes05.rpc, labels = types, col = type.cols[types], cex = 0.7, font = 2) 





### some plots, dims 3 and 4

### LRA
model.lra.rpc <- -model.lra$rowpcoord 
model.lra.ccc <- -model.lra$colcoord * sqrt(model.lra$colmass)
100*model.lra$sv^2/sum(model.lra$sv^2)
[1] 20.471480 15.902407  9.676252  7.758739  5.541147  5.038894  3.174871  2.970315  2.828779  2.662129  2.621921
[12]  2.488058  2.430173  2.116851  1.793059  1.778180  1.707557  1.668930  1.578414  1.538921  1.507372  1.450327
[23]  1.295223

rescale <- 20
dim <- c(3,4)
col <- c("blue","red")
perc.hor <- 9.7; perc.ver <- 7.8
par(mar=c(4.2,4,4,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, mfrow=c(1,3))
plot(1.05 * rbind(model.lra.rpc[,dim], rescale*model.lra.ccc[,dim]), type = "n", asp = 1, 
     xlab = paste("LRA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("LRA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "LRA of 3569 genes")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*model.lra.ccc[,dim], pch = 19, col = "gray80", cex = 0.5)
text(model.lra.rpc[,dim], labels = types, col = type.cols[types], cex = 0.7, font = 2) 



   
legend("bottomleft", legend=ABnames, 
       pch=tellus.pch, col=tellus.col, text.col=tellus.col, pt.cex=0.6, cex=0.8)


### ALR
model.alr <- PCA(ALR(model0, denom=2138, weight=FALSE)$LR, weight=FALSE)
model.alr.rpc <- -model.alr$rowpcoord 
model.alr.ccc <- -model.alr$colcoord * sqrt(model.alr$colmass)
100*model.alr$sv^2/sum(model.alr$sv^2)
 [1] 2.043512e+01 1.805437e+01 9.333898e+00 7.286892e+00 5.888677e+00 5.067924e+00 3.162319e+00 2.780785e+00
 [9] 2.650742e+00 2.523144e+00 2.477763e+00 2.378629e+00 2.284329e+00 2.013036e+00 1.801611e+00 1.674334e+00
[17] 1.654793e+00 1.561244e+00 1.515620e+00 1.441875e+00 1.434793e+00 1.360053e+00 1.218051e+00 3.588396e-29

rescale <- 20
dim <- c(3,4)
col <- c("blue","red")
perc.hor <- 9.3; perc.ver <- 7.3

plot(1.05 * rbind(model.alr.rpc[,dim], rescale*model.alr.ccc[,dim]), type = "n", asp = 1, 
     xlab = paste("PCA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("PCA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "PCA of 3568 ALRs w.r.t. reference #2138")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*model.alr.ccc[,dim], pch = 19, col = "gray80", cex = 0.5)
text(model.alr.rpc[,dim], labels = types, col = type.cols[types], cex = 0.7, font = 2) 

   
### CA with 0.5 power transformation
genes05.ca <- CA(CLOSE(genes0^0.5))
genes05.rpc <- -genes05.ca$rowpcoord/0.5
genes05.ccc <- -genes05.ca$colcoord * sqrt(genes05.ca$colmass)
100*genes05.ca$sv^2/sum(genes05.ca$sv^2)
[1] 22.093012 14.848020  8.329700  8.191394  3.483925  3.137939  2.981124  2.877830  2.695056  2.642760  2.586955
[12]  2.528250  2.479204  2.422827  2.363432  2.286030  2.233785  2.189060  2.134450  2.032936  1.997756  1.805325
[23]  1.659232

rescale <- 15
dim <- c(3,4)
col <- c("blue","red")
perc.hor <- 8.3; perc.ver <- 8.2

plot(1.05 * rbind(genes05.rpc[,dim], rescale*genes05.ccc[,dim]), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "CA of square-root compositions (with zeros)")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
    col = "black", col.ticks = col[2], col.axis = col[2])
points(rescale*genes05.ccc[,dim], pch = 19, col = "gray80", cex = 0.5)
text(genes05.rpc[,dim], labels = types, col = type.cols[types], cex = 0.7, font = 2) 

genes05.cm <- colMeans(CLOSE(genes0^0.5))
hist(100*genes05.cm, col="lightblue", main="Gene weights in CA", xlab="CA column masses (%)")

### Ionas claims that gene #2138 is the closest to the geometric mean
model <- model0
model.gm <- exp(apply(log(model), 1, mean)) 
cor(model[,2138], model.gm)
[1] 0.2789223

cors.gm <- rep(0, ncol(model))
for(j in 1:ncol(model)) cors.gm[j] <- cor(model[,j], model.gm)
hist(cors.gm, main="Correlations of parts with gm", xlab="Correlations", col="lightblue")


var(log(model[,2]/model.gm))
[1] 1.221584

vars.gm <- rep(0, ncol(model))
for(j in 1:ncol(model)) vars.gm[j] <- var(log(model[,j]/model.gm))
hist(vars.gm, main="Variances of log(part/gm)", xlab="Variances", col="lightblue", breaks=seq(0,10,0.5))



###################################################################################################################
### correlation exercise with genes
apply(genes0==0, 1, sum) 
 1_1  1_2  1_3  1_4  1_5  1_6  2_1  2_2  2_3  2_4  2_5  2_6  3_1  3_2  3_3  3_4  3_5  3_6  4_1  4_2  4_3  5_1  5_2 
1078  791 1118 1095 1066  890  884  889  953 1028 1156 1099 1038 1043 1188 1026 1037  994 1149 1479  765  887  732 
 5_3 
 819 

genes0.order <- order(colSums(genes0), decreasing=TRUE)
par(mfrow=c(2,4), font.lab=2)
for(k in c(0.1,0.2,0.3,0.4,0.5,0.6,1,5)) {
  foo <- genes0[, genes0.order][,1:100*k]
  foo.cor <- as.dist(cor(foo))
  foo.norm <- CLOSE(foo)
  foo.norm.cor <- as.dist(cor(foo.norm))
  plot(foo.cor, foo.norm.cor, xlab="Correlations in counts", ylab="Correlations in compositions", 
       main=paste("Number of parts=",100*k, sep="") )
  abline(v=0, h=0, col="red", lwd=2)
  c00 <- sum(foo.cor<0 & foo.norm.cor<0)
  c11 <- sum(foo.cor>0 & foo.norm.cor>0)
  c01 <- sum(foo.cor<0 & foo.norm.cor>0)
  c10 <- sum(foo.cor>0 & foo.norm.cor<0)
  text(-0.25,-0.25, label=paste(round(100*c00/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
  text(-0.25,0.5, label=paste(round(100*c01/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
  text( 0.5,-0.25, label=paste(round(100*c10/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
  text( 0.5, 0.5, label=paste(round(100*c11/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
}


genes0.norm.order <- order(colSums(CLOSE(genes0)), decreasing=TRUE)
par(mfrow=c(2,4), font.lab=2)
for(k in c(1,5,10,15,20,25,30,35)) {
  foo     <- CLOSE(genes0[,genes0.norm.order])
  foo.cor <- as.dist(cor(foo)[1:100*k,1:100*k])
  foo.sub <- CLOSE(genes0[,genes0.norm.order][,1:100*k])
  foo.sub.cor <- as.dist(cor(foo.sub))
  plot(foo.cor, foo.sub.cor, xlab="Correlations in composition", ylab="Correlations in subcomposition", 
       main=paste("Number of parts=",100*k, sep="") )
  abline(v=0, h=0, col="red", lwd=2)
  c00 <- sum(foo.cor<0 & foo.sub.cor<0)
  c11 <- sum(foo.cor>0 & foo.sub.cor>0)
  c01 <- sum(foo.cor<0 & foo.sub.cor>0)
  c10 <- sum(foo.cor>0 & foo.sub.cor<0)
  text(-0.25,-0.25, label=paste(round(100*c00/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
  text(-0.25,0.5, label=paste(round(100*c01/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
  text( 0.5,-0.25, label=paste(round(100*c10/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
  text( 0.5, 0.5, label=paste(round(100*c11/length(foo.cor),1),"%", sep=""), col="red", font=2, cex=1.2)
}


#########################################################################################################
### incoherence in CA 
incoh.CA   <- matrix(0, nrow=100, ncol=44)
procr.CA.foo   <- matrix(0, nrow=100, ncol=44)
#procr.CA.05 <- matrix(0, nrow=10, ncol=44)
 
# genes.pro.cm.order <- order(genes.pro.cm, decreasing=TRUE)

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

### genes.pro.cm are the gene averages
D.chi      <- as.matrix(chidist(genes.pro, 2))
# D.chi.05   <- as.matrix(chidist(CLOSE(genes.pro^0.5), 2))

set.seed(1234567)
for(j in seq(24,4,-2)) {
  nparts <- round((j/52)*3569)
  for(i in 1:100) {
# find the subcompositional parts  
    jparts <- sample(1:3569, nparts)
    foo <- CLOSE(genes.pro[,jparts])
# remove parts all zeros
    allzero <- which(colSums(foo)==0)
    if(length(allzero)>0) {
      jparts <- jparts[-allzero]
      foo <- CLOSE(genes.pro[,jparts])
    }
 
    D <- as.dist(D.chi[jparts, jparts])
#    D.05 <- as.dist(D.chi.05[jparts, jparts])
    D.rpc <- cmdscale(D, eig=TRUE, k=24)$points
#    D.rpc.05 <- cmdscale(D.05, eig=TRUE, k=length(jparts)-1)$points

# from samples point of view
#    tellus.sca <- ca(tellus.pro, subsetcol=jparts)
#    tellus.sca.05 <- ca(tellus.pro^0.5, subsetcol=jparts)
#    Ds.rpc <- tellus.sca$rowcoord %*% diag(tellus.sca$sv)
#    Ds.rpc.05 <- tellus.sca.05$rowcoord %*% diag(tellus.sca.05$sv)
#    Ds <- dist(Ds.rpc)
#    Ds.05 <- dist(Ds.rpc.05)
# remove samples all zeros
    allzero <- which(rowSums(foo)==0)
    if(length(allzero)>0) {
      foo <- foo[-allzero,]
    }
    D2 <- chidist(foo, 2)
#    D2.05 <- chidist(foo^0.5, 2)
    D2.rpc <- cmdscale(D2, eig=TRUE, k=24)$points
#    D2.rpc.05 <- cmdscale(D2.05, eig=TRUE, k=length(jparts)-1)$points 
#    D2s <- chidist(foo)
#    alpha <- sum(D2s*Ds)/sum(Ds*2)
    procr.CA.foo[i,j] <- protest(D2.rpc, D.rpc, permutations=0)$t0
#    procr.CA.05[i,j] <- protest(D2.rpc.05, D.rpc.05, permutations=0)$t0
#    incoh.CA[i,j] <- sum((Ds-alpha*alpha*D2s)^2)/sum(Ds^2)
  }
}

procr.CA.quants <- apply(procr.CA, 2, quantile, c(0.025,0.975), na.rm=TRUE)
round(procr.CA.quants[,seq(4,44,2)],4)

procr.CA.ones <- rep(0,44)
for(j in seq(4,44,2)) procr.CA.ones[j] <- sum(procr.CA[,j]>0.999)
procr.CA.ones[seq(4,44,2)]

procr.CA.med <- rep(0,44)
for(j in 1:24) procr.CA.med[j] <- median(procr.CA.foo[,j])
for(j in 25:44) procr.CA.med[j] <- median(procr.CA[,j])

### Figure 17
# pdf(file="incoherence3.pdf", width=7.5, height=4, useDingbats=FALSE, family="ArialMT")
par(mar=c(5,5,1,1), mgp=c(3.5,0.7,0), font.lab=2, las=1, mfrow=c(1,1))
plot(rep(seq(4,44,2), each=2), as.numeric(procr.CA.quants[,seq(4,44,2)]), xlab="Number of parts in subcomposition",
     ylab="Procrustes correlation", bty="n", xaxt="n", ylim=c(0.95, 1.005), type="n", font.lab=2, xlim=c(4,45))
axis(1, at=seq(4,44,2), labels=round((seq(4,44,2)/52)*3569))
for(j in seq(4,44,2)) segments(j, procr.CA.quants[1,j], j, procr.CA.quants[2,j], col="blue", lwd=2)
eps <- 0.2
for(j in seq(4,44,2)) segments(j-eps, procr.CA.quants[1,j], j+eps, procr.CA.quants[1,j], col="blue", lwd=2, lend=2)
for(j in seq(4,44,2)) segments(j-eps, procr.CA.quants[2,j], j+eps, procr.CA.quants[2,j], col="blue", lwd=2, lend=2)

points(seq(4,44,2), procr.CA.med[seq(4,44,2)], pch=21, col="blue", bg="white", cex=0.7)
text(seq(4,44,2), rep(1.0025, 21), labels=procr.CA.ones[seq(4,44,2)], font=2, cex=0.8)

dev.off()


procr.CA.comp <- matrix(0,100,11)
set.seed(12345)
for(j in 1:11) {
  nparts <- round(((j*10-8)/100)*3569)
  if(j==11) nparts <- 3569
  for(k in 1:100) {
    foo <- CLOSE(genes0.pro[,sample(1:3569, nparts)])
### 50% sample in subcomposition
    subsample <- sample(1:nparts, round(nparts/5))
    foo.sub <- CLOSE(foo[,subsample])
    foo.cpc <- CA(foo)$colpcoord[subsample,]
    foo.sub.cpc <- CA(foo.sub)$colpcoord
    procr.CA.comp[k,j] <- protest(foo.cpc, foo.sub.cpc, permutations=0)$t0  
  } 
}
procr.CA.comp.quant <- apply(procr.CA.comp, 2, quantile, c(0.025,0.5,0.975))  
round(procr.CA.comp.quant, 3)
nparts <- c(round((((1:10)*10-8)/100)*3569), 3569)
procr.CA.comp.ones <- rep(0,11)
for(j in 1:11) procr.CA.comp.ones[j] <- sum(procr.CA.comp[,j]>0.999)
procr.CA.comp.ones

### Figure 17
# pdf(file="incoherence4.pdf", width=7.5, height=4, useDingbats=FALSE, family="ArialMT")
par(mar=c(5,5,1,1), mgp=c(2.5,0.7,0), font.lab=2, las=1, mfrow=c(1,1), cex.axis=0.9)
plot(rep(1:11, each=3), as.numeric(procr.CA.comp.quant), xlab="Number of parts in composition",
     ylab="Procrustes correlation (coherence)", bty="n", xaxt="n", ylim=c(0.95, 1.005), type="n", font.lab=2,
     main="20% subcompositions of compositions of increasing sizes")
axis(1, at=1:11, labels=nparts)
for(j in 1:11) segments(j, procr.CA.comp.quant[1,j], j, procr.CA.comp.quant[3,j], col="blue", lwd=2)
eps <- 0.05
for(j in 1:11) segments(j-eps, procr.CA.comp.quant[1,j], j+eps, procr.CA.comp.quant[1,j], col="blue", lwd=2, lend=2)
for(j in 1:11) segments(j-eps, procr.CA.comp.quant[3,j], j+eps, procr.CA.comp.quant[3,j], col="blue", lwd=2, lend=2)

points(1:11, procr.CA.comp.quant[2,], pch=21, col="blue", bg="white", cex=0.8)
text(1:11, rep(1.0025, 11), labels=procr.CA.comp.ones, font=2, cex=0.8)

dev.off()

procr.CA.comp05 <- matrix(0,100,11)
set.seed(12345)
for(j in 1:11) {
  nparts <- round(((j*10-8)/100)*3569)
  if(j==11) nparts <- 3569
  for(k in 1:100) {
    foo <- genes0.pro[,sample(1:3569, nparts)]
### 50% sample in subcomposition
    subsample <- sample(1:nparts, round(nparts/5))
    foo.sub <- foo[,subsample]
    foo <- CLOSE(foo^0.5)
    foo.sub <- CLOSE(foo.sub^0.5)
    foo.cpc <- CA(foo)$colpcoord[subsample,]
    foo.sub.cpc <- CA(foo.sub)$colpcoord
    procr.CA.comp05[k,j] <- protest(foo.cpc, foo.sub.cpc, permutations=0)$t0  
  } 
}
procr.CA.comp05.quant <- apply(procr.CA.comp05, 2, quantile, c(0.025,0.5,0.975))  
round(procr.CA.comp05.quant, 3)
nparts <- c(round((((1:10)*10-8)/100)*3569), 3569)
procr.CA.comp05.ones <- rep(0,11)
for(j in 1:11) procr.CA.comp05.ones[j] <- sum(procr.CA.comp05[,j]>0.999)
procr.CA.comp05.ones

### Figure 17
# pdf(file="incoherence5.pdf", width=7.5, height=4, useDingbats=FALSE, family="ArialMT")
par(mar=c(5,5,1,1), mgp=c(2.5,0.7,0), font.lab=2, las=1, mfrow=c(1,1), cex.axis=0.9)
plot(rep(1:11, each=3), as.numeric(procr.CA.comp05.quant), xlab="Number of parts in composition",
     ylab="Procrustes correlation (coherence)", bty="n", xaxt="n", ylim=c(0.95, 1.005), type="n", font.lab=2,
     main="20% subcompositions of compositions^0.5 of increasing sizes")
axis(1, at=1:11, labels=nparts)
for(j in 1:11) segments(j, procr.CA.comp05.quant[1,j], j, procr.CA.comp05.quant[3,j], col="blue", lwd=2)
eps <- 0.05
for(j in 1:11) segments(j-eps, procr.CA.comp05.quant[1,j], j+eps, procr.CA.comp05.quant[1,j], col="blue", lwd=2, lend=2)
for(j in 1:11) segments(j-eps, procr.CA.comp05.quant[3,j], j+eps, procr.CA.comp05.quant[3,j], col="blue", lwd=2, lend=2)

points(1:11, procr.CA.comp05.quant[2,], pch=21, col="blue", bg="white", cex=0.8)
text(1:11, rep(1.0025, 11), labels=procr.CA.comp05.ones, font=2, cex=0.8)

dev.off()
  
