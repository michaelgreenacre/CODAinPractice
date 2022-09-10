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

  
