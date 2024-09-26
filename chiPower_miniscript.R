### Script to demonstrate chiPower transformation on Crohn data

### packages 'easyCODA' and 'vegan'required, install if necessary
# install.packages("easyCODA")
# install.packages("vegan")
library(easyCODA)
library(vegan)

### The Crohn data set is supplied on GitHub
load("Crohn_modified.RData")
ls()
# [1] "Crohn.x" "Crohn.y"
### x_Crohn is 975 x 48 table of counts
dim(x_Crohn)
# [1] 975  48
### y_Crohn is vector of binary responses (0 = no Crohn, 1 = Crohn)
table(y_Crohn)
# y_Crohn
#  CD  no 
# 662 313
### convert to factor, with CD as the second 
y_Crohn <- factor(y_Crohn, levels=c("no","CD")) 
table(y_Crohn)
# y_Crohn
#  no  CD 
# 313 662 
### x_Crohn is the modified Crohn dataset, which has 1 added to it
head(x_Crohn[,1:5])
                g__Turicibacter g__Parabacteroides g___Ruminococcus_ g__Sutterella f__Peptostreptococcaceae_g__
1939.SKBTI.0175               1               1191              2258         31712                            4
1939.SKBTI.1068               1                 17              1913            23                           32
1939.SKBTI047                 3                152               868           142                            1
1939.SKBTI051                 1                 48               370           846                          326
1939.SKBTI063                 1                 18              2830            38                            2
1939.SKBTI072                 1               6608               279          1920                            1
### (so, in above the 1s are originally 0s, and other counts incremented by 1)
### clounting the 1s gives the number of 0s in original data set
sum(x_Crohn == 1)
# [1] 13474
sum(x_Crohn == 1)/(975*48)
# [1] 0.287906
### (hence, more than 28% zeros in original Crohn dataset)

### The modified Crohn data set in in package coda4microbiome 
### To get the original dataset subtract 1 from all the elements
### Let's call that x_Crohn0
x_Crohn0 <- x_Crohn - 1

### --------------------------------------------------------------------------------
### chiPower function
chiPower<- function(data, close=TRUE, power=1, chi=TRUE, BoxCox=TRUE, center=TRUE) {
# function to compute chiPower transformation with various options
# close: close data
# power: tranformation (i.e., value of lambda)
# chi: chi-square standardization
# BoxCox: apply 1/power rescaling (without the subtraction of 1)
# center: center the final result by column means
   foo <- data^power
   if(close)   foo <- foo / rowSums(foo)
   if(chi)     foo <- sweep(foo, 2, sqrt(colMeans(foo)), FUN="/")
   if(BoxCox)  foo <- (1/power) * foo
   if(center)  foo <- sweep(foo, 2, colMeans(foo))
   if(!center) foo <- foo - 1/sqrt(ncol(data))
   data.chiPower <- foo
   data.chiPower
}
### ----------------------------------------------------------------------------------
### Various chiPower distances now compared to logratio distances
### Logratio distances have to be computed on x_Crohn
### whereas chiPower distances computed on original data with zeros (x_Crohn0)

### logratio distances (with column weight constant in this case, i.e. weight=FALSE)
### Use CLRs to get logratio distances
Crohn.clr <- CLR(CLOSE(x_Crohn), weight=FALSE)$LR
### LRA is PCA of CLRs (or you can get it directly using function LRA)
Crohn.lra <- PCA(Crohn.clr, weight=FALSE)
### Principal row coordinates define logratio geometry (last dimension has no variance)
Crohn.lra.rpc <- Crohn.lra$rowpcoord[,-(ncol(x_Crohn)-1)]
### Logratio distances from CLRs
Crohn.clr.dist <- dist(Crohn.clr)/sqrt(ncol(x_Crohn))

### loop on powers from 0.001 to 1 and find optimal using Procrustes
### comparison with chiPower transforms in two different ways
### 1. Spearman rank correlation between distances
### 2. Procrustes correlation between geometries
npower <- 100
Crohn.cor <- rep(0,npower)
Crohn.pro <- rep(0,npower)  
### Original composition on data set with 0s
Crohn0.comp <- CLOSE(x_Crohn0)    
### start of loop (takes a while) -----------
for(power in 1:npower) {
  foo.chiPower <- chiPower(Crohn0.comp, power=power/npower)
  foo.dist     <- dist(foo.chiPower)
  foo.rpc      <- PCA(foo.chiPower, weight=FALSE)$rowpcoord[,1:(ncol(Crohn.comp)-1)]
  Crohn.cor[power] <- cor(Crohn.clr.dist, foo.dist, method="spearman")
  Crohn.pro[power] <- protest(Crohn.lra.rpc, foo.rpc, permutations=0)$t0
}
### end of loop -----------------------------
max(Crohn.pro)
# [1] 0.9006722
which(Crohn.pro==max(Crohn.pro))
# [1] 25
max(Crohn.cor)
# [1] 0.7533164
which(Crohn.cor==max(Crohn.cor))
# [1] 30

### plot the Procrustes (black) and Spearman correlations (blue)
plot(Crohn.pro, type="l", lwd=2, xaxt="n", bty="n", ylim=c(0.3,0.9),
     xlab="Power of chiPower transformation", ylab="Procrustes correlation / Spearman")
axis(1, at=seq(0,100,10), labels=seq(0, 1, 0.1))
### blue curve for Spearman correlation
lines(Crohn.cor, col="blue")
### (Procrustes correlation maximum of 0.901 at power=0.25
###  Spearman correlation maximum of 0.753 at power=0.30)




