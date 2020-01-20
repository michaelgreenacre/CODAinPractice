### Aar Massif data set is in data object 'Aar'
### in the package 'compositions'
library(compositions)
data(Aar)
### Columns 3-12 are the oxide data
aar <- Aar[,3:12]
### close the data set
aar <- aar / rowSums(aar)
round(head(aar),5)
#      SiO2    TiO2   Al2O3     MnO     MgO     CaO    Na2O     K2O    P2O5  Fe2O3t
# 1 0.73638 0.00221 0.14236 0.00037 0.00482 0.01043 0.04344 0.04194 0.00060 0.01746
# 2 0.73958 0.00260 0.14019 0.00038 0.00511 0.01142 0.04326 0.03955 0.00070 0.01722
# 3 0.74027 0.00301 0.13883 0.00041 0.00551 0.01153 0.04320 0.03829 0.00070 0.01824
# 4 0.78219 0.00320 0.11374 0.00043 0.00571 0.01021 0.03414 0.03164 0.00060 0.01812
# 5 0.78943 0.00391 0.10313 0.00061 0.00862 0.01072 0.02816 0.02926 0.00070 0.02546
# 6 0.74817 0.00543 0.11942 0.00073 0.00945 0.02232 0.02965 0.03217 0.00292 0.02975

### compute the (unweighted) logratio variance three different ways
require(easyCODA)
### (1) as the sum of all pairwise logratio variances divided by 100 (=J squared)
### option weight=FALSE implies equal weights, each logratio gets weight (1/10)x(1/10)
### notice that the weights are always passed on in the logratio object
aar.LR <- LR(aar, weight=FALSE)
LR.VAR(aar.LR)
# [1] 0.1522912
### (2) as the sum of all CLR variances divided by 10 (=J)
aar.CLR <- CLR(aar, weight=FALSE)
LR.VAR(aar.CLR)
# [1] 0.1522912
### (3) as the average of all squared elements of double-centred log-transformed data matrix
aarlog <- log(aar)
aarlog.DC <- sweep(aarlog, 1, apply(aarlog, 1, mean))
aarlog.DC <- sweep(aarlog.DC, 2, apply(aarlog.DC, 2, mean))
sum(aarlog.DC^2) /(nrow(aar)*ncol(aar))
# [1] 0.1522912

### Compute the three amalgamations
mafic      <- rowSums(aar[,c(5,10,4)])
felsic     <- rowSums(aar[,c(7,1,3,8)])
carb_apat  <- rowSums(aar[,c(6,9)])
### Add the three amalgamations to the existing parts
aar.amalg  <- cbind(aar, mafic, felsic, carb_apat)

### Perform the stepwise logratio selection using function
### STEP in the easyCODA package: logratios are formed from the
### columns in aar.amalg, with the target aar to be explained
### (see Option 4 in Section 2)
aar.step   <- STEP(aar.amalg, aar, weight=FALSE)
### Names of the logratios in $names of the STEP object
aar.step$names
# [1] "MgO/Na2O"  "K2O/P2O5"          "SiO2/K2O"        "TiO2/Na2O"       
# [5] "SiO2/Na2O" "felsic/carb_apat"  "MnO/carb_apat"   "Al2O3/MgO"       
# [9] "TiO2/Fe2O3t"
ratio.names <- aar.step$names
### The ones involving amalgamations have to be renamed
### with the definitions of the amalgamations, as required by invSLR
ratio.names[6] <- "Na2O&SiO2&Al2O3&K2O/CaO&P2O5"
ratio.names[7] <- "MnO/CaO&P2O5"
### The names of the parts (for matching with logratio names)
part.names <- colnames(aar)

### Transform the logratios back to the original parts using invSLR
### The logratios are in $logratios of the STEP object 
SLRinverses <- invSLR(aar.step$logratios, part.names, ratio.names)
round(head(SLRinverses), 5)
 #         SiO2    TiO2   Al2O3     MnO     MgO     CaO    Na2O     K2O    P2O5  Fe2O3t
 # [1,] 0.73638 0.00221 0.14236 0.00037 0.00482 0.01043 0.04344 0.04194 0.00060 0.01746
 # [2,] 0.73958 0.00260 0.14019 0.00038 0.00511 0.01142 0.04326 0.03955 0.00070 0.01722
 # [3,] 0.74027 0.00301 0.13883 0.00041 0.00551 0.01153 0.04320 0.03829 0.00070 0.01824
 # [4,] 0.78219 0.00320 0.11374 0.00043 0.00571 0.01021 0.03414 0.03164 0.00060 0.01812
 # [5,] 0.78943 0.00391 0.10313 0.00061 0.00862 0.01072 0.02816 0.02926 0.00070 0.02546
 # [6,] 0.74817 0.00543 0.11942 0.00073 0.00945 0.02232 0.02965 0.03217 0.00292 0.02975
### The values check with the part values shown at the start

### The stepwise search for the best pairwise logratios
### (see Option 3 in Section 2)
require(easyCODA)
aar.step2 <- STEP(aar, weight=FALSE)
ratio.names2 <- aar.step2$names
part.names2  <- colnames(aar)

### The inverse transformation
SLR.inverses2 <- invSLR(aar.step2$logratios, part.names2, ratio.names2)
round(head(SLR.inverses2), 5)
#        SiO2    TiO2   Al2O3     MnO     MgO     CaO    Na2O     K2O    P2O5  Fe2O3t
# [1,] 0.73638 0.00221 0.14236 0.00037 0.00482 0.01043 0.04344 0.04194 0.00060 0.01746
# [2,] 0.73958 0.00260 0.14019 0.00038 0.00511 0.01142 0.04326 0.03955 0.00070 0.01722
# [3,] 0.74027 0.00301 0.13883 0.00041 0.00551 0.01153 0.04320 0.03829 0.00070 0.01824
# [4,] 0.78219 0.00320 0.11374 0.00043 0.00571 0.01021 0.03414 0.03164 0.00060 0.01812
# [5,] 0.78943 0.00391 0.10313 0.00061 0.00862 0.01072 0.02816 0.02926 0.00070 0.02546
# [6,] 0.74817 0.00543 0.11942 0.00073 0.00945 0.02232 0.02965 0.03217 0.00292 0.02975
