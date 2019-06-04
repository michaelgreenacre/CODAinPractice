invSLR <- function (SLRmatrix, parts=NA, ratios=colnames(SLRmatrix))
# SLRmatrix contains the I x (J-1) matrix of SLRs
# parts contains the names of the parts (required)
# ratios contains the definition of the ratios, using the part names
# SLRdef is the intermediate J-1 x J matric of 1s and -1s indicating the 
#   numerator and denominator parts
{
    J <- ncol(SLRmatrix)+1                  # number of parts
    inv <- matrix(0, nrow(SLRmatrix), J)    # resultant inverses
    for(i in 1:nrow(SLRmatrix)) {
      A <- matrix(0, J, J)
      for(j in 1:(J-1)) {
# numerator and denominator of the ratio
        split <- strsplit(ratios[j], "/")[[1]]
# numerator parts
        numer <- strsplit(split[1], "&")[[1]]
# denominator parts
        denom <- strsplit(split[2], "&")[[1]]
        for(k in 1:length(numer)) A[j, which(parts==numer[k])]  <- 1
        for(k in 1:length(denom)) A[j, which(parts==denom[k])] <- -exp(SLRmatrix[i,j])
      }
      A[J,] <- rep(1, J)
      b <- rep(0, J)
      b[J] <- 1
      inv[i,] <- solve(A, b)
    }
    colnames(inv) <- parts
    inv
}


### -------------- TEST --------------
# 3-part example of Martin-Fernandez
### SIMULATION EXAMPLE
# means ilr1, ilr2 from M.Greenacre example
set.seed(1)
rilr1=rnorm(50,0.5653601,0.2)
rilr2=rnorm(50,-1.010777,0.06)
############ ilr data set
rilr=cbind(rilr1,rilr2)
cov(rilr)
cor(rilr) # independent ILRs
########### contrast matrix "wine/beer" and "spirits/(wine&beer)"
FI=rbind(c(sqrt(0.5),-sqrt(0.5),0),c(-sqrt(1/6),-sqrt(1/6),sqrt(2/3)))
########## clr coordinates
rclr=rilr%*%FI
######### raw data
rX=exp(rclr)/apply(exp(rclr),1,sum)
colnames(rX)<-c("wine", "beer" , "spirits")
head(rX)
          wine      beer    spirits
[1,] 0.5696636 0.3057274 0.12460898
[2,] 0.6219330 0.2654298 0.11263728
[3,] 0.5574888 0.3174289 0.12508223
[4,] 0.6997466 0.2003298 0.09992368
[5,] 0.6189221 0.2534688 0.12760907
[6,] 0.5496008 0.3115986 0.13880069

SvsBW <- SLR(rX, 3, 1:2)$LR
BvsW  <- SLR(rX, 2, 1)$LR

rX.SLR <- cbind(SvsBW, BvsW)
colnames(rX.SLR) <- c("spirits/beer&wine", "beer/wine")
rX.def <- matrix(0, 2, 3)
rX.def[1,3] <- 1;  rX.def[1,1] <- rX.def[1,2] <- -1
rX.def[2,2] <- 1;  rX.def[2,1] <- -1
rX.def
     [,1] [,2] [,3]
[1,]   -1   -1    1
[2,]   -1    1    0

invSLR(rX.SLR, rX.def)[1:6,]
          [,1]      [,2]       [,3]
[1,] 0.5696636 0.3057274 0.12460898
[2,] 0.6219330 0.2654298 0.11263728
[3,] 0.5574888 0.3174289 0.12508223
[4,] 0.6997466 0.2003298 0.09992368
[5,] 0.6189221 0.2534688 0.12760907
[6,] 0.5496008 0.3115986 0.13880069

head(rX)
          [,1]      [,2]       [,3]
[1,] 0.5696636 0.3057274 0.12460898
[2,] 0.6219330 0.2654298 0.11263728
[3,] 0.5574888 0.3174289 0.12508223
[4,] 0.6997466 0.2003298 0.09992368
[5,] 0.6189221 0.2534688 0.12760907
[6,] 0.5496008 0.3115986 0.13880069

# simple logratios
rX.LR <- LR(rX)$LR
head(rX.LR)
     wine/beer wine/spirits beer/spirits
[1,] 0.6223520     1.519865    0.8975132
[2,] 0.8514821     1.708660    0.8571775
[3,] 0.5631885     1.494471    0.9312826
[4,] 1.2507535     1.946312    0.6955581
[5,] 0.8927388     1.579008    0.6862692
[6,] 0.5674764     1.376153    0.8086767

# just need two of them, let's take wine/spirits and beer/spirits
rX.LR2 <- rX.LR[,2:3]

ratios <- colnames(rX.SLR)   # "spirits/beer&wine" "beer/wine"
parts <- colnames(rX)        # "wine"    "beer"    "spirits"

SLRinverses <- invSLR(rX.SLR, parts=parts, ratios=ratios)
head(SLRinverses)


### on Aar Massif
mafic      <- apply(aar[,c(5,10,4)], 1, sum)
felsic     <- apply(aar[,c(7,1,3,8)], 1, sum)
carbonate  <- apply(aar[,c(6,9)], 1, sum)
aar.amalg  <- cbind(aar, mafic, felsic, carbonate)

#   Perform the stepwise analysis
aar.step   <- STEP(aar.amalg, aar, weight=FALSE)

aar.step$names
[1] "MgO/Na2O"         "K2O/P2O5"         "SiO2/K2O"         "TiO2/Na2O"       
[5] "SiO2/Na2O"        "felsic/carbonate" "MnO/carbonate"    "Al2O3/MgO"       
[9] "TiO2/Fe2O3t"

ratios <- aar.step$names
ratios[6] <- "Na2O&SiO2&Al2O3&K2O/CaO&P2O5"
ratios[7] <- "MnO/CaO&P2O5"
parts <- colnames(aar)

SLR.inverses <- invSLR(aar.step$logratios, parts, ratios)
round(head(SLR.inverses), 5)

sum((aar-SLRinverses)^2)
[1] 3.647801e-30


aar.step2 <- STEP(aar, weight=FALSE)
ratios2 <- aar.step2$names
parts2  <- colnames(aar)
SLR.inverses2 <- invSLR(aar.step2$logratios, parts2, ratios2)
round(head(SLR.inverses2), 5)
       SiO2    TiO2   Al2O3     MnO     MgO     CaO    Na2O     K2O    P2O5  Fe2O3t
[1,] 0.73638 0.00221 0.14236 0.00037 0.00482 0.01043 0.04344 0.04194 0.00060 0.01746
[2,] 0.73958 0.00260 0.14019 0.00038 0.00511 0.01142 0.04326 0.03955 0.00070 0.01722
[3,] 0.74027 0.00301 0.13883 0.00041 0.00551 0.01153 0.04320 0.03829 0.00070 0.01824
[4,] 0.78219 0.00320 0.11374 0.00043 0.00571 0.01021 0.03414 0.03164 0.00060 0.01812
[5,] 0.78943 0.00391 0.10313 0.00061 0.00862 0.01072 0.02816 0.02926 0.00070 0.02546
[6,] 0.74817 0.00543 0.11942 0.00073 0.00945 0.02232 0.02965 0.03217 0.00292 0.02975


sum((aar-SLR.inverses2)^2)
[1] 2.484784e-30

### is the set of logratio amalgamations a basis?

logratios <- as.data.frame(aar.step$logratios)
for(j in 1:ncol(aar)) {
  part.lm <- lm(aar.CLR[,j] ~ ., data=logratios)
  print(summary(part.lm)$r.squared)
}
[1] 0.9999948
[1] 0.9999855
[1] 0.9999837
[1] 0.9998461
[1] 0.9999962
[1] 0.9995555
[1] 0.9999952
[1] 0.999989
[1] 0.9999936
[1] 0.9999895

logratios.ilr <- logratios
# felsic     <- apply(aar[,c(7,1,3,8)], 1, sum)
# carbonate  <- apply(aar[,c(6,9)], 1, sum)
# ratios[6] <- "Na2O&SiO2&Al2O3&K2O/CaO&P2O5"
# ratios[7] <- "MnO/CaO&P2O5"
logratios.ilr[,6] <- ILR(aar, c(7,1,3,8), c(6,9), weight=FALSE)$LR 
logratios.ilr[,7] <- ILR(aar, 4, c(6,9), weight=FALSE)$LR 

for(j in 1:ncol(aar)) {
  part.lm <- lm(aar.CLR[,j] ~ ., data=logratios.ilr)
  print(summary(part.lm)$r.squared)
}
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1


logratios.lr <- as.data.frame(aar.step2$logratios)
for(j in 1:ncol(aar)) {
  part.lm <- lm(aar.CLR[,j] ~ ., data=logratios.lr)
  print(summary(part.lm)$r.squared)
}

[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1
[1] 1

