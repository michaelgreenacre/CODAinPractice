### read data from github site
cancer <- read.table("https://raw.githubusercontent.com/michaelgreenacre/CODAinPractice/master/Baxter_OTU_table.txt", 
                     header=TRUE, check.names=FALSE)
dim(cancer)
# [1] 490 338

### remove first three columns to get the OTU dataset
cancer <- cancer[,-c(1:3)]
cancer.no0 <- cancer+1
# remove strong outlier, possibly an error
cancer.pro <- cancer.no0[,-265]/rowSums(cancer.no0[,-265])

### load easyCODA package
require(easyCODA)

### unweighted (i.e. equally weighted) option
starttime <- Sys.time()
cancer.alr <- FINDALR(cancer.pro, weight=FALSE)
endtime <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 182.8436 secs on Toshiba latop

plot(cancer.alr$var.log,cancer.alr$procrust.cor)
cancer.alr$tot.var       # [1] 1.530197
cancer.alr$procrust.max  # [1] 0.9355615
cancer.alr$procrust.ref  # [1] 269
cancer.alr$var.min       # [1] 0.3082664
cancer.alr$var.ref       # [1] 320

### weighted option
starttime <- Sys.time()
cancer.alrw <- FINDALR(cancer.pro, weight=TRUE)
endtime <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 190.9042 secs on Toshiba laptop

cancer.alrw$tot.var      # [1] 2.709184
cancer.alrw$procrust.max # [1] 0.9525636
cancer.alrw$procrust.ref # [1] 269
cancer.alrw$var.min      # [1] 0.3082664
cancer.alrw$var.ref      # [1] 320


### -----------------------------------------------------------------------------------------
### read meta data of Baxter microbiome study
meta <- read.delim("https://raw.githubusercontent.com/michaelgreenacre/CODAinPractice/master/Baxter_Metadata.txt", 
                   header=TRUE, check.names=FALSE, sep="\t")
dim(meta)
# [1] 490  27
colnames(meta)
#  [1] "sample"       "fit_result"   "Site"         "Dx_Bin"       "dx"           "Hx_Prev"      "Hx_of_Polyps"
#  [8] "Age"          "Gender"       "Smoke"        "Diabetic"     "Hx_Fam_CRC"   "Height"       "Weight"      
# [15] "BMI"          "White"        "Native"       "Black"        "Pacific"      "Asian"        "Other"       
# [22] "Ethnic"       "NSAID"        "Abx"          "Diabetes_Med" "stage"        "Location" 

# the group labels, also convert to numbers
dx <- meta[,"dx"]
table(dx)
# adenoma  cancer  normal 
#     198     120     172 
dx.num <- as.numeric(as.factor(dx))
table(dx.num)
#   1   2   3 
# 198 120 172 

### perform RDA of CLRs of the OTUS on the categorical variable dx (groups)
cancer.rda <- rda(CLR(cancer.pro, weight=FALSE)$LR ~ factor(dx))

# Inertia Proportion Rank
# Total         5.121e+02  1.000e+00     
# Constrained   4.194e+00  8.189e-03    2     <- 0.82% of variance due to group differences
# Unconstrained 5.079e+02  9.918e-01  333

### permutation test of significance (9999 permutations)
set.seed(123)
anova(cancer.rda, permutations=9999)
# Df Variance      F Pr(>F)    
# Model      2     4.19 2.0106  1e-04 ***     <- nevertheless, highly significant, p<0.0001

### variances explained in constrained and full spaces
100*cancer.rda$CCA$eig/cancer.rda$CCA$tot.chi
#     RDA1     RDA2 
# 74.62057 25.37943 
100*cancer.rda$CCA$eig/cancer.rda$tot.chi
#      RDA1      RDA2 
# 0.6110919 0.2078403 


### row coordinates in exact geometry in restricted 2-d space of group differences
cancer.rda.wa <- cancer.rda$CCA$wa
cancer.procrust.cor <- rep(0,nrow=ncol(cancer.pro)) 

### loop through the reference components but fit in the reduced space
starttime <- Sys.time()
for(j in 1:ncol(cancer.pro)) {
  foo.alr <- ALR(cancer.pro, denom=j, weight=FALSE)$LR
  foo.rda <- rda(foo.alr ~ factor(dx))
  cancer.procrust.cor[j] <- protest(foo.rda$CCA$wa,cancer.rda.wa, permutations=0)$t0  
}
endtime <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 149.3339 secs on Toshiba laptop

max(cancer.procrust.cor)
# [1] 0.9996624

which(cancer.procrust.cor==max(cancer.procrust.cor))
# [1] 312
colnames(cancer.pro)[312]
# [1] "Otu000363"

### compute ALRs with this reference
cancer.alr312     <- ALR(cancer.pro, denom=312, weight=FALSE)$LR
cancer.alr312.rda <- rda(cancer.alr312 ~ factor(dx))
cancer.alr312.wa  <- cancer.alr312.rda$CCA$wa

### variances explained in constrained and full spaces
100*cancer.alr312.rda$CCA$eig/cancer.alr312.rda$CCA$tot.chi
#     RDA1     RDA2 
# 74.61953 25.38047 
100*cancer.alr312.rda$CCA$eig/cancer.alr312.rda$tot.chi
#      RDA1      RDA2 
# 0.4804437 0.1634142  

### plot 2-D configuration using all logratios
cancer.cols <- c("blue","red","forestgreen")

par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2)
plot(cancer.rda.wa, type="n", asp=1, main="Constrained LRA of OTUs", 
     xlab="LRA dimension 1 (74.6% / 0.61%)", ylab="LRA dimension 2 (25.4% / 0.21%)")
abline(v=0, h=0, col="gray", lty=2)
text(cancer.rda.wa, labels=substr(dx,1,1), col=cancer.cols[as.numeric(factor(dx))], cex=0.6)
set.seed(123)
CIplot_biv(cancer.rda.wa[,1],cancer.rda.wa[,2],
           group=factor(dx), shade=TRUE, add=TRUE, groupcols=cancer.cols,
           groupnames=c("A","C","N"))
set.seed(123)
CIplot_biv(cancer.rda.wa[,1],cancer.rda.wa[,2],
           group=factor(dx), add=TRUE, groupcols=cancer.cols,
           shownames=FALSE)

### plot 2-D configurations using best set of ALRs
par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2)
plot(cancer.alr312.wa, type="n", asp=1, main="RDA of ALRs w.r.t. 312", 
     xlab="RDA dimension 1 (74.6% / 0.48%)", ylab="RDA dimension 2 (25.4% / 0.16%)")
abline(v=0, h=0, col="gray", lty=2)
text(cancer.alr312.wa, labels=substr(dx,1,1), col=cancer.cols[as.numeric(factor(dx))], cex=0.6)
set.seed(123)
CIplot_biv(cancer.alr312.wa[,1],cancer.alr312.wa[,2],
           group=factor(dx), shade=TRUE, add=TRUE, groupcols=cancer.cols,
           groupnames=c("A","C","N"))
set.seed(123)
CIplot_biv(cancer.alr312.wa[,1],cancer.alr312.wa[,2],
           group=factor(dx), add=TRUE, groupcols=cancer.cols,
           shownames=FALSE)


### ----------------------------------------------------------------------------
### do all of above again with weighted components
### perform RDA of CLRs of the OTUS on the categorical variable dx (groups)
cancer.rdaw <- rda(CLR(cancer.pro)$LR%*%diag(sqrt(c)) ~ factor(dx))

# Inertia Proportion Rank
# Total         2.714724   1.000000     
# Constrained   0.017192   0.006333    2      <- 0.63% explained by group differences
# Unconstrained 2.697533   0.993667  333

### permutation test of significance (9999 permutations)
set.seed(123)
anova(cancer.rdaw, permutations=9999)
Df Variance      F Pr(>F)  
Model      2  0.01719 1.5519 0.0156 *         <- still significant

### variances explained in constrained and full spaces
100*cancer.rdaw$CCA$eig/cancer.rdaw$CCA$tot.chi
#     RDA1     RDA2 
# 66.42767 33.57233 
100*cancer.rdaw$CCA$eig/cancer.rdaw$tot.chi
#      RDA1      RDA2 
# 0.4206694 0.2126049  


### row coordinates in exact geometry in restricted 2-d space of group differences
cancer.rdaw.wa <- cancer.rdaw$CCA$wa
cancer.procrustw.cor <- rep(0,nrow=ncol(cancer.pro)) 

### loop through the reference components but fit in the reduced space
starttime <- Sys.time()
for(j in 1:ncol(cancer.pro)) {
  cc <- c*c[j]
  cc <- cc[-j]
  foo.alr <- ALR(cancer.pro, denom=j, weight=FALSE)$LR
  foo.rda <- rda(foo.alr%*%diag(sqrt(cc)) ~ factor(dx))
  cancer.procrustw.cor[j] <- protest(foo.rda$CCA$wa,cancer.rdaw.wa, permutations=0)$t0  
}
endtime <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 161.5026 secs on Toshiba laptop

max(cancer.procrustw.cor)
# [1] 0.9982797

which(cancer.procrustw.cor==max(cancer.procrustw.cor))
# [1] 241
colnames(cancer.pro)[241]
# [1] "Otu000262"

### compute ALRs with this reference
cc <- c*c[241]
cc <- cc[-241]
cancer.alr241     <- ALR(cancer.pro, denom=241, weight=FALSE)$LR
cancer.alr241.rda <- rda(cancer.alr241%*%diag(sqrt(cc)) ~ factor(dx))
cancer.alr241.wa  <- cancer.alr241.rda$CCA$wa

### variances explained in constrained and full spaces
100*cancer.alr241.rda$CCA$eig/cancer.alr241.rda$CCA$tot.chi
#     RDA1     RDA2 
# 67.72902 32.27098  
100*cancer.alr241.rda$CCA$eig/cancer.alr241.rda$tot.chi
#      RDA1      RDA2 
# 0.3983381 0.1897969  

### plot 2-D configuration using all logratios
cancer.cols <- c("blue","red","forestgreen")

# invert second dimension to agree with previous plots
cancer.rdaw.wa[,2] <- -cancer.rdaw.wa[,2]

par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2)
plot(cancer.rdaw.wa, type="n", asp=1, main="Constrained weighted LRA of OTUs", 
     xlab="LRA dimension 1 (66.4% / 0.42%)", ylab="LRA dimension 2 (33.6% / 0.21%)")
abline(v=0, h=0, col="gray", lty=2)
text(cancer.rdaw.wa, labels=substr(dx,1,1), col=cancer.cols[as.numeric(factor(dx))], cex=0.6)
set.seed(123)
CIplot_biv(cancer.rdaw.wa[,1],cancer.rdaw.wa[,2],
           group=factor(dx), shade=TRUE, add=TRUE, groupcols=cancer.cols,
           groupnames=c("A","C","N"))
set.seed(123)
CIplot_biv(cancer.rdaw.wa[,1],cancer.rdaw.wa[,2],
           group=factor(dx), add=TRUE, groupcols=cancer.cols,
           shownames=FALSE)

### plot 2-D configurations using best set of ALRs
# invert second dimension to agree with previous plots
cancer.alr241.wa[,2] <- -cancer.alr241.wa[,2]

par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2)
plot(cancer.alr241.wa, type="n", asp=1, main="RDA of weighted ALRs w.r.t. 241", 
     xlab="RDA dimension 1 (67.7% / 0.40%)", ylab="RDA dimension 2 (32.3% / 0.19%)")
abline(v=0, h=0, col="gray", lty=2)
text(cancer.alr241.wa, labels=substr(dx,1,1), col=cancer.cols[as.numeric(factor(dx))], cex=0.6)
set.seed(123)
CIplot_biv(cancer.alr241.wa[,1],cancer.alr241.wa[,2],
           group=factor(dx), shade=TRUE, add=TRUE, groupcols=cancer.cols,
           groupnames=c("A","C","N"))
set.seed(123)
CIplot_biv(cancer.alr241.wa[,1],cancer.alr241.wa[,2],
           group=factor(dx), add=TRUE, groupcols=cancer.cols,
           shownames=FALSE)

### vaginal microbiome data set by Deng et al. (2018), cited and analysed by Wu et al. (2021)

vagina <- read.table("clipboard", check.names=FALSE)
vagina <- t(vagina)
dim(vagina)
# [1]  40 103
sum(vagina==0)/(nrow(vagina)*ncol(vagina))  #13% zeros
vagina1 <- vagina+1
vagina.pro <- vagina1/rowSums(vagina1)

starttime <- Sys.time()
vagina.alr <- FINDALR(vagina.pro)
endtime <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 0.7315788 secs on Toshiba.laptop

vagina.alr$tot.var
# [1] 3.257595
vagina.alr$procrust.max
# [1] 0.968543
vagina.alr$procrust.ref
# [1] 51

starttime <- Sys.time()
vagina.alrw <- FINDALR(vagina.pro, weight=TRUE)
endtime <- Sys.time()
difftime(endtime, starttime, units="secs")
# Time difference of 0.8090138 secs on Toshiba laptop

vagina.alrw$tot.var
# [1] 7.710076
vagina.alrw$procrust.max
# [1] 0.9825666
vagina.alr$procrust.ref
# [1] 12

### this illustrates getting ALR weights from function ALR
vagina.alr12 <- ALR(vagina.pro, denom=12)
vagina.alr12.LR <- vagina.alr12$LR
vagina.alr12.LRwt <- vagina.alr12$LR.wt

### exact weighted logratio geometry
vagina.lra      <- LRA(vagina.pro)
vagina.lra.rpc <- vagina.lra$rowpcoord
100*vagina.lra$sv[1:2]^2/sum(vagina.lra$sv^2)
# [1] 63.39400 24.44925

par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2)
plot(vagina.lra.rpc, type="n", asp=1, main="LRA of vaginal microbiome", 
     xlab="LRA dimension 1 (63.4%)", ylab="LRA dimension 2 (24.4%)")
abline(v=0, h=0, col="gray", lty=2)
text(vagina.lra.rpc, labels=rownames(vagina), col="blue", cex=0.6)


### plot 2-D configuration using best set of ALRs
### note that weights in the ALR object are used in the PCA
vagina.pca <- PCA(vagina.alr12)
vagina.pca.rpc <- vagina.pca$rowpcoord
vagina.lra.rpc <- vagina.lra$rowpcoord
100*vagina.pca$sv[1:2]^2/sum(vagina.pca$sv^2)
# [1] 76.67587 14.67200

par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2)
plot(vagina.pca.rpc, type="n", asp=1, main="PCA of ALRs w.r.t. 12", 
     xlab="PCA dimension 1 (74.7%)", ylab="RDA dimension 2 (14.7%)")
abline(v=0, h=0, col="gray", lty=2)
text(vagina.pca.rpc, labels=rownames(vagina), col="blue", cex=0.6)



