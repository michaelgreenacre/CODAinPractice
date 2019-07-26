# Small data set of ratios from Chapter 1
B6	1.576  6.333  4.019	  0.884  0.561  0.140
B7	1.224  4.800  3.920	  0.466  0.381  0.097
D4	2.184 16.462  7.538	 25.680 11.760  1.560
D5	3.127 35.650 11.400	 17.825  5.700  0.500
H5	0.504  0.418  0.830	  3.540  7.020  8.460
H6	0.946  0.409  0.432	  3.542  3.746  8.661

# Copy the above lines of data and read from clipboard
# Using Windows
ratios <- read.table("clipboard", row.names=1)
# Using MacOs
ratios <- read.table(pipe("pbpaste"), row.names=1
)
# colorspace package for HCL colour palette
require(colorspace)

# DOT function in easyCODA makes dot plot of ratios and logratios
require(easyCODA)
DOT(ratios, names=1:6, cols=rainbow_hcl(ncol(ratios)), 
    pch=rep(19, ncol(ratios)), ylim=c(0,40), ylab="Ratios")
DOT(log(ratios), names=1:6, cols=rainbow_hcl(ncol(ratios)), 
    pch=rep(19, ncol(ratios)), ylim=c(-3,4), ylab="Logratios")

# Data set 'veg' and function BAR in easyCODA
data(veg)
# Compositional bar plot
# Option order.column specifies which column orders the samples
# (if omitted, samples are plotted in orginal order)
BAR(veg, order.column=2, ylab="Vegetables", cols=rainbow_hcl(ncol(veg)))

# Function 'ternaryplot' in package vcd
require(vcd)
# Ternary plot of Fig. 2.5
par(mar=c(1,1,1,1))
ternaryplot(veg, cex=0.6, col="blue", id=rownames(veg), main="", 
            dimnames_color="red", labels="outside")

# Computing additive logratios
ALR(veg, denom=2)

# To extract the ALRs to an object
veg.ALR <- ALR(veg, denom=2)$LR

# Computing centred logratios
CLR(veg)

# Computing an isometric logratio
# numer = numerator set, denom = denominator set
ILR(veg, numer=1, denom=2:3)

# Computing an amalgamation logratio
# numer = numerator set, denom = denominator set
# e.g. this SLR of first part relative to sum of parts 2 and 3
# (for both ILR and SLR, the two subsets can be simply indicated in order)
SLR(veg, 1, 2:3)

# Computing a set of pivot logratios
# ordering = original order of parts by default
# but can be reordered as illustrated below
PLR(veg, ordering=c(1,3,2))

# Computing all pairwise logratios
# order as in columns of upper triangle above diagonal of square matrix
# i.e. 1st/2nd, 1st/3rd, 2nd/3rd, 1st/4th, 2nd/4th, 3rd/4th, 1st/5th etc...
LR(veg)

# Plot of two amalgamation logratios and corresponding isometric logratios
# (cf. Figs 3.2 (a) and (b) that plots the ratios on a log-scale
# whereas here the logratio values are plotted directly)
# First simplify the names of beans and potatos
rownames(veg)[c(2,9)] <- c("Beans", "Potatoes")
amalgs <- cbind(SLR(veg, 3, 1:2)$LR, SLR(veg, 1,2)$LR)
par(mar=c(4.2,4,3,1), mgp=c(2,0.7,0), font.lab=2)
plot(amalgs, type="n", xlab="log(Fat/[Protein+Carbohydrate])", 
     ylab="log(Protein/Carbohydrate)",  main="Amalgamations")
text(amalgs, labels=rownames(veg))
ilrs <- cbind(ILR(veg, 3, 1:2)$LR, ILR(veg, 1,2)$LR)
plot(ilrs, type="n", xlab="sqrt(2/3)*log(Fat/[Protein+Carbohydrate]^0.5)", 
     ylab="sqrt(1/2)*log(Protein/Carbohydrate)", main="ILRs")
text(ilrs, labels=rownames(veg))

# Inverting ALRs with respect to second part
veg.alr <- ALR(veg, denom=2)$LR
invALR(veg.alr, denom=2, part.names=colnames(veg))

# Computing total logratio variance and optionally individual variances
# e.g. total logratio variance of RomanCups data set, using 11 CLRs
# Data available in easyCODA package
data(cups)
LR.VAR(CLR(cups))

# Computing same total logratio variance, but using all 55 pairwise LRs
LR.VAR(LR(cups))

# Computing total variance and individual variances
LR.VAR(CLR(cups), vars=TRUE)

# Percentage contributions to total logratio 
cups.var <- LR.VAR(CLR(cups), vars=TRUE)
round(100 * cups.var$LRvars / cups.var$LRtotvar, 2)

# Percentages for unweighted data, showing the dominance of the rarest elements,
# particularly Mn which only has values (in percent) 0.01, 0.02 and 0.03
cups.var <- LR.VAR(CLR(cups, weight=FALSE), vars=TRUE)
round(100 * cups.var$LRvars / cups.var$LRtotvar, 2)

# Predicting sex of fish in data set FishMorphology
# fish data frame in easyCODA has sex, habitat and mass in first three columns
data(fish)
sex     <- fish[,1]
habitat <- fish[,2]
mass    <- fish[,3]
# remaining columns contain fish morphometric data, in object 'fishm' (75x26)
fishm <- as.matrix(fish[,4:29])
# convert fishm to a compositional data matrix
fishm <- fishm/apply(fishm, 1, sum)

# compute all logratios
fish.LR <- LR(fishm)$LR
# look for best logratio predictor of sex
# save deviance measures of lack of fit in 'dev'
# numLR = total number of logratios
numLR <- ncol(fishm) * (ncol(fishm)-1) / 2
dev   <- rep(0, numLR)
for(j in 1:numLR) dev[j] <- summary(glm(factor(sex) ~ fish.LR[,j], 
                                    family="binomial"))$aic
j.min <- which(dev == min(dev)) 
colnames(fish.LR)[j.min]
# result: [1] "Faw/Fdl"
fish.glm <- glm(factor(sex) ~ fish.LR[,j.min], family="binomial") 
summary(fish.glm)

# (no other logratios add significantly to this model)

# Predictions of log-odds: positive predictions (i.e. p>0.5) are for males
fish.pred <- predict(fish.glm)
table(fish.pred>0, factor(sex))

# Classification tree, predicting sex from morphometric logratios
require(rpart)
fish.tree <- rpart(factor(sex) ~ fish.LR)
plot(fish.tree, margin=0.1)
text(fish.tree, use.n=T)
fish.pred.tree <- predict(fish.tree)
table(fish.pred.tree[,1]<0.5, factor(sex))

# fishm contains morphometric measurements
# Scale the CLRs with square root of their weights 
# to ensure weighted logratio variance is explained
fish.CLR <- CLR(fishm)
fish.CLRw <- fish.CLR$LR %*% diag(sqrt(fish.CLR$LR.wt)) 

# logmass contains the logarithm of fish masses
logmass <- log(mass)
require(vegan)
fish.rda <- rda(fish.CLRw ~ logmass)
set.seed(123)
anova(fish.rda, by="terms")

# LRA, weighted and unweighted, of the RomanCups data set
# plotting by ca package function (default is symmetric plot)
# the mass option shows column symbols related to their weight
require(ca)
data(cups)
cups.wLRA <- LRA(cups)
cups.uLRA <- LRA(cups, weight = FALSE)
par(mar=c(3.5,3.3,4,2), font.lab=2, mgp=c(2,0.5,0), cex.axis=0.8)
plot(cups.uLRA, labels=c(1,2), main="Unweighted LRA")
plot(cups.wLRA, labels=c(1,2), mass=c(FALSE,TRUE), main="Weighted LRA") 

# The TimeBudget data set is in the R object 'time', available in easyCODA
# Usually, weighted LRA will be used 
# Here PLOT.LRA is used, not the plot function in ca
data(time)
time.LRA <- LRA(time)
par(mar=c(3.5,3.3,4,2), font.lab=2, mgp=c(2,0.5,0), cex.axis=0.8)
PLOT.LRA(time.LRA, map="contribution", main="No rescaling")
PLOT.LRA(time.LRA, map="contribution", rescale=0.5, 
         main="Rescaling columns by 0.5") 

# summary function from ca package
summary(time.LRA)

# LRA with amalgamation specified
LRA(time, amalg = list(nowork=c(4:6)))

# compute logratios of Vegetables data and perform unweighted PCA
veg.LR <- LR(veg)
veg.PCA <- PCA(veg.LR$LR, weight = FALSE)
PLOT.PCA(veg.PCA, map="asymmetric")

# perform weighted PCA on the ALRs, with denominator silicon
cups.ALR <- ALR(cups, denom = 1)
# ALR column weights automatically picked up by the PCA
cups.PCA <- PCA(cups.ALR)
# reverse second axis to agree with orientation of results in Chap. 9
PLOT.PCA(cups.PCA, map="contribution", rescale=0.2, axes.inv=c(1,-1))

# Perform RDA on the weighted CLRs of FishMorphology data set
# constraining variables are in matrix object vars
# assume fishm, mass, habitat and sex already available objects
# (see previous code)
logmass <- log(mass)
sexhab <- 2*(sex-1) + habitat
sexhab.names <- c("FL","FP","ML","MP")
rownames(fishm) <- sexhab.names[sexhab]

# create dummy variables for sexhab and create vars
sexhab.Z <- DUMMY(sexhab, catnames=sexhab.names)
vars <- cbind(logmass, sexhab.Z)

# perform RDA
fish.RDA <- RDA(CLR(fish), cov=vars)

# plot RDA - note that option indcat indicates which columns
# in the constraining variables object 'vars' are dummies
require(colorspace)
colcat <- rainbow_hcl(4, c=80)
par(mar=c(3.5,3.3,2.5,2.5), font.lab=2, mgp=c(2,0.5,0), cex.axis=0.8)
PLOT.RDA(fish.RDA, map="contribution", rescale=0.02, axes.inv=c(1,-1),
         indcat=2:5, colrows=colcat[sex_hab], colcats=colcat, 
         pchrows=rep(19, nrow(fishm)), cexs=c(1,1,1), fonts=c(2,4,2))
legend("topright", legend=c("female litoral (FL)", "male litoral (ML)", 
       "female pelagic (FP)", "male pelagic (MP)"), pch=19, 
       col=colcat, text.col=colcat, text.font=2)

# First compute row points as supplementary weighted averages
# of column standard coordinates, then substitute them in RDA object
fish.spc <- fishm %*% fish.RDA$colcoord
fish.RDA$rowpcoord <- fish.spc
par(mar=c(3.5,3.3,2.5,2.5), font.lab=2, mgp=c(2,0.5,0), cex.axis=0.8)
PLOT.RDA(fish.RDA, map="asymmetric", rescale=0.02, axes.inv=c(1,-1), 
         indcat=2:5, colrows=colcat[sexhab], colcats=colcat, 
         pchrows=rep(19, nrow(fishm)), cexs=c(1,1,1), fonts=c(2,4,2))

# add confidence ellipses (needs package ellipse)
require(ellipse)
CIplot_biv(fish.spc[,1], -fish.spc[,2],  alpha=0.99, group=sexhab,
           groupcols=colcat, shade=TRUE, add=TRUE, shownames=FALSE)

# Use of Procrustes analysis:
# three sets of principal coordinates: unweighted LRA, weighted LRA and CA
data(cups)
cups.ulra.rpc <- LRA(cups, weight=FALSE)$rowpcoord cups.wlra.rpc <- LRA(cups)$rowpcoord
cups.ca.rpc   <- CA(cups)$rowpcoord

# Procrustes correlations in full 10-D space
require(vegan)
# ... between unweighted and weighted LRA
protest(cups.ulra.rpc, cups.wlra.rpc, permutations=0)$t0
# [1] 0.7377072
# ... between weighted LRA and CA
protest(cups.wlra.rpc, cups.ca.rpc, permutations=0)$t0
# [1] 0.9962435

# Procrustes correlations in reduced 2-D space
# ... between unweighted and weighted LRA
protest(cups.ulra.rpc[,1:2], cups.wlra.rpc[,1:2], permutations=0)$t0
# [1] 0.6979616
# ... between weighted LRA and CA
protest(cups.wlra.rpc[,1:2], cups.ca.rpc[,1:2], permutations=0)$t0
# [1] 0.9938629

# Stepwise variable selection (random breaking of ties)
set.seed(2872)
cups.step <- STEP(cups, random=TRUE)
cups.step

# Stepwise variable selection (Procrustes breaking of ties)
STEP(cups)$names

# Including two-part amalgations as logratios
# compute amalgamations
amalgs <- matrix(0, nrow(cups), ncol=0.5*ncol(cups)*(ncol(cups)-1))
colnames(amalgs) <- 1:ncol(amalgs)
cupnames <- colnames(cups)
k <- 1
for(jj in 2:ncol(cups)) {
  for(j in 1:(jj-1)) {
    amalgs[,k] <- cups[,j]+cups[,jj]
    colnames(amalgs)[k] <- paste(cupnames[j], cupnames[jj], sep="+")
    k <- k+1
 }
}

# combine amalgamations with original parts
cups_amalgs <- cbind(cups, amalgs)

# perform stepwise selection and show some results
cups.step2 <- STEP(cups_amalgs, cups)
cups.step2$R2max
#  [1] 0.6454359 0.7755263 0.8828579 0.9428255 0.9732152 0.9859494 
#  [7] 0.9914720 0.9950723 0.9976510 0.9994871
cups.step2$names
#  [1] "Si/Ca+Na"    "Na/P+Sb"     "Al+Mg/Na+Sb" "Na/Ca+P"     "K+P/Mg+Sb"    
#  [6] "Al+Fe/Mg+Ti" "Al+P/Fe+Mn"  "Fe+K/Ti+P"   "Mn/Fe+Ti"    "Ti/P"



