### ====================================================================================
### function to compute concordance measure (closeness to a 45 degree line through zero)
### i.e., closeness to exact numerical agreement
CONCORD <- function(x1,x2) {
# Lin's coefficient of concordance
  s1 <- sd(x1)
  s2 <- sd(x2)
  corr <- cor(x1,x2)
  concord <- 2*corr*s1*s2 / ((mean(x1) - mean(x2))^2 + s1^2 + s2^2)
  concord
}

### ====================================================================================
SUBCO <- function(data, nsub=100, propsub=1/3, weight=TRUE, power=NA, chi=FALSE,
                  procrustes=TRUE, distances=TRUE, correlations=TRUE, covariances=TRUE,
                  shannon=TRUE, simpson=TRUE, seed=123)  {             # removed MRM=TRUE, 

### requires package easyCODA, functions CONCORD and chiPower

### study of subcompositional coherence (unsupervised)
### data         = compositional data set (checked to be nonnegative), closed in the function
### nsub         = number of subcompositions tested (default 100)
### propsub      = proportion of the data used in subcomposition (default 1/3)
### weight       = weighted subcomposition selection using average part values (=TRUE, default)
### power        = power lambda for chiPower transform (default NA, no transformation)
###                set power=1 if chisquare normalization is required on raw compositions
###                otherwise raw compositional values are used
### The following are all TRUE by default, so need to be switched off if not required
### procrustes   = Procrustes correlations between subcompositional and compositional part geometries
### distances    = distance concordance between subcompositional and compositional part geometries
### correlations = concordance of correlations between subcompositional and compositional part geometries
### covariances  = concordance of covariances between subcompositional and compositional part geometries
### shannon      = Shannon indices (only for one subcomposition, the last one)
### simpson      = Simpson indices (only for one subcomposition, the last one)
### MRM          = MRM cross-product similarity (only for one subcomposition, the last one)
### seed         = seed for random number generator (subcompositions selected randomly)
 
procr.sub <- rep(NA,nsub)                                  # Procrustes correlations 
procr.chi.sub <- rep(NA,nsub)                              # Procrustes correlations (with chi-normalization)
concord.dist.sub <- rep(NA,nsub)                           # concordance of distances
concord.dist.chi.sub <- rep(0,nsub)                        # concordance of distances (with chi-normalization)
concord.cor.sub <- rep(NA,nsub)                            # concordance of correlations
concord.cov.sub <- rep(NA,nsub)                            # concordance of covariances
concord.cov.chi.sub <- rep(NA,nsub)                        # concordance of covariances (with chi-normalization)
correl.cov.sub  <- rep(NA,nsub)                            # correlation of covariances 
correl.cov.chi.sub  <- rep(NA,nsub)                        # correlation of covariances (with chi-normalization) 
shannon.sub <- matrix(NA, nrow(data), nsub)                # Shannon indices, as many as samples
simpson.sub <- matrix(NA, nrow(data), nsub)                # Simpson indices, as many as samples
# concord.MRM.sub <- rep(0,nsub)                             # concordance of values in MRM matrix
# MRM.sub <- array(0, c(nrow(data),nrow(data),nsub))         # microbiome relative matrices, as many square matrices as samples

### check and close the data
data.pro <- CLOSE(data)
data.cm <- colMeans(data.pro) 

### transformed data
data.trans <- data.pro
if(chi) data.trans <- sweep(data.pro, 2, sqrt(data.cm), FUN="/")
if(!is.na(power)) { 
### chiPower transformed (suppress the translation and rescaling in Box-Cox
### because here we are not trying to approximate the CLR)
  data.trans <- chiPower(data.pro, power=power, BoxCox=FALSE) 
}
### the following are needed for comparing with full compositions
### column principal coordinates for Procrustes
data.pro.PCA   <- PCA(data.pro, weight=FALSE)
if(chi) data.trans.PCA <- PCA(data.trans, weight=FALSE)
data.pro.cpc   <- data.pro.PCA$colpcoord
if(chi) data.trans.cpc <- data.trans.PCA$colpcoord
### distances between columns of data.pro and data.trans
data.pro.dist   <- as.matrix(dist(t(data.pro)))
if(chi) data.trans.dist <- as.matrix(dist(t(data.trans)))
### correlations between data, i.e. normalized by standard deviations (so chi-normalization makes no difference)
data.trans.corr <- cor(data.trans)
### covariances between data, using chi-square normalization on relative abundances
data.pro.cov <- cov(data.pro)
if(chi) data.trans.cov <- cov(data.trans)
### For Shannon and Simpson we need relative abundances
### Shannon on full composition
shannon     <- apply(CLOSE(data.pro), 1, diversity, index="shannon")
if(!is.na(power)) shannon <- apply(CLOSE(data.pro^power), 1, diversity, index="shannon")
shannon.max <- diversity(rep(1/ncol(data.pro), ncol(data.pro)))
shannon.adj <- shannon/shannon.max
### Simpson on full composition
simpson     <- apply(CLOSE(data.pro), 1, diversity, index="simpson")
if(!is.na(power)) simpson <- apply(CLOSE(data.pro^power), 1, diversity, index="simpson")
simpson.max <- diversity(rep(1/ncol(data.pro), ncol(data.pro)), index="simpson")
simpson.adj <- simpson/simpson.max
### cross-product similarity (MRM matrix)
# X <- data.trans
# XXt <- X %*% t(X)
# MRM <- XXt   # no division by p, the chi-square nromalization already includes adjusting for number of variables

### sampling will use weights proportional to average data relative abundances (in data.cm)

### -------------------------------------------------------------------------------------
### start of loop on nsub subcompositions and storing results 
### we use subcompositions with 1/3 the number of data
### sub.rand contains column numbers of selected subcompositional data 
set.seed(seed)
for(isub in 1:nsub) {	
  if(weight)  sub.rand <- sample(1:ncol(data.pro), round(propsub*ncol(data.pro)), prob=data.cm)   
  if(!weight) sub.rand <- sample(1:ncol(data.pro), round(propsub*ncol(data.pro)))
# The following would be for unweighted random sampling 
#  sub.rand <- sample(1:ncol(data.pro), sub.prop*ncol(data.pro))  
# create the closed subcomposition
  data.sub.pro <- CLOSE(data.pro[,sub.rand])
  data.sub.cm <- colMeans(data.sub.pro)
  data.sub.trans <- sweep(data.sub.pro, 2, sqrt(data.sub.cm), FUN="/")
  if(!is.na(power)) { 
# chiPower transformed (suppress the translation and rescaling in Box-Cox
# because here we are not trying to approximate the CLR)
  data.sub.trans <- chiPower(data.sub.pro, power=power, BoxCox=FALSE) 
  }
  data.sub.pro.PCA <- PCA(data.sub.pro, weight=FALSE)
  if(chi) data.sub.trans.PCA <- PCA(data.sub.trans, weight=FALSE)
  data.sub.pro.cpc <- data.sub.pro.PCA$colpcoord
  if(chi) data.sub.trans.cpc <- data.sub.trans.PCA$colpcoord
# Procrustes, omit last dimension corresponding to 0 variance
  procr.sub[isub] <- protest(data.pro.cpc[sub.rand, -ncol(data.trans.cpc)], data.sub.pro.cpc[,-ncol(data.sub.pro.cpc)], permutations=0)$t0                              
  if(chi) procr.chi.sub[isub] <- protest(data.trans.cpc[sub.rand, -ncol(data.trans.cpc)], data.sub.trans.cpc[,-ncol(data.sub.trans.cpc)], permutations=0)$t0
# distances, measure is Spearman correlation between compositional and subcompositional distances
  data.sub.pro.dist <- as.matrix(dist(t(data.sub.pro)))
  if(chi) data.sub.trans.dist <- as.matrix(dist(t(data.sub.trans)))
#  dist.sub[isub]  <- cor(as.dist(data.trans.dist[sub.rand,sub.rand]), as.dist(data.sub.trans.dist), method="spearman")
  concord.dist.sub[isub] <- CONCORD(as.dist(data.pro.dist[sub.rand,sub.rand]), as.dist(data.sub.pro.dist))
  if(chi) concord.dist.chi.sub[isub] <- CONCORD(as.dist(data.trans.dist[sub.rand,sub.rand]), as.dist(data.sub.trans.dist))
# correlations, and then concordance between correlations
  data.sub.trans.corr <- cor(data.sub.trans)
  concord.cor.sub[isub] <- CONCORD(as.dist(data.trans.corr[sub.rand,sub.rand]),as.dist(data.sub.trans.corr))
# covariances (with chi-square normalization), and then concordance between covariances
  data.sub.pro.cov <- cov(data.sub.pro)
  if(chi) data.sub.trans.cov <- cov(data.sub.trans)
  data.pro.cov.vec <- c(diag(data.pro.cov[sub.rand,sub.rand]),as.dist(data.pro.cov[sub.rand,sub.rand]))
  if(chi) data.trans.cov.vec <- c(diag(data.trans.cov[sub.rand,sub.rand]),as.dist(data.trans.cov[sub.rand,sub.rand]))
  data.sub.pro.cov.vec <- c(diag(data.sub.pro.cov),as.dist(data.sub.pro.cov))
  if(chi) data.sub.trans.cov.vec <- c(diag(data.sub.trans.cov),as.dist(data.sub.trans.cov))
  concord.cov.sub[isub] <- CONCORD(data.pro.cov.vec, data.sub.pro.cov.vec)
  if(chi) concord.cov.chi.sub[isub] <- CONCORD(data.trans.cov.vec, data.sub.trans.cov.vec)
  correl.cov.sub[isub] <- cor(data.pro.cov.vec, data.sub.pro.cov.vec)
  if(chi) correl.cov.chi.sub[isub] <- cor(data.trans.cov.vec, data.sub.trans.cov.vec)
# shannon
  shannon.sub.max    <- diversity(rep(1/ncol(data.sub.pro), ncol(data.sub.pro), index="shannon")) 
  shannon.sub.adj    <-  apply(data.sub.pro, 1, diversity, index="shannon") / shannon.sub.max
  shannon.sub[,isub] <- shannon.sub.adj
# simpson
  simpson.sub.max    <- diversity(rep(1/ncol(data.sub.pro), ncol(data.sub.pro)), index="simpson") 
  simpson.sub.adj    <- apply(data.sub.pro, 1, diversity, index="simpson") / simpson.sub.max
  simpson.sub[,isub] <- simpson.sub.adj
# MRM matrix, using chi-square normalized subcomposition, weights in data.sub.cm, then concordance
#   Xsub <- data.sub.trans
#   MRM.sub[,,isub] <- Xsub %*% t(Xsub) # no division by p, the chi-square nromalization already includes adjusting for number of variables 
#   concord.MRM.sub[isub] <- CONCORD(as.numeric(MRM), as.numeric(Xsub %*% t(Xsub)))
}
### ---------end of loop-----------------------------------------------------------------
  return(list(procr            = procr.sub,
              procr.chi        = procr.chi.sub,
              concord.dist     = concord.dist.sub,
              concord.dist.chi = concord.dist.chi.sub,
              concord.cor      = concord.cor.sub,
              concord.cov      = concord.cov.sub,
              concord.cov.chi  = concord.cov.chi.sub,
              correl.cov       = correl.cov.sub,
              correl.cov.chi   = correl.cov.chi.sub,
              shannon          = shannon.adj,
              shannon.sub      = shannon.sub,
              simpson          = simpson.adj,
              simpson.sub      = simpson.sub))
#             MRM          = MRM,
#             MRM.sub      = MRM.sub,
#             concord.MRM  = concord.MRM.sub
}


