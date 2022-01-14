### script for paper by Coenders & Greenacre on stepwise selection of logratios
### Crohn data is in package selbal (so needs prior installation)

library(selbal)
data(Crohn)
dim(Crohn)
# [1] 975  49

### last column is the response variable, extract it then remove
y <- Crohn[,49]
Crohn <- Crohn[,-49]
table(y)
# CD  no
# 662 313

### relevel so "no" is the reference level
y <- relevel(y, ref="no")
table(y)
# no  CD
# 313 662

### zeros problem, and replacement using zCompositions package
sum(Crohn==0)
# [1] 13474
require(zCompositions)
Crohn.no0 <- cmultRepl(Crohn, output = "p-counts")
sum(Crohn.no0==0)
# [1] 0

### closure wasn't necessary for computing logratios, but we do it anyway 
### (.pro stands for profile, a term from correspondence analysis)
Crohn.no0.pro <- Crohn.no0 / rowSums(Crohn.no0)

### The timings below are on a Toshiba Satellite Pro
### ------------------------------------------------------------------------------------------------------
### Method 1
time0 <- Sys.time()
Crohn.step1 <- STEPR(Crohn.no0.pro, y, method=1, family="binomial")
# [1] "Criterion increases when 9-th ratio enters"
time1 <- Sys.time()
difftime(time1, time0, units="secs")
# Time difference of 53.00978 sec

Crohn.step1$logLik
# [1] 1064.7252  999.0034  966.1923  939.5918  916.7610  900.1576  883.5464  870.6115  863.2098
Crohn.step1$deviance
# [1] 1064.7252  999.0034  966.1923  939.5918  916.7610  900.1576  883.5464  870.6115  863.2098
Crohn.step1$AIC
# [1] 1068.7252 1005.0034  974.1923  949.5918  928.7610  914.1576  899.5464  888.6115  883.2098
Crohn.step1$BIC
# [1] 1078.4901 1019.6508  993.7221  974.0040  958.0557  948.3347  938.6059  932.5534  932.0341
Crohn.step1$Bonferroni
#  1086.1513 1031.1426 1009.0445  993.1570  981.0393  975.1489  969.2507  967.0289  970.3402
Crohn.step1$names
#  [1] "g__Roseburia/g__Streptococcus"             "f__Peptostreptococcaceae_g__/g__Dialister"
#  [3] "g__Bacteroides/g__Dorea"                   "g__Prevotella/g__Aggregatibacter"        
#  [5] "g__Adlercreutzia/g__Lachnospira"           "o__Lactobacillales_g__/g__Roseburia"      
#  [7] "g__Oscillospira/o__Clostridiales_g__"      "g__Sutterella/g__Bilophila"                  
 
summary(glm(as.factor(y) ~ Crohn.step1$logratios, family="binomial"))
# Coefficients:
#                                                                  Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                                                     3.10968    0.32698   9.510  < 2e-16 ***
#   Crohn.step1$logratiosg__Roseburia/g__Streptococcus             -0.30220    0.03147  -9.603  < 2e-16 ***
#   Crohn.step1$logratiosf__Peptostreptococcaceae_g__/g__Dialister -0.16175    0.02181  -7.417 1.20e-13 ***
#   Crohn.step1$logratiosg__Bacteroides/g__Dorea                   -0.23930    0.03724  -6.426 1.31e-10 ***
#   Crohn.step1$logratiosg__Prevotella/g__Aggregatibacter          -0.10080    0.02195  -4.591 4.41e-06 ***
#   Crohn.step1$logratiosg__Adlercreutzia/g__Lachnospira            0.11579    0.02734   4.235 2.28e-05 ***
#   Crohn.step1$logratioso__Lactobacillales_g__/g__Streptococcus    0.14816    0.03638   4.072 4.66e-05 ***
#   Crohn.step1$logratiosg__Oscillospira/o__Clostridiales_g__       0.16881    0.04255   3.967 7.26e-05 ***
#   Crohn.step1$logratiosg__Sutterella/g__Bilophila                 0.08729    0.02462   3.545 0.000392 ***
#   ---

# Null deviance: 1223.90  on 974  degrees of freedom
# Residual deviance:  870.61  on 966  degrees of freedom
# AIC: 888.61
 
BIC(glm(as.factor(y) ~ Crohn.step1$logratios, family="binomial"))
# [1] 932.5534

-2*logLik(glm(as.factor(y) ~ Crohn.step1$logratios, family="binomial"))
# 'log Lik.' 870.6115 (df=9)

### ------------------------------------------------------------------------------------------------------
### Method 2
time0 <- Sys.time()
Crohn.step2 <- STEPR(Crohn.no0.pro, y, method=2, family="binomial")
# [1] "Criterion increases when 9-th ratio enters"
time1 <- Sys.time()
difftime(time1, time0, units="secs")
# Time difference of 65.47664 secs

Crohn.step2$logLik
# [1] 1064.7252  999.0034  966.1923  939.5918  916.7610  902.0850  888.9433  877.6979  869.1014
Crohn.step2$deviance
# [1] 1064.7252  999.0034  966.1923  939.5918  916.7610  902.0850  888.9433  877.6979  869.1014
Crohn.step2$AIC
# [1] 1068.7252 1005.0034  974.1923  949.5918  928.7610  916.0850  904.9433  895.6979  889.1014
Crohn.step2$BIC
# [1] 1078.4901 1019.6508  993.7221  974.0040  958.0557  950.2621  944.0028  939.6398  937.9257
Crohn.step2$Bonferroni
# [1] 1086.1512 1031.1424 1009.0443  993.1568  981.0390  977.0760  974.6473  974.1149  976.2314
Crohn.step2$names
# [1] "g__Roseburia/g__Streptococcus"               "f__Peptostreptococcaceae_g__/g__Dialister"  
# [3] "g__Dorea/g__Bacteroides"                     "g__Aggregatibacter/g__Prevotella"          
# [5] "g__Lachnospira/g__Adlercreutzia"             "o__Clostridiales_g__/f__Ruminococcaceae_g__"
# [7] "g__Sutterella/g__Bilophila"                  "g__Oscillospira/g__Faecalibacterium"        
         
summary(glm(as.factor(y) ~ Crohn.step2$logratios, family="binomial"))
# Coefficients:
#                                                                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                                       2.26931    0.24982   9.084  < 2e-16 ***
# Crohn.step2$logratiosg__Roseburia/g__Streptococcus               -0.24440    0.02908  -8.405  < 2e-16 ***
# Crohn.step2$logratiosf__Peptostreptococcaceae_g__/g__Dialister   -0.17017    0.02165  -7.860 3.85e-15 ***
# Crohn.step2$logratiosg__Dorea/g__Bacteroides                      0.22717    0.03709   6.124 9.12e-10 ***
# Crohn.step2$logratiosg__Aggregatibacter/g__Prevotella             0.10867    0.02222   4.890 1.01e-06 ***
# Crohn.step2$logratiosg__Lachnospira/g__Adlercreutzia             -0.11390    0.02756  -4.133 3.58e-05 ***
# Crohn.step2$logratioso__Clostridiales_g__/f__Ruminococcaceae_g__ -0.25529    0.06415  -3.979 6.91e-05 ***
# Crohn.step2$logratiosg__Sutterella/g__Bilophila                   0.08440    0.02449   3.447 0.000568 ***
# Crohn.step2$logratiosg__Oscillospira/g__Faecalibacterium          0.10879    0.03326   3.271 0.001073 **
# ---
#     Null deviance: 1223.9  on 974  degrees of freedom
# Residual deviance:  877.7  on 966  degrees of freedom
# AIC: 895.7

BIC(glm(as.factor(y) ~ Crohn.step2$logratios, family="binomial"))
# [1] 939.6398

### -------------------------------------------------------------------------------------------------------
### Method 3
time0 <- Sys.time()
Crohn.step3 <- STEPR(Crohn.no0.pro, y, method=3, family="binomial")
# [1] "Criterion increases when 10-th ratio enters"
time1 <- Sys.time()
difftime(time1, time0, units="secs")
# Time difference of 6.293941 secs

Crohn.step3$logLik
#  [1] 1064.7252 1021.5936  998.5788  978.2177  965.1117  942.0897  926.4469  911.4926  898.6042  890.9743
Crohn.step3$deviance
#  [1] 1064.7252 1021.5936  998.5788  978.2177  965.1117  942.0897  926.4469  911.4926  898.6042  890.9743
Crohn.step3$AIC
# [1] 1068.7252 1027.5936 1006.5788  988.2177  977.1117  956.0897  942.4469  929.4926  918.6042  912.9743
Crohn.step3$BIC
# [1] 1078.4901 1042.2410 1026.1085 1012.6299 1006.4064  990.2668  981.5064  973.4345  967.4286  966.6811
Crohn.step3$Bonferroni
# [1] 1086.151 1053.733 1041.431 1031.783 1029.390 1017.081 1012.151 1007.910 1005.734 1008.817
Crohn.step3$names
# [1] "g__Roseburia/g__Streptococcus"                 "g__Dialister/g__Streptococcus"                
# [3] "f__Peptostreptococcaceae_g__/g__Streptococcus" "o__Lactobacillales_g__/g__Streptococcus"      
# [5] "g__Bacteroides/g__Streptococcus"               "g__Dorea/g__Streptococcus"                    
# [7] "g__Adlercreutzia/g__Streptococcus"             "g__Aggregatibacter/g__Streptococcus"          
# [9] "g__Prevotella/g__Streptococcus"
   
summary(glm(as.factor(y) ~ Crohn.step3$logratios, family="binomial"))
# Coefficients:
#                                                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                                 3.72618    0.46039   8.094 5.80e-16 ***
# foo$logratiosg__Roseburia/g__Streptococcus                 -0.33749    0.03833  -8.805  < 2e-16 ***
# foo$logratiosg__Dialister/g__Streptococcus                  0.14069    0.02624   5.362 8.25e-08 ***
# foo$logratiosf__Peptostreptococcaceae_g__/g__Streptococcus -0.20652    0.03238  -6.377 1.81e-10 ***
# foo$logratioso__Lactobacillales_g__/g__Streptococcus        0.14197    0.03965   3.581 0.000342 ***
# foo$logratiosg__Bacteroides/g__Streptococcus               -0.27915    0.04810  -5.804 6.49e-09 ***
# foo$logratiosg__Dorea/g__Streptococcus                      0.20213    0.04393   4.601 4.20e-06 ***
# foo$logratiosg__Adlercreutzia/g__Streptococcus              0.15111    0.03600   4.197 2.70e-05 ***
# foo$logratiosg__Aggregatibacter/g__Streptococcus            0.13782    0.03284   4.196 2.71e-05 ***
# foo$logratiosg__Prevotella/g__Streptococcus                -0.09203    0.02581  -3.565 0.000363 ***
# ---

#     Null deviance: 1223.9  on 974  degrees of freedom
# Residual deviance:  898.6  on 965  degrees of freedom
# AIC: 918.6

BIC(glm(as.factor(y) ~ Crohn.step3$logratios, family="binomial"))
# [1] 967.4286

### -------------------------------------------------------------------------------------------------------
### Method 3 again, with Roseburia in denominator
### (overwrites Crohn.step3, this is the ALR we want)

Crohn.step3 <- STEPR(Crohn.no0.pro, y, method=3, family="binomial", denom=32)
# [1] "Criterion increases when 10-th ratio enters"
summary(glm(as.factor(y) ~ Crohn.step3$logratios, family="binomial"))
# Coefficients:
#                                                                Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                                     3.72618    0.46039   8.094 5.80e-16 ***
# Crohn.step3$logratiosg__Streptococcus/g__Roseburia              0.14147    0.04378   3.232 0.001231 **
# Crohn.step3$logratiosg__Dialister/g__Roseburia                  0.14069    0.02624   5.362 8.25e-08 ***
# Crohn.step3$logratiosf__Peptostreptococcaceae_g__/g__Roseburia -0.20652    0.03238  -6.377 1.81e-10 ***
# Crohn.step3$logratioso__Lactobacillales_g__/g__Roseburia        0.14197    0.03965   3.581 0.000342 ***
# Crohn.step3$logratiosg__Bacteroides/g__Roseburia               -0.27915    0.04810  -5.804 6.49e-09 ***
# Crohn.step3$logratiosg__Dorea/g__Roseburia                      0.20213    0.04393   4.601 4.20e-06 ***
# Crohn.step3$logratiosg__Adlercreutzia/g__Roseburia              0.15111    0.03600   4.197 2.70e-05 ***
# Crohn.step3$logratiosg__Aggregatibacter/g__Roseburia            0.13782    0.03284   4.196 2.71e-05 ***
# Crohn.step3$logratiosg__Prevotella/g__Roseburia                -0.09203    0.02581  -3.565 0.000363 ***
# ---

# Null deviance: 1223.9  on 974  degrees of freedom
# Residual deviance:  898.6  on 965  degrees of freedom
# AIC: 918.6

BIC(glm(as.factor(y) ~ Crohn.step3$logratios, family="binomial"))
# [1] 967.4286

### -------------------------------------------------------------------------------------------------------
### Method 3 checking how Egge would enter at Step 9 - here with Stre as the denominator
### Top 10 asked for at 9th step
time0 <- Sys.time()
Crohn.step3a <- STEPR(Crohn.no0.pro, y, nsteps=1, top=10, method=3, family="binomial",
                      previous=Crohn.step3$logratios[,1:8], denom=39)
time1 <- Sys.time()
difftime(time1, time0, units="secs")
# Time difference of 0.374877 secs

Crohn.step3a$deviance.top
# [1] 898.6042 903.3561 904.4669 904.4930 905.2139 905.8167 906.1235 906.6195 906.7808 906.8320
Crohn.step3a$AIC.top
# [1] 918.6042 923.3561 924.4669 924.4930 925.2139 925.8167 926.1235 926.6195 926.7808 926.8320
Crohn.step3a$Bonferroni.top
# [1] 1005.734 1010.486 1011.597 1011.623 1012.344 1012.947 1013.253 1013.749 1013.911 1013.962

Crohn.step3a$ratios.top
#                                       row col
# g__Prevotella/g__Streptococcus         31  39
# g__Eggerthella/g__Streptococcus         9  39   <--- Egge/Stre is second from top
# g__Lachnospira/g__Streptococcus        33  39
# g__Oscillospira/g__Streptococcus       26  39
# g__Bilophila/g__Streptococcus          48  39
# g__Anaerostipes/g__Streptococcus       20  39
# g__Actinomyces/g__Streptococcus        23  39
# o__Clostridiales_g__/g__Streptococcus  34  39
# g__Sutterella/g__Streptococcus          4  39
# g__Faecalibacterium/g__Streptococcus   10  39

### do GLM with 8 logratios and Egge/Stre as the 9th
foo <- as.data.frame(list(previous=Crohn.step3$logratios[,1:8], Egge_Stre=Crohn.step3a$logratios.top[,2]))  
summary(glm(as.factor(y) ~ ., family="binomial", data=foo))
# Coefficients:
#                                                        Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                             4.31839    0.47121   9.164  < 2e-16 ***
# previous.g__Roseburia.g__Streptococcus                 -0.35449    0.03818  -9.284  < 2e-16 ***
# previous.g__Dialister.g__Streptococcus                  0.12704    0.02627   4.836 1.33e-06 ***
# previous.f__Peptostreptococcaceae_g__.g__Streptococcus -0.20447    0.03213  -6.364 1.97e-10 ***
# previous.o__Lactobacillales_g__.g__Streptococcus        0.12651    0.03899   3.244 0.001177 **
# previous.g__Bacteroides.g__Streptococcus               -0.31649    0.04859  -6.513 7.36e-11 ***
# previous.g__Dorea.g__Streptococcus                      0.16845    0.04449   3.786 0.000153 ***
# previous.g__Adlercreutzia.g__Streptococcus              0.13286    0.03644   3.646 0.000267 ***
# previous.g__Aggregatibacter.g__Streptococcus            0.12457    0.03218   3.871 0.000108 ***
# Egge_Stre                                               0.09049    0.03197   2.831 0.004640 **
# ---

#     Null deviance: 1223.90  on 974  degrees of freedom
# Residual deviance:  903.36  on 965  degrees of freedom
# AIC: 923.36

### -------------------------------------------------------------------------------------------------------
### Method 3 with Egge (9th part) entering at Step 9 - have to specify Rose as denominator
### (use the step3 analysis where Rose was the denominator reference)
### Again ask for top 10 at step 9
Crohn.step3b <- STEPR(Crohn.no0.pro, y, nsteps=1, top=10, method=3, family="binomial",
                      previous=Crohn.step3$logratios[,1:8], denom=32)
Crohn.step3b$deviance.top
#  [1] 898.6042 903.3561 904.4669 904.4930 905.2139 905.8167 906.1235 906.6195 906.7808 906.8320

Crohn.step3b$Bonferroni.top
#  [1] 1005.734 1010.486 1011.597 1011.623 1012.344 1012.947 1013.253 1013.749 1013.911 1013.962

Crohn.step3b$ratios.top
#                                   row col
# g__Prevotella/g__Roseburia         31  32
# g__Eggerthella/g__Roseburia         9  32   <--- Egge/Rose is second from top
# g__Lachnospira/g__Roseburia        33  32
# g__Oscillospira/g__Roseburia       26  32
# g__Bilophila/g__Roseburia          48  32
# g__Anaerostipes/g__Roseburia       20  32
# g__Actinomyces/g__Roseburia        23  32
# o__Clostridiales_g__/g__Roseburia  34  32
# g__Sutterella/g__Roseburia          4  32
# g__Faecalibacterium/g__Roseburia   10  32

### do GLM with 8 logratios (w.r.t. Rose) and Egge/Rose as the 9th
foo <- as.data.frame(list(previous=Crohn.step3$logratios[,1:8], Egge_Rose=Crohn.step3b$logratios.top[,2]))  
summary(glm(as.factor(y) ~ ., family="binomial", data=foo))
# Coefficients:
#                                                    Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                         4.31839    0.47121   9.164  < 2e-16 ***
# previous.g__Streptococcus.g__Roseburia              0.10552    0.04384   2.407 0.016084 *  
# previous.g__Dialister.g__Roseburia                  0.12704    0.02627   4.836 1.33e-06 ***
# previous.f__Peptostreptococcaceae_g__.g__Roseburia -0.20447    0.03213  -6.364 1.97e-10 ***
# previous.o__Lactobacillales_g__.g__Roseburia        0.12651    0.03899   3.244 0.001177 **
# previous.g__Bacteroides.g__Roseburia               -0.31649    0.04859  -6.513 7.36e-11 ***
# previous.g__Dorea.g__Roseburia                      0.16845    0.04449   3.786 0.000153 ***
# previous.g__Adlercreutzia.g__Roseburia              0.13286    0.03644   3.646 0.000267 ***
# previous.g__Aggregatibacter.g__Roseburia            0.12457    0.03218   3.871 0.000108 ***
# Egge_Rose                                           0.09049    0.03197   2.831 0.004640 **
# ---

#     Null deviance: 1223.90  on 974  degrees of freedom
# Residual deviance:  903.36  on 965  degrees of freedom
# AIC: 923.36

BIC(glm(as.factor(y) ~ ., family="binomial", data=foo))
# [1] 972.1805

### -------------------------------------------------------------------------------------------------------
### Method 3 moving on with Egge/Rose included in previous, always with Rose in denominator
### Just ask for 3 steps more
previous <- as.data.frame(list(first8=Crohn.step3$logratios[,1:8], Egge_Rose=Crohn.step3b$logratios.top[,2]))
Crohn.step3c <- STEPR(Crohn.no0.pro, y, nsteps=3, top=1, method=3, family="binomial", criterion=NA,
                      previous=previous, denom=32)
Crohn.step3c$BIC
# [1] 969.2846 967.5029 967.5860
Crohn.step3c$ratios
#                            row col
# g__Prevotella/g__Roseburia  31  32
# g__Sutterella/g__Roseburia   4  32
# g__Bilophila/g__Roseburia   48  32

### Do GLM with original first 8 logratios, then Egge/Rose and finally the first 2 steps above
foo <- as.data.frame(list(previous=Crohn.step3$logratios[,1:8], Egge_Rose=Crohn.step3b$logratios.top[,2], Crohn.step3c$logratios[,1:2]))  
summary(glm(as.factor(y) ~ ., family="binomial", data=foo))
# Coefficients:
#                                                    Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                         4.18945    0.48982   8.553  < 2e-16 ***
# previous.g__Streptococcus.g__Roseburia              0.13616    0.04510   3.019 0.002533 **
# previous.g__Dialister.g__Roseburia                  0.12830    0.02668   4.809 1.52e-06 ***
# previous.f__Peptostreptococcaceae_g__.g__Roseburia -0.20869    0.03264  -6.393 1.63e-10 ***
# previous.o__Lactobacillales_g__.g__Roseburia        0.13630    0.03992   3.415 0.000639 ***
# previous.g__Bacteroides.g__Roseburia               -0.36592    0.05538  -6.608 3.91e-11 ***
# previous.g__Dorea.g__Roseburia                      0.18475    0.04502   4.104 4.07e-05 ***
# previous.g__Adlercreutzia.g__Roseburia              0.13409    0.03665   3.659 0.000253 ***
# previous.g__Aggregatibacter.g__Roseburia            0.12489    0.03331   3.750 0.000177 ***
# Egge_Rose                                           0.08228    0.03294   2.498 0.012485 *  
# g__Prevotella.g__Roseburia                         -0.09226    0.02683  -3.439 0.000583 ***
# g__Sutterella.g__Roseburia                          0.09635    0.03292   2.927 0.003424 **
# ---

#     Null deviance: 1223.90  on 974  degrees of freedom
# Residual deviance:  884.91  on 963  degrees of freedom
# AIC: 908.91

BIC(glm(as.factor(y) ~ ., family="binomial", data=foo))
# [1] 967.5029
