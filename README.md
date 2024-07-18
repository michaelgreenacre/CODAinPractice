# CODAinPractice

## **COMPOSITIONAL DATA ANALYSIS IN PRACTICE**

This repository contains data files and R scripts for the book ***Compositional Data Analysis in Practice*** (Michael Greenacre, Chapman & Hall / CRC Press, 2018):

  https://www.crcpress.com/Compositional-Data-Analysis-in-Practice/Greenacre/p/book/9781138316430

as well as some other files related to articles on compositional data analysis.

  \
The **easyCODA** R package accompanies the book and is available on CRAN, presently version 0.34.
The package is still under development and the latest version can always be found on R-Forge, installing as follows from R:

  install.packages("easyCODA", repos="http://R-Forge.R-project.org")

The package does include the data sets, but the data files are given here as well in Excel or character format.

  \
**ERRATA in "Compositional Data Analysis in Practice":**

**CDAiP_typos.pdf**: For the readers of my book, there is a provisional list of errors that I and others have found since its publication.

**CDAiP_typos.rtf**: rich text format of the above

  \
**DATA SETS for Compositional Data Analysis in Practice:**

**Vegetables.txt**: Vegetables data set, data object 'veg' in easyCODA

**TimeBudget.txt**: Time budget data set, data object 'time' in easyCODA

**RomanCups.xls**:  Archaeometric data on Roman glass cups, data object 'cups' in easyCODA

**FishMorphology.txt**: Fish morphometric data, with three additional variables, data object 'fish' in easyCODA 

  \
**R SCRIPTS for Compositional Data Analysis in Practice:**

**easyCODA_script.R**: file of all the R commands in Appendix C, also with some slight corrections.

PLEASE REPORT ANY BUGS OR DIFFICULTIES WITH THE PACKAGE TO Michael Greenacre at michael.greenacre@gmail.com

-----------------------------------------------------------------------------------------------------------

## **R scripts and data sets from other publications**

  \
**ARTICLE**: ***"A comparison of amalgamation logratio balances and isometric logratio balances in compositional data analysis"***, by Michael Greenacre, Eric Grunsky and John Bacon-Shone, Computers and Geosciences (2020)

**CAGEOscript.R**: script related Greenacre, Grunsky & Bacon-Shone (2020) 

  \
**ARTICLE**: ***"Amalgamations are valid in compositional data analysis, can be used in agglomerative clustering, and their logratios have an inverse transformation"***, by Michael Greenacre, Applied Computing and Geosciences (2020)

**SLRscript.R**: script related to article 

  \
**ARTICLE**: ***"The selection and analysis of fatty acid ratios: A new approach for the univariate and multivariate analysis of fatty acid trophic markers in marine pelagic organisms"***, by Martin Graeve and Michael Greenacre, Limnology & Oceanography Methods (2020)

**amphipod_ratios.R**: script related to the article by Graeve & Greenacre (2020)

**amphipods.csv**: data set for R script above

**copepod_ratios.R**: script related to the article "The selection and analysis of fatty acid ratios: A new approach for the univariate and multivariate analysis of fatty acid trophic markers in marine pelagic organisms", by Martin Graeve and Michael Greenacre, Limnology & Oceanography Methods (2020)

**copepods.csv**: data set for R script above

  \
**ARTICLE**: ***"Compositional Data Analysis"***, by Michael Greenacre, Annual Reviews in Statistics and its Application (2021)

**ANNUALREVIEWSscript.R**: R script for  article by Greenacre (2021)

**copepods_TL.csv**: data set for R script above (same data set as for Graeve & Greenacre (2020), with added variable Total Lipids (TL))

**Baxter_OTU_table.txt**: microbiome data set for R script above. From Baxter et al. (2016)

**Baxter_metadata.txt**: metadata that goes with the OTU table. From Baxter et al. (2016)

  \
**ARTICLE**: ***"Making the most of expert knowledge to analyse archaeologicaldata: a case study on Parthian and Sasanian glazed pottery"***, by Jonathan Wood and Michael Greenacre, Archaeological and Anthropological Sciences (2021)

**Supplementary material** for the article by Wood & Greenacre (2021). Two zip files:

**Wood&Greenacre_CSV**: data files 

**Wood&Greenacre_CODE&FUNCTIONS**: R code and additional R functions

  \
**ARTICLE**: ***"Compositional data analysis of microbiome and any-omics datasets: a validation of the additive logratio transformation"***, by Michael Greenacre, Marina Martinez-Alvaro and Agustin Blasco, Frontiers in Microbiology (2021)

**Frontiers_ALR.R**: script for Greenacre et al. (2021)

**Frontiers_ALR_supplementary.R**: script for supplementary material of Greenacre et al. (2021)

**FINDALR.R**: function FINDALR to identify optimal ALR reference, without or with weights

**Rabbits.xlsx**: Rabbits data set (89 rows, 3937 columns), rabbits x microbial genes. The different case ID label style refers to the two laboratories that did the testing

**Baxter_OTU_table.txt**: microbiome data set (Baxter et al., 2016), used in supplementary material of Greenacre et al. (2021)

**Baxter_metadata.txt**: metadata of above data set (Baxter et al., 2016), used in supplementary material of Greenacre et al. (2021)

**Deng_vaginal_microbiome.txt**: Vaginal microbiome data (Deng et al., 2018), cited by Wu et al. (2021)

  \
**ARTICLE:** ***"Three approaches to supervised learning for compositional data with pairwise logratios"***, by Germà Coenders and Michael Greenacre (2022)

**Coenders&Greenacre_CODE.R**: R code for analysis of Crohn data (Crohn data available in R package selbal, as shown in code)

**STEPR.R**: function for stepwise selection of logratios for GLM models (this function is in the pre-release of easyCODA on RForge)

  \
**ARTICLE:** ***"Aitchison's Compositional Data Analysis 40 Years On: A Reappraisal"***, by Michael Greenacre, Eric Grunsky, John Bacon-Shone, Ionas Erb and Thomas Quinn (2022). There are two versions: the original Version 1, and the new revised Version 2

**tellus_Appendix.R**: script for reproducing the analysis of the Tellus geochemical data set in Appendix of Version 1

**tellus_CoDA_script.R**: script for reproducing the analysis of the Tellus geochemical data set in Version 2 (essentially the same as before)

**tellus.xrf.a.cation.txt**: Tellus cation data set

**singlecell_CoDA_script.R**: script for reproducing the analysis of the single cell genetic data set in Version 2

**SingleCell.RData**: R workspace containing all the data files for the single cell application in Section 6 of Version 2

  \
  **ARTICLE:** ***"GeoCoDA: Recognizing and Validating Structural Processes in Geochemical Data. A Workflow on Compositional Data Analysis in Lithogeochemistry"***, by Eric Grunsky, Michael Greenacre, and Bruce Kjarsgaard (2023). 

**kimberlite.cation.closed.txt**: kimberlite cation data (Grunsky and Kjarsgaard, 2008), samples closed to sum to 100%

**GeoCoDA_script.R**: script for reproducing the GeoCoDA workflow

 \
  **ARTICLE:** ***"A Comprehensive Workflow for Compositional Data Analysis in Archaeometry, with code in R"***, by Michael Greenacre and Jonathan Wood (2024). 

**Bronze.csv**: dataset of Chines ritual bronzes

**Bronze_script.R**: script for reproducing the workflow
 
