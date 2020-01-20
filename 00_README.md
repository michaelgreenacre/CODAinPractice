# CODAinPractice
COMPOSITIONAL DATA ANALYSIS IN PRACTICE

This repository contains data files and R scripts for the book Compositional Data Analysis in Practice (Michael Greenacre, Chapman & Hall / CRC Press, 2018):

  https://www.crcpress.com/Compositional-Data-Analysis-in-Practice/Greenacre/p/book/9781138316430

as well as some other files related to articles on compositional data analysis.

The easyCODA R package accompanies the book and is available on CRAN, presently version 0.31.
The package is still under development and the latest version can always be found on R-Forge, installing as follows from R:

  install.packages("easyCODA", repos="http://R-Forge.R-project.org")

(the latest version 0.31 has a few important bug fixes, especially to the function STEP, so please reinstall this version if you have an earlier one).

The package does include the data sets, but the data files are given here as well in Excel or character format.


ERRATA in "Compositional Data Analysis in Practice":

CDAiP_typos.pdf: For the readers of my book, there is a provisional list of errors that I and others have found since its publication.

CDAiP_typos.rtf: rich text format of the above


DATA SETS:

Vegetables.txt: Vegetables data set, data object 'veg' in easyCODA

TimeBudget.txt: Time budget data set, data object 'time' in easyCODA

RomanCups.xls:  Archaeometric data on Roman glass cups, data object 'cups' in easyCODA

FishMorphology.txt: Fish morphometric data, with three additional variables, data object 'fish' in easyCODA 


R SCRIPTS:

easyCODA_script.R: file of all the R commands in Appendix C, also with some slight corrections.

PLEASE REPORT ANY BUGS OR DIFFICULTIES WITH THE PACKAGE TO Michael Greenacre at michael.greenacre@gmail.com

-----------------------------------------------------------------------------------------------------------

R scripts and data sets from other publications

CAGEOscript.R: script related to the article "A comparison of amalgamation logratio balances and isometric logratio balances in compositional data analysis", by Michael Greenacre, Eric Grunsky and John Bacon-Shone, Computers and Geosciences (2020)

SLRscript.R: script related to article " Amalgamations are valid in compositional data analysis, can be used in agglomerative clustering, and their logratios have an inverse transformation", by Michael Greenacre, Applied Computing and Geosciences (2020)

amphipod_ratios.R: script related to the article "The selection and analysis of fatty acid ratios: A new approach for the univariate and multivariate analysis of fatty acid trophic markers in marine pelagic organisms", by Martin Graeve and Michael Greenacre, Limnology & Oceanography Methods (2020)

amphipods.csv: data set for R script above
