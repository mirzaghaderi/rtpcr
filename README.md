


[![](https://cranlogs.r-pkg.org/badges/grand-total/rtpcr)](https://cran.rstudio.com/web/packages/rtpcr/index.html)

  


# rtpcr: a package for statistical analysis and graphical presentation of qPCR data in R
### Author: Ghader Mirzaghaderi

# Overview

'rtpcr' handles amplification efficiency calculation, statistical analysis, and graphical representation of quantitative real-time PCR (qPCR) data based on any number of specified reference genes. By accounting for amplification efficiency values, 'rtpcr' was developed using a general calculation method described by Ganger et al. (2017) and Taylor et al. (2019), covering the Livak and Pfaffl methods. Based on the experimental conditions, the functions of the 'rtpcr' package use t-test (for experiments with a two-level factor), analysis of variance (ANOVA), analysis of covariance (ANCOVA) or analysis of repeated measure data to analyse the relative expression (Delta Delta Ct method or Delta Ct method). The functions also provide standard errors and confidence intervals for means, apply statistical mean comparisons, and present significance.

# Installing and loading

The rtpcr package and source code are available for download from CRAN website (https://www.r-project.org) under GPL-3 license. The `rtpcr` package can be installed by running the following code in R:

```r
# Installing from CRAN
install.packages("rtpcr")

# Loading the package
library(rtpcr)
```


Further information about how to use rtpcr package can be found 
<a href="https://cran.r-project.org/web/packages/rtpcr/vignettes/vignette.html">here </a>.



# Contact 
Email: gh.mirzaghaderi at uok.ac.ir


# References
Livak Kenneth J, and Schmittgen TD. 2001. Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR and the Double Delta CT Method. Methods 25 (4). <a href="https://doi.org/10.1006/meth.2001.1262">doi.org/10.1006/meth.2001.1262</a>.


Ganger MT, Dietz GD, and Ewing SJ. 2017. A common base method for analysis of qPCR data and the application of simple blocking in qPCR experiments. BMC bioinformatics 18, 1-11. <a href="https://doi.org/10.1186/s12859-017-1949-5">doi.org/10.1186/s12859-017-1949-5</a>.

Pfaffl MW, Horgan GW, Dempfle L. 2002. Relative expression software tool (REST©) for group-wise comparison and statistical analysis of relative expression results in real-time PCR. Nucleic acids research 30, e36-e36. <a href="https://doi.org/10.1093/nar/30.9.e36">doi.org/10.1093/nar/30.9.e36</a>.


Taylor SC, Nadeau K, Abbasi M, Lachance C, Nguyen M, Fenrich, J. 2019. The ultimate qPCR experiment: producing publication quality, reproducible data the first time. Trends in Biotechnology, 37(7), 761-774<a href="https://doi.org/10.1016/j.tibtech.2018.12.002">doi.org/10.1016/j.tibtech.2018.12.002</a>.


Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006. Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85). <a href="https://doi.org/10.1186/1471-2105-7-85">doi.org/10.1186/1471-2105-7-85</a>.


.
