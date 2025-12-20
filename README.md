


[![](https://cranlogs.r-pkg.org/badges/grand-total/rtpcr)](https://cran.rstudio.com/web/packages/rtpcr/index.html)

  


# Multi-target multi-reference qPCR data analysis using 'rtpcr' package

# Overview

'rtpcr' handles amplification efficiency calculation, statistical analysis, and graphical representation of quantitative real-time PCR (qPCR) data using any number of references for any number of target genes. By accounting for amplification efficiency values, 'rtpcr' was developed using a general calculation method described by Ganger et al. (2017) and Taylor et al. (2019), covering the Livak and Pfaffl methods. Based on the experimental conditions, the functions of the 'rtpcr' package use t-test (for experiments with a two-level factor), analysis of variance (ANOVA), analysis of covariance (ANCOVA) or analysis of repeated measure data to analyse the relative expression (Delta Delta Ct method or Delta Ct method). The functions also provide standard errors and confidence intervals for means, apply statistical mean comparisons, and present significance.

![My Plot](inst/Rplot01.png)

![My Plot](inst/Rplot02.png)

# Installing and loading

The rtpcr package and source code are available for download from CRAN website (https://www.r-project.org) under GPL-3 license. The `rtpcr` package can be installed by running the following code in R:

```r
# Installing from CRAN
install.packages("rtpcr")

# Loading the package
library(rtpcr)
```


Or install the development version from GitHub with:
```r
devtools::install_github("mirzaghaderi/rtpcr", build_vignettes = TRUE)
```

Further information about how to use rtpcr package can be found 
<a href="https://cran.r-project.org/web/packages/rtpcr/vignettes/vignette.html">here </a>.



# Contact 
Email: gh.mirzaghaderi at uok.ac.ir

# Citation
```r
citation("rtpcr")

To cite the package ‘rtpcr’ in publications, please use:

  Ghader Mirzaghaderi (2025). rtpcr: A package for statistical analysis and graphical
  presentation of qPCR data in R. PeerJ 13:e20185. https://doi.org/10.7717/peerj.20185

A BibTeX entry for LaTeX users is

  @Article{,
    title = {rtpcr: A package for statistical analysis and graphical presentation of qPCR data in R},
    author = {Ghader Mirzaghaderi},
    journal = {PeerJ},
    volume = {13},
    pages = {e20185},
    year = {2025},
    doi = {10.7717/peerj.20185},
  }
```

