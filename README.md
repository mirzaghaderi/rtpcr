


[![](https://cranlogs.r-pkg.org/badges/grand-total/rtpcr)](https://cran.rstudio.com/web/packages/rtpcr/index.html)
[![R package version](https://img.shields.io/github/r-package/v/mirzaghaderi/rtpcr)](
https://github.com/mirzaghaderi/rtpcr
)
[![CRAN version](https://www.r-pkg.org/badges/version/rtpcr)](
https://CRAN.R-project.org/package=rtpcr
)



# Multi-target multi-reference qPCR data analysis using 'rtpcr' package

# Overview

'rtpcr' handles amplification efficiency calculation, statistical analysis, and graphical representation of quantitative real-time PCR (qPCR) data using any number of references for any number of target genes. By accounting for amplification efficiency values, 'rtpcr' was developed using a general calculation method described by Ganger et al. (2017) and Taylor et al. (2019), covering the Livak and Pfaffl methods. Based on the experimental conditions, the functions of the 'rtpcr' package use t-test (for experiments with a two-level factor), analysis of variance (ANOVA), analysis of covariance (ANCOVA) or analysis of repeated measure data to analyse the relative expression (Delta Delta Ct method or Delta Ct method). The functions also provide standard errors and confidence intervals for means, apply statistical mean comparisons, and present significance.

## Improvement in the current GitHub version rtpcr_v2.1.1:

### Same data structure for all functions.
### No restriction for the number of target and reference genes in data.
### Analysis of all or a subset of genes.
### Graphic enhancement.



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

Sample data structure
```r
library(rtpcr)
data <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
data

factor1	factor2	block	Rep	PO	PO	PAO5	PAO5	ref1	ref1	ref2	ref2
R	L1	1	1	2	33.30	2	31.53	2	36.81	2	27.01
R	L1	1	2	2	33.39	2	31.57	2	30.71	2	27.17
R	L1	2	3	2	33.34	2	31.53	2	30.03	2	28.53
R	L1	2	4	2	30.34	2	30.50	2	30.03	2	28.95
R	L2	1	1	2	32.73	2	31.30	2	29.91	2	27.64
R	L2	1	2	2	   NA	2	32.55	2	30.05	2	27.36
R	L2	2	3	2	32.64	2	31.92	2	29.99	2	28.46
R	L2	2	4	2	30.65	2	30.81	2	30.55	2	28.48
R	L3	1	1	2	33.48	2	33.3	2	30.1	2	27.72
R	L3	1	2	2	33.27	2	33.37	2	30.17	2	28.14
R	L3	2	3	2	33.32	2	33.35	2	30.51	2	29.61
R	L3	2	4	2	30.32	2	29.35	2	30.51	2	29.51
S	L1	1	1	2	26.85	2	26.94	2	30.73	2	27.7
S	L1	1	2	2	28.17	2	27.69	2	30.68	2	27.68
S	L1	2	3	2	27.99	2	27.34	2	30	2	   NA
S	L1	2	4	2	28.99	2	27.36	2	30.71	2	27.11
S	L2	1	1	2	30.41	2	28.77	2	30.58	2	28.7
S	L2	1	2	2	29.49	2	28.66	2	30.03	2	28.92
S	L2	2	3	2	   NA	2	28.71	2	31.06	2	28.43
S	L2	2	4	2	28.33	2	28.72	2	31.06	2	27.43
S	L3	1	1	2	29.03	2	30.61	2	31.14	2	27.42
S	L3	1	2	2	28.77	2	30.22	2	30.24	2	27.81
S	L3	2	3	2	28.84	2	30.49	2	30.11	2	28.29
S	L3	2	4	2	27.81	2	29.34	2	30.11	2	27.24
```



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

