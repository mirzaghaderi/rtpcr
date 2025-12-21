


[![](https://cranlogs.r-pkg.org/badges/grand-total/rtpcr)](https://cran.rstudio.com/web/packages/rtpcr/index.html)
[![R package version](https://img.shields.io/github/r-package/v/mirzaghaderi/rtpcr)](
https://github.com/mirzaghaderi/rtpcr
)
[![CRAN version](https://www.r-pkg.org/badges/version/rtpcr)](
https://CRAN.R-project.org/package=rtpcr
)



# Multi-target multi-reference qPCR data analysis using 'rtpcr' package

# Overview

'rtpcr' handles amplification efficiency calculation, ΔΔCt or ΔCt expression analysis, and graphical representation of quantitative real-time PCR (qPCR) data using any number of references for any number of target genes. 'rtpcr' uses a general calculation method described by Ganger et al. (2017) and Taylor et al. (2019), covering the Livak and Pfaffl methods. Based on the experimental conditions, the functions of the 'rtpcr' package use t-test (for experiments with a two-level factor), analysis of variance (ANOVA, or ANCOVA), or analysis of repeated measure data. to analyze the relative expression. The functions also provide standard errors and confidence intervals for means, apply statistical mean comparisons.

## Improvement in the current GitHub version rtpcr_v2.1.1:

<span style="color: green;">**Same data structure for all functions.**</span>
 
<span style="color: green;">**No restriction for the number of target and reference genes in data.**</span>

<span style="color: green;">**Expression analysis of all or a subset of genes.**</span>
 
<span style="color: green;">**Graphic enhancement.**</span>


# Quick start
### Installing and loading

The current version of the `rtpcr` package can be installed from GitHub by running the following code in R:

```r
devtools::install_github("mirzaghaderi/rtpcr", build_vignettes = FALSE)

# Loading the package
library(rtpcr)
```

Or from CRAN which may install the previous version.
```r
# Installing from CRAN
install.packages("rtpcr")
```


Further information about how to use rtpcr package can be found 
<a href="https://cran.r-project.org/web/packages/rtpcr/vignettes/vignette.html">here </a>.


### Data structure
Input data structure is important and should be in wide format:
For analysis using `TTEST_DDCt`, `ANOVA_DCt`, and `ANOVA_DDCt`, the required column structure is:

1. Experimental condition columns (Factors, and block if available) 
2. Biological replicate information (if applicable)  
3. Target genes efficiency and Ct values (a pair column for each target gene)
5. Reference genes efficiency and Ct values (a pair column for each reference gene)

The package supports **one or more target gene(s) and reference gene(s)**, supplied as efficiency–Ct column pairs.  
**Reference gene columns must always appear last.** Each row represents a single biological replicate, corresponding to a non-repeated measures design.

### Sample input data
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
R	L3	1	1	2	33.48	2	33.3	2	30.10	2	27.72
R	L3	1	2	2	33.27	2	33.37	2	30.17	2	28.14
R	L3	2	3	2	33.32	2	33.35	2	30.51	2	29.61
R	L3	2	4	2	30.32	2	29.35	2	30.51	2	29.51
S	L1	1	1	2	26.85	2	26.94	2	30.73	2	27.70
S	L1	1	2	2	28.17	2	27.69	2	30.68	2	27.68
S	L1	2	3	2	27.99	2	27.34	2	30		2	   NA
S	L1	2	4	2	28.99	2	27.36	2	30.71	2	27.11
S	L2	1	1	2	30.41	2	28.77	2	30.58	2	28.78
S	L2	1	2	2	29.49	2	28.66	2	30.03	2	28.92
S	L2	2	3	2	   NA	2	28.71	2	31.06	2	28.43
S	L2	2	4	2	28.33	2	28.72	2	31.06	2	27.43
S	L3	1	1	2	29.03	2	30.61	2	31.14	2	27.42
S	L3	1	2	2	28.77	2	30.22	2	30.24	2	27.81
S	L3	2	3	2	28.84	2	30.49	2	30.11	2	28.29
S	L3	2	4	2	27.81	2	29.34	2	30.11	2	27.24
```

### Data structure for `REPEATED_DDCt` function

The `REPEATED_DDCt` function is intended for experiments with repeated observations (e.g. time-course data).  
The input data frame for `REPEATED_DDCt` must follow this structure:

1. The **first column** is `id`, a unique identifier for each individual  
2. Factor and block columns (if available), and the `time` variable 
3. Remaining columns contain efficiency and Ct values for target and reference genes.

Each row corresponds to one observation at a specific time point for a given individual.
### Sample input data for `REPEATED_DDCt` function
```r
data <- read.csv(system.file("extdata", "data_repeated_measure_2.csv", package = "rtpcr"))
data

id	treatment	time	Target	Ct_target	E_Ref	Ct_Ref
1	untreated	1	2	19.24	2	33.73
1	untreated	2	2	19.95	2	34.20
1	untreated	3	2	19.16	2	33.90
2	untreated	1	2	20.11	2	32.56
2	untreated	2	2	20.91	2	33.98
2	untreated	3	2	20.91	2	33.16
3	untreated	1	2	20.63	2	33.72
3	untreated	2	2	19.16	2	34.51
3	untreated	3	2	19.91	2	34.33
4	treated		1	2	18.92	2	32.77
4	treated		2	2	19.46	2	33.03
4	treated		3	2	15.73	2	32.95
5	treated		1	2	15.82	2	32.45
5	treated		2	2	17.56	2	33.24
5	treated		3	2	17.21	2	33.64
6	treated		1	2	19.84	2	31.62
6	treated		2	2	19.74	2	32.08
6	treated		3	2	18.09	2	33.40
```

### Functions
Different functions for DDCt and DCt analysis, and efficiency calculation!

```r
# Example

res <- ANOVA_DDCt(
  x = data,
  mainFactor.column = 1,
  NumOfFactors = 2,
  numberOfrefGenes = 1,
  block = "block",
  analyseAllTarget = TRUE)

df <- res$combinedFoldChange
```
### Output
As lot of outputs including lm models, ANOVA table, residuals, raw and expression tables are returned.
```r
Relative Expression
  contrast      RE  log2FC pvalue sig    LCL     UCL     se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  gene
1        R  1.0000  0.0000 1.0000     0.0000  0.0000 0.5506      0.6828      1.4647          0.0000          0.0000    PO
2   S vs R 11.6130  3.5377 0.0001 *** 4.4233 30.4888 0.2286      9.9115     13.6066          3.0193          4.1450    PO
3        R  1.0000  0.0000 1.0000     0.0000  0.0000 0.4815      0.7162      1.3962          0.0000          0.0000 GAPDH
4   S vs R  6.6852  2.7410 0.0001 *** 3.0687 14.5641 0.3820      5.1301      8.7118          2.1034          3.5719 GAPDH
5        R  1.0000  0.0000 1.0000     0.0000  0.0000 0.6928      0.6186      1.6164          0.0000          0.0000  ref2
6   S vs R  0.9372 -0.0936 0.9005     0.3145  2.7929 0.2414      0.7927      1.1079         -0.1107         -0.0792  ref2
```

### Sample plot output
```r
data <- read.csv(system.file("extdata", "data_2factorBlock.csv", package = "rtpcr"))
 a <- ANOVA_DCt(data, 
                NumOfFactors = 2,
                block = "Block", 
                numberOfrefGenes = 1)
 df <- a$combinedResults
 
 p1 <- plotTwoFactor(
   data = df,
   x_col = "factor2",
   y_col = "RE",
   group_col = "factor1",
   Lower.se_col = "Lower.se.RE",
   Upper.se_col = "Upper.se.RE",
   letters_col = "sig",
   letters_d = 0.2,
   fill_colors = c("aquamarine4", "gold2"),
   color = "black",
   alpha = 1,
   col_width = 0.7,
   dodge_width = 0.7,
   base_size = 16, 
   legend_position = c(0.2, 0.8))
 
p1
```

<p align="center">
<img src="inst/Rplot01.png" width="100%">
</p>

<p align="center">
<img src="inst/Rplot02.png" width="100%">
</p>


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

