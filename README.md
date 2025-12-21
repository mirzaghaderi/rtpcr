


[![](https://cranlogs.r-pkg.org/badges/grand-total/rtpcr)](https://cran.rstudio.com/web/packages/rtpcr/index.html)
[![R package version](https://img.shields.io/github/r-package/v/mirzaghaderi/rtpcr)](
https://github.com/mirzaghaderi/rtpcr
)
[![CRAN version](https://www.r-pkg.org/badges/version/rtpcr)](
https://CRAN.R-project.org/package=rtpcr
)
![GitHub Release Downloads](https://img.shields.io/github/downloads/mirzaghaderi/rtpcr/total?label=GitHub%20downloads)



# Multi-target multi-reference qPCR data analysis using 'rtpcr' package

# Overview

`rtpcr` handles amplification efficiency calculation, ΔΔCt or ΔCt expression analysis, and graphical representation of quantitative real-time PCR (qPCR) data using Ct (Cq) and amplification efficiency values as input. `rtpcr` uses a general calculation method described by Ganger et al. (2017) and Taylor et al. (2019), covering the Livak and Pfaffl methods. Based on the experimental conditions, the functions of the `rtpcr` package use t-test (for experiments with a two-level factor), analysis of variance (ANOVA, or ANCOVA), or analysis of repeated measure data. to analyze the relative expression. The functions also provide standard errors and confidence intervals for means, apply statistical mean comparisons.

## Improvement in the current GitHub version rtpcr_v2.1.1:

<span style="color: green;">**Same data structure for all functions.**</span>
 
<span style="color: green;">**No restriction for the number of target and reference genes in data.**</span>

<span style="color: green;">**Expression analysis of all or a subset of genes.**</span>
 
<span style="color: green;">**Graphic enhancement.**</span>


# Functions (Brief explanations)

| Function            | Description                                                  |
|---------------------|--------------------------------------------------------------|
| `ANOVA_DCt`         | ΔCt ANOVA analysis                                           |
| `ANOVA_DDCt`        | ΔΔCt ANOVA analysis                                          |
| `REPEATED_DDCt`     | ΔΔCt ANOVA analysis for repeated-measures data               |
| `TTEST_DDCt`        | ΔΔCt method *t*-test analysis                                |
| `plotOneFactor`     | Bar plot of gene expression for single-factor experiments    |
| `plotTwoFactor`     | Bar plot of gene expression for two-factor experiments       |
| `plotThreeFactor`   | Bar plot of gene expression for three-factor experiments     |
| `efficiency`        | Amplification efficiency statistics and standard curves      |
| `meanTech`          | Calculate mean of technical replicates                       |
| `multiplot`         | Combine multiple ggplot objects into a single layout          |




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


### Amplification Efficiency

The `efficiency` function calculates the amplification efficiency (E), slope, and $R^2$ statistics for genes. It takes a data frame where the first column contains the dilution ratios followed by the $C_t$ value columns of genes.

```{r eval= T}
# Applying the efficiency function
data <- read.csv(system.file("extdata", "data_efficiency.csv", package = "rtpcr"))
data
dilutions	Gene1	Gene2	Gene3
1.00	25.58	24.25	22.61
1.00	25.54	24.13	22.68
1.00	25.50	24.04	22.63
0.50	26.71	25.56	23.67
0.50	26.73	25.43	23.65
0.50	26.87	26.01	23.70
0.20	28.17	27.37	25.11
0.20	28.07	26.94	25.12
0.20	28.11	27.14	25.11
0.10	29.20	28.05	26.17
0.10	29.49	28.89	26.15
0.10	29.07	28.32	26.15
0.05	30.17	29.50	27.12
0.05	30.14	29.93	27.14
0.05	30.12	29.71	27.16
0.02	31.35	30.69	28.52
0.02	31.35	30.54	28.57
0.02	31.35	30.04	28.53
0.01	32.55	31.12	29.49
0.01	32.45	31.29	29.48
0.01	32.28	31.15	29.26

# Analysis
efficiency(data)

$Efficiency
   Gene     Slope        R2        E
1 Gene1 -3.388094 0.9965504 1.973110
2 Gene2 -3.528125 0.9713914 1.920599
3 Gene3 -3.414551 0.9990278 1.962747

$Slope_compare
$contrasts
 contrast          estimate    SE df t.ratio p.value
 C2H2.26 - C2H2.01   0.1400 0.121 57   1.157  0.4837
 C2H2.26 - GAPDH     0.0265 0.121 57   0.219  0.9740
 C2H2.01 - GAPDH    -0.1136 0.121 57  -0.938  0.6186
```



### Data structure for `ANOVA_DDCt`, `ANOVA_DCt` and `TTEST_DDCt` functions
Input data structure is important and should be in wide format:
For analysis using `TTEST_DDCt`, `ANOVA_DCt`, and `ANOVA_DDCt`, the required column structure is:

1. Experimental condition columns (Factors, and block if available) 
2. Biological replicate information (if applicable)  
3. Target genes efficiency and Ct values (a pair column for each target gene)
5. Reference genes efficiency and Ct values (a pair column for each reference gene)

The package supports **one or more target gene(s) and reference gene(s)**, supplied as efficiency–Ct column pairs.  
**Reference gene columns must always appear last.** Each row represents a single biological replicate, corresponding to a non-repeated measures design. A sample input data is presented below.

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

Each row corresponds to one observation at a specific time point for a given individual. Below is an example:
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

### Analysis 
Different functions for ΔΔCt and ΔCt analysis, and efficiency calculation!

```r
# Example
data <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))

res <- ANOVA_DDCt(
  x = data,
  mainFactor.column = 1,
  NumOfFactors = 2,
  numberOfrefGenes = 1,
  block = "block",
  analyseAllTarget = TRUE)
```


### Output
A lot of outputs including relative expression table, lm models, residuals, raw data and ANOVA table for each gene can be accessed.
The expression table of all genes is returned by `res$combinedFoldChange`. Other outpus for each gene (i.e. for the PO gene) can be obtained as followe:

| Per_gene Output    | Code                                              |
|--------------------|-------------------------------------------------------|
| ANOVA table        | `res$perGene[["PO"]]$ANOVA_table`                     |
| lm ANOVA           | `res$perGene[["PO"]]$lm_ANOVA`                        |
| lm ANCOVA          | `res$perGene[["PO"]]$lm_ANCOVA`                       |
| Residuals          | `resid(res$perGene[["PO"]]$lm_ANOVA)`                 |


```r
df <- res$combinedFoldChange
df
Relative Expression
gene   contrast	      RE  log2FC pvalue sig    LCL     UCL     se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
PO            R	  1.0000  0.0000 1.0000     0.0000  0.0000 0.5506      0.6828      1.4647          0.0000          0.0000
PO       S vs R  11.6130  3.5377 0.0001 *** 4.4233 30.4888 0.2286      9.9115     13.6066          3.0193          4.1450
GAPDH         R	  1.0000  0.0000 1.0000     0.0000  0.0000 0.4815      0.7162      1.3962          0.0000          0.0000
GAPDH    S vs R	  6.6852  2.7410 0.0001 *** 3.0687 14.5641 0.3820      5.1301      8.7118          2.1034          3.5719
ref2          R	  1.0000  0.0000 1.0000     0.0000  0.0000 0.6928      0.6186      1.6164          0.0000          0.0000
ref2     S vs R	  0.9372 -0.0936 0.9005     0.3145  2.7929 0.2414      0.7927      1.1079         -0.1107         -0.0792
```

### Sample plot output
```r
data <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
#Perform analysis first
res <- ANOVA_DCt(
  data,
  NumOfFactors = 3,
  numberOfrefGenes = 1,
  block = NULL)
  
df <- res$combinedResults
 df
 # Generate three-factor bar plot
 p <- plotThreeFactor(
  df,
  x_col = "SA",       
  y_col = "log2FC",       
  group_col = "Type",   
  facet_col = "Conc",    
  Lower.se_col = "Lower.se.log2FC",
  Upper.se_col = "Upper.se.log2FC",
  letters_col = "sig",
  letters_d = 0.3,
  col_width = 0.7, 
  dodge_width = 0.7,
  fill_colors = c("palegreen3", "skyblue"),
  color = "black",
  base_size = 14, 
  alpha = 1,
  legend_position = c(0.1, 0.2))

library(ggplot2)
p + theme(
  panel.border = element_rect(color = "black", linewidth = 0.5))
```

<p align="center">
<img src="inst/Rplot02.png" width="100%">
</p>



# How to edit ouptput graphs?
the `rtpcr` plot functions (`plotOneFactor`, `plotTwoFactor`, and `plotThreeFactor`) create ggplot objects that can furtherbe edited by adding new layers:

| Task | Example Code |
|------|--------------|
| **Change y-axis label** | `p + ylab("Relative expression (ΔΔCt method)")` |
| **Add a horizontal reference line** | `p + geom_hline(yintercept = 0, linetype = "dashed")` |
| **Change y-axis limits** | `p + scale_y_continuous(limits = c(0, 20))` |
| **Relabel x-axis** | `p + scale_x_discrete(labels = c("A" = "Control", "B" = "Treatment"))` |
| **Change fill colors** | `p + scale_fill_brewer(palette = "Set2")` |


## A full graph Example
```{r eval= F, warning = F}
# Example 1
data <- read.csv(system.file("extdata", "data_2factorBlock.csv", package = "rtpcr"))
res <- ANOVA_DCt(data, 
      NumOfFactors = 2,
      block = "block",
      numberOfrefGenes = 1)

df <- res$combinedResults

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

library(ggplot2)
p1 + scale_y_continuous(expand = c(-1.5, +1.5)) + 
  theme(axis.text.x = element_text(size = 14, color = "black", angle = 45),
        axis.text.y = element_text(size = 14,color = "black", angle = 0, hjust = 0.5)) +
  theme(legend.text = element_text(colour = "black", size = 14),
        legend.background = element_rect(fill = "transparent"))
```

<p align="center">
<img src="inst/Rplot01.png" width="100%">
</p>


# Checking normality of residuals

If the residuals from a `t.test` or an `lm` or and `lmer` object are not normally distributed, the significance results might be violated. In such cases, non-parametric tests can be used. For example, the Mann–Whitney test (also known as the Wilcoxon rank-sum test), implemented via `wilcox.test()`, is an alternative to t.test, and `kruskal.test()` is an alternative to one-way analysis of variance. These tests assess differences between population medians using independent samples. However, the `t.test` function (along with the `TTEST_DDCt` function described above) includes the `var.equal` argument. When set to `FALSE`, perform `t.test` under the unequal variances hypothesis. Residuals for `lm` (from `ANOVA_DCt` and `ANOVA_DDCt` functions) and `lmer` (from `REPEATED_DDCt` function) objects can be extracted and plotted as follow:

```r
data1 <- read.csv(system.file("extdata", "data_1factor.csv", package = "rtpcr"))
res <- ANOVA_DCt(data1,
                 NumOfFactors = 1,
                 numberOfrefGenes = 1, 
                 block = NULL)

Extracting residuals
residuals <-  resid(res$perGene[["PO"]]$lmCRD)

shapiro.test(residuals) 
par(mfrow = c(1,2))
plot(residuals)
qqnorm(residuals)
qqline(residuals, col = "red")


data2 <- read.csv(system.file("extdata", "data_repeated_measure_1.csv", package = "rtpcr"))
res3 <- REPEATED_DDCt(
  data2,
  NumOfFactors = 1,
  numberOfrefGenes = 1,
  factor = "time",
  calibratorLevel = "1",
  block = NULL
)
residuals <- resid(res3$perGene[["Target"]]$lm)
```

# Mean of technical replicates
Calculating the mean of technical replicates and generating an output table suitable for subsequent ANOVA analysis can be accomplished using the `meanTech` function. The input dataset must follow the column structure illustrated in the example data below. Columns used for grouping should be explicitly specified via the `groups` argument of the `meanTech` function.

```{r eval= T}
# See example input data frame:
data <- read.csv(system.file("extdata", "data_withTechRep.csv", package = "rtpcr"))
data

# Calculating mean of technical replicates
meanTech(data, groups = 1:4)
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

