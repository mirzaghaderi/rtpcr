# Delta Delta Ct method t-test analysis

The `TTEST_DDCt` function performs fold change expression analysis based
on the \\\Delta \Delta C_T\\ method using Student's t-test. It supports
analysis of one or more target genes evaluated under two experimental
conditions (e.g. control vs treatment).

## Usage

``` r
TTEST_DDCt(
  x,
  numberOfrefGenes,
  Factor.level.order = NULL,
  paired = FALSE,
  var.equal = TRUE,
  p.adj = "none",
  order = "none"
)
```

## Arguments

- x:

  A data frame containing experimental conditions, biological
  replicates, and amplification efficiency and Ct values for target and
  reference genes. The number of biological replicates must be equal
  across genes. If this is not true, or there are `NA` values use
  `ANODA_DDCt` function for independent samples or `REPEATED_DDCt` for
  paired samples. See the package vignette for details on the required
  data structure.

- numberOfrefGenes:

  Integer specifying the number of reference genes used for
  normalization.

- Factor.level.order:

  Optional character vector specifying the order of factor levels. If
  `NULL`, the first level of the factor column is used as the
  calibrator.

- paired:

  Logical; if `TRUE`, a paired t-test is performed.

- var.equal:

  Logical; if `TRUE`, equal variances are assumed and a pooled variance
  estimate is used. Otherwise, Welch's t-test is applied.

- p.adj:

  Method for p-value adjustment. One of `"holm"`, `"hochberg"`,
  `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`, or `"none"`. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- order:

  Optional character vector specifying the order of genes in the output
  plot.

## Value

A list with the following components:

- Result:

  Table containing RE values, log2FC, p-values, significance codes,
  confidence intervals, standard errors, and lower/upper SE limits.

## Details

Relative expression values are computed using one or more reference
genes for normalization. Both paired and unpaired experimental designs
are supported.

Paired samples in quantitative PCR refer to measurements collected from
the same individuals under two different conditions (e.g. before vs
after treatment), whereas unpaired samples originate from different
individuals in each condition. Paired designs allow within-individual
comparisons and typically reduce inter-individual variability.

The function returns numerical summaries as well as bar plots based on
either relative expression (RE) or log2 fold change (log2FC).

All the functions for relative expression analysis (including
\`TTEST_DDCt()\`, \`WILCOX_DDCt()\`, \`ANOVA_DDCt()\`,
\`ANCOVA_DDCt()\`, \`REPEATED_DDCt()\`, and \`ANOVA_DCt()\`) return the
relative expression table which include fold change and corresponding
statistics. The output of \`ANOVA_DDCt()\`, \`ANCOVA_DDCt()\`,
\`ANCOVA_DDCt()\`, \`REPEATED_DDCt()\`, and \`ANOVA_DCt()\` also include
lm models, residuals, raw data and ANOVA table for each gene.

The expression table returned by \`TTEST_DDCt()\`, \`WILCOX_DDCt()\`,
\`ANOVA_DDCt()\`, \`ANCOVA_DDCt()\`, and \`REPEATED_DDCt()\` functions
include these columns: gene (name of target genes), contrast (calibrator
level and contrasts for which the relative expression is computed), RE
(relative expression or fold change), log2FC (log(2) of relative
expression or fold change), pvalue, sig (per-gene significance), LCL
(95% lower confidence level), UCL (95% upper confidence level), se
(standard error of mean calculated from the weighted delta Ct values of
each of the main factor levels), Lower.se.RE (The lower limit error bar
for RE which is 2^(log2(RE) - se)), Upper.se.RE (The upper limit error
bar for RE which is 2^(log2(RE) + se)), Lower.se.log2FC (The lower limit
error bar for log2 RE), and Upper.se.log2FC (The upper limit error bar
for log2 RE)

## References

LivakKJ, Schmittgen TD (2001). Analysis of Relative Gene Expression Data
Using Real-Time Quantitative PCR and the Double Delta CT Method.
*Methods*, 25(4), 402–408. doi:10.1006/meth.2001.1262

Ganger MT, Dietz GD, and Ewing SJ (2017). A common base method for
analysis of qPCR data and the application of simple blocking in qPCR
experiments. *BMC Bioinformatics*, 18, 1–11.

Taylor SC, Nadeau K, Abbasi M, Lachance C, Nguyen M, Fenrich, J. (2019).
The ultimate qPCR experiment: producing publication quality,
reproducible data the first time. *Trends in Biotechnology*, 37,
761-774.

Yuan JS, Reed A, Chen F, Stewart N (2006). Statistical Analysis of
Real-Time PCR Data. *BMC Bioinformatics*, 7, 85.

## Author

Ghader Mirzaghaderi

## Examples

``` r
# Example data structure
data1 <- read.csv(system.file("extdata", "data_ttest18genes.csv", package = "rtpcr"))

# Unpaired t-test
TTEST_DDCt(
  data1,
  paired = FALSE,
  var.equal = TRUE,
  numberOfrefGenes = 1)
#> *** 18 target(s) using 1 reference gene(s) was analysed!
#> *** The Control level was used as calibrator.
#>      gene     ddCt      RE     LCL      UCL  pvalue      se   log2FC
#> 1  ANGPT1  0.22175 0.85752 0.16876  4.35732 0.82472 0.28410 -0.22175
#> 2  ANGPT2 -1.45575 2.74299 1.69011  4.45178 0.00222 0.25051  1.45575
#> 3    CCL2  3.14850 0.11277 0.07430  0.17116 0.00001 0.23541 -3.14850
#> 4    CCL5  1.45925 0.36368 0.08320  1.58964 0.14436 0.53962 -1.45925
#> 5    CSF2  1.33800 0.39557 0.06001  2.60767 0.27417 0.80423 -1.33800
#> 6    FGF2 -1.15125 2.22106 1.46374  3.37021 0.00339 0.20688  1.15125
#> 7    IL1A -3.07550 8.42981 5.07767 13.99494 0.00005 0.24035  3.07550
#> 8    IL1B  0.59650 0.66136 0.30574  1.43058 0.23771 0.45245 -0.59650
#> 9     IL6 -0.21300 1.15910 0.70866  1.89583 0.49048 0.23038  0.21300
#> 10    IL8  2.62550 0.16205 0.03805  0.69020 0.02186 0.79514 -2.62550
#> 11  PDGFA  1.34350 0.39406 0.06212  2.49970 0.26354 0.88569 -1.34350
#> 12  PDGFB  0.27450 0.82674 0.21596  3.16494 0.74057 0.58071 -0.27450
#> 13   TGFA  2.29900 0.20320 0.07659  0.53911 0.00715 0.39389 -2.29900
#> 14   TGFB  1.56325 0.33839 0.13092  0.87464 0.03149 0.53238 -1.56325
#> 15    TNF -0.02975 1.02084 0.42354  2.46045 0.95612 0.44239  0.02975
#> 16  VEGFA  1.69900 0.30800 0.12631  0.75105 0.01785 0.40141 -1.69900
#> 17  VEGFB -0.29875 1.23008 0.63993  2.36445 0.46755 0.15382  0.29875
#> 18  VEGFC -1.89400 3.71664 1.40612  9.82382 0.01631 0.31851  1.89400
#>    Lower.se.RE Upper.se.RE sig Lower.se.log2FC Upper.se.log2FC
#> 1      0.70425     1.04416            -0.27001        -0.18211
#> 2      2.30576     3.26313  **         1.22370         1.73180
#> 3      0.09579     0.13276 ***        -3.70654        -2.67448
#> 4      0.25020     0.52864            -2.12115        -1.00390
#> 5      0.22653     0.69075            -2.33644        -0.76623
#> 6      1.92434     2.56353  **         0.99745         1.32876
#> 7      7.13615     9.95799 ***         2.60353         3.63304
#> 8      0.48332     0.90497            -0.81623        -0.43592
#> 9      0.98802     1.35979             0.18156         0.24988
#> 10     0.09339     0.28119   *        -4.55588        -1.51304
#> 11     0.21328     0.72809            -2.48231        -0.72714
#> 12     0.55279     1.23645            -0.41054        -0.18354
#> 13     0.15465     0.26700  **        -3.02073        -1.74971
#> 14     0.23397     0.48941   *        -2.26095        -1.08085
#> 15     0.75125     1.38716             0.02189         0.04043
#> 16     0.23319     0.40680   *        -2.24403        -1.28635
#> 17     1.10568     1.36848             0.26854         0.33236
#> 18     2.98036     4.63482   *         1.51879         2.36190

# With amplification efficiencies
data2 <- read.csv(system.file("extdata", "data_1factor_one_ref_Eff.csv", package = "rtpcr"))

TTEST_DDCt(
  data2,
  numberOfrefGenes = 1)
#> *** 1 target(s) using 1 reference gene(s) was analysed!
#> *** The control level was used as calibrator.
#>     gene    ddCt      RE    LCL     UCL  pvalue      se   log2FC Lower.se.RE
#> 1 target 1.28974 0.40902 0.2388 0.70059 0.00994 0.19848 -1.28974     0.35645
#>   Upper.se.RE sig Lower.se.log2FC Upper.se.log2FC
#> 1     0.46935  **        -1.47996        -1.12397

# Two reference genes
data3 <- read.csv(system.file("extdata", "data_1factor_Two_ref.csv", package = "rtpcr"))
TTEST_DDCt(
  data3,
  numberOfrefGenes = 2)
#> *** 1 target(s) using 2 reference gene(s) was analysed!
#> *** The control level was used as calibrator.
#>   gene    ddCt      RE     LCL     UCL  pvalue      se   log2FC Lower.se.RE
#> 1 DER5 1.20386 0.43411 0.19154 0.98387 0.04727 0.42036 -1.20386     0.32438
#>   Upper.se.RE sig Lower.se.log2FC Upper.se.log2FC
#> 1     0.58095   *        -1.61108        -0.89957
```
