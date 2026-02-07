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
  set_missing_target_Ct_to_40 = FALSE
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

- set_missing_target_Ct_to_40:

  If `TRUE`, missing target gene Ct values become 40; if `FALSE`
  (default), they become NA.

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
#>      gene     ddCt      RE   log2FC     LCL      UCL      se Lower.se.RE
#> 1  ANGPT1  0.22175 0.85752 -0.22175 0.16876  4.35732 0.95843     0.44130
#> 2  ANGPT2 -1.45575 2.74299  1.45575 1.69011  4.45178 0.28552     2.25048
#> 3    CCL2  3.14850 0.11277 -3.14850 0.07430  0.17116 0.24599     0.09509
#> 4    CCL5  1.45925 0.36368 -1.45925 0.08320  1.58964 0.86965     0.19904
#> 5    CSF2  1.33800 0.39557 -1.33800 0.06001  2.60767 1.11192     0.18302
#> 6    FGF2 -1.15125 2.22106  1.15125 1.46374  3.37021 0.24586     1.87306
#> 7    IL1A -3.07550 8.42981  3.07550 5.07767 13.99494 0.29888     6.85245
#> 8    IL1B  0.59650 0.66136 -0.59650 0.30574  1.43058 0.45490     0.48250
#> 9     IL6 -0.21300 1.15910  0.21300 0.70866  1.89583 0.29009     0.94797
#> 10    IL8  2.62550 0.16205 -2.62550 0.03805  0.69020 0.85438     0.08963
#> 11  PDGFA  1.34350 0.39406 -1.34350 0.06212  2.49970 1.08923     0.18521
#> 12  PDGFB  0.27450 0.82674 -0.27450 0.21596  3.16494 0.79148     0.47765
#> 13   TGFA  2.29900 0.20320 -2.29900 0.07659  0.53911 0.57527     0.13638
#> 14   TGFB  1.56325 0.33839 -1.56325 0.13092  0.87464 0.55990     0.22955
#> 15    TNF -0.02975 1.02084  0.02975 0.42354  2.46045 0.51868     0.71255
#> 16  VEGFA  1.69900 0.30800 -1.69900 0.12631  0.75105 0.52555     0.21397
#> 17  VEGFB -0.29875 1.23008  0.29875 0.63993  2.36445 0.38528     0.94178
#> 18  VEGFC -1.89400 3.71664  1.89400 1.40612  9.82382 0.57308     2.49825
#>    Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1      1.66634        -0.43090        -0.11412 0.82472    
#> 2      3.34328         1.19437         1.77434 0.00222  **
#> 3      0.13374        -3.73383        -2.65493 0.00001 ***
#> 4      0.66453        -2.66636        -0.79862 0.14436    
#> 5      0.85495        -2.89185        -0.61906 0.27417    
#> 6      2.63373         0.97087         1.36515 0.00339  **
#> 7     10.37026         2.50002         3.78345 0.00005 ***
#> 8      0.90651        -0.81762        -0.43518 0.23771    
#> 9      1.41725         0.17420         0.26044 0.49048    
#> 10     0.29298        -4.74685        -1.45217 0.02186   *
#> 11     0.83841        -2.85844        -0.63146 0.26354    
#> 12     1.43096        -0.47512        -0.15859 0.74057    
#> 13     0.30277        -3.42541        -1.54300 0.00715  **
#> 14     0.49884        -2.30448        -1.06043 0.03149   *
#> 15     1.46250         0.02077         0.04262 0.95612    
#> 16     0.44336        -2.44568        -1.18028 0.01785   *
#> 17     1.60663         0.22873         0.39020 0.46755    
#> 18     5.52925         1.27311         2.81770 0.01631   *

# With amplification efficiencies
data2 <- read.csv(system.file("extdata", "data_1factor_one_ref_Eff.csv", package = "rtpcr"))

TTEST_DDCt(
  data2,
  numberOfrefGenes = 1)
#> *** 1 target(s) using 1 reference gene(s) was analysed!
#> *** The control level was used as calibrator.
#>     gene    ddCt      RE   log2FC    LCL     UCL      se Lower.se.RE
#> 1 target 1.28974 0.40902 -1.28974 0.2388 0.70059 0.27963     0.33695
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1     0.49651         -1.5656        -1.06249 0.00994  **

# Two reference genes
data3 <- read.csv(system.file("extdata", "data_1factor_Two_ref.csv", package = "rtpcr"))
TTEST_DDCt(
  data3,
  numberOfrefGenes = 2)
#> *** 1 target(s) using 2 reference gene(s) was analysed!
#> *** The control level was used as calibrator.
#>   gene    ddCt      RE   log2FC     LCL     UCL      se Lower.se.RE Upper.se.RE
#> 1 DER5 1.20386 0.43411 -1.20386 0.19154 0.98387 0.42515     0.32331     0.58289
#>   Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1        -1.61644        -0.89659 0.04727   *
```
