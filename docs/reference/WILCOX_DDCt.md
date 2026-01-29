# Delta Delta Ct method wilcox.test analysis

The `WILCOX_DDCt` function performs fold change expression analysis
based on the \\\Delta \Delta C_T\\ method using wilcox.test. It supports
analysis of one or more target genes evaluated under two experimental
conditions (e.g. control vs treatment).

## Usage

``` r
WILCOX_DDCt(
  x,
  numberOfrefGenes,
  Factor.level.order = NULL,
  paired = FALSE,
  p.adj = "none"
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

  Logical; if `TRUE`, a paired wilcox.test is performed.

- p.adj:

  Method for p-value adjustment. One of `"holm"`, `"hochberg"`,
  `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`, or `"none"`. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A table containing RE values, log2FC, p-values, significance, confidence
intervals, standard errors, and lower/upper SE limits.

## Details

Relative expression values are computed using reference gene(s) for
normalization. Both paired and unpaired experimental designs are
supported.

Paired samples in quantitative PCR refer to measurements collected from
the same individuals under two different conditions (e.g. before vs
after treatment), whereas unpaired samples originate from different
individuals in each condition. Paired designs allow within-individual
comparisons and typically reduce inter-individual variability.

The function returns expression table. The expression table returned by
\`TTEST_DDCt()\`, \`WILCOX_DDCt()\`, \`ANOVA_DDCt()\`,
\`ANCOVA_DDCt()\`, and \`REPEATED_DDCt()\` functions include these
columns: gene (name of target genes), contrast (calibrator level and
contrasts for which the relative expression is computed), RE (relative
expression or fold change), log2FC (log(2) of relative expression or
fold change), pvalue, sig (per-gene significance), LCL (95% lower
confidence level), UCL (95% upper confidence level), se (standard error
of mean calculated from the weighted delta Ct values of each of the main
factor levels), Lower.se.RE (The lower limit error bar for RE which is
2^(log2(RE) - se)), Upper.se.RE (The upper limit error bar for RE which
is 2^(log2(RE) + se)), Lower.se.log2FC (The lower limit error bar for
log2 RE), and Upper.se.log2FC (The upper limit error bar for log2 RE)

## References

Yuan, J. S., Reed, A., Chen, F., and Stewart, N. (2006). Statistical
Analysis of Real-Time PCR Data. *BMC Bioinformatics*, 7, 85.

## Author

Ghader Mirzaghaderi

## Examples

``` r
# Example data structure
data <- read.csv(system.file("extdata", "data_Yuan2006PMCBioinf.csv", package = "rtpcr"))

# Unpaired t-test
WILCOX_DDCt(
  data,
  paired = FALSE,
  numberOfrefGenes = 1)
#> *** 1 target(s) using 1 reference gene(s) was analysed!
#> *** The control level was used as calibrator.
#>     gene   ddCt      RE  log2FC     LCL     UCL      se Lower.se.RE Upper.se.RE
#> 1 target 0.6354 0.64376 -0.6354 0.54318 0.74603 0.07964     0.60919      0.6803
#>   Lower.se.log2FC Upper.se.log2FC pvalue sig
#> 1        -0.67146        -0.60127      0 ***


# Two reference genes
data2 <- read.csv(system.file("extdata", "data_1factor_Two_ref.csv", package = "rtpcr"))
WILCOX_DDCt(
  data2,
  numberOfrefGenes = 2,
  p.adj = "none")
#> *** 1 target(s) using 2 reference gene(s) was analysed!
#> *** The control level was used as calibrator.
#>   gene    ddCt      RE   log2FC     LCL     UCL      se Lower.se.RE Upper.se.RE
#> 1 DER5 1.36522 0.38817 -1.36522 0.25981 0.81119 0.42036     0.29006     0.51948
#>   Lower.se.log2FC Upper.se.log2FC pvalue sig
#> 1        -1.82702        -1.02015    0.1   .
```
