# Delta Delta Ct repeated measure analysis

`REPEATED_DDCt` function performs Delta Delta Ct (ddCt) method analysis
of observations repeatedly taken over different time courses. Data may
be obtained over time from a uni- or multi-factorial experiment.

## Usage

``` r
REPEATED_DDCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  mainFactor.column,
  block,
  mainFactor.level.order = NULL,
  p.adj = "none",
  analyseAllTarget = TRUE
)
```

## Arguments

- x:

  The input data frame containing experimental design columns, target
  gene E/Ct column pairs, and reference gene E/Ct column pairs.
  Reference gene columns must be located at the end of the data frame.

- numOfFactors:

  Integer. Number of experimental factor columns (excluding `rep` and
  optional `block`).

- numberOfrefGenes:

  Integer. Number of reference genes.

- mainFactor.column:

  Integer. Column index of the factor (commonly `"time"`) for which
  relative expression is calculated.

- block:

  Character or `NULL`. Name of the blocking factor column. When a qPCR
  experiment is done in multiple qPCR plates, variation resulting from
  the plates may interfere with the actual amount of gene expression.
  One solution is to conduct each plate as a randomized block so that at
  least one replicate of each treatment and control is present on a
  plate. Block effect is usually considered as random and its
  interaction with any main effect is not considered.

- mainFactor.level.order:

  Optional character vector specifying the order of levels for the main
  factor. If `NULL`, the first observed level is used as the calibrator.
  If provided, the first element of the vector is used as the calibrator
  level.

- p.adj:

  Method for p-value adjustment. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- analyseAllTarget:

  Logical or character. If `TRUE` (default), all target genes are
  analysed. Alternatively, a character vector specifying the names
  (names of their Efficiency columns) of target genes to be analysed.

## Value

An object containing expression table, lm model, residuals, raw data and
ANOVA table for each gene:

- ddCt expression table along with per-gene statistical comparison
  outputs:

  `object$relativeExpression`

- ANOVA table:

  `object$perGene$gene_name$ANOVA_table`

- lm ANOVA:

  `object$perGene$gene_name$lm`

- lm_formula:

  `object$perGene$gene_name$lm_formula`

- Residuals:

  `resid(object$perGene$gene_name$lm)`

## Details

ddCt analysis of repeated measure data is performed for the
`mainFactor.column` based on a full model factorial experiment.

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
level and contrasts for which the relative expression is computed), ddCt
(mean of weighted delta delta Ct values), RE (relative expression or
fold change = 2^-ddCt), log2FC (log(2) of relative expression or fold
change), pvalue, sig (per-gene significance), LCL (95% lower confidence
level), UCL (95% upper confidence level), se (standard error of mean
calculated from the weighted delta Ct values of each of the main factor
levels), Lower.se.RE (The lower limit error bar for RE which is
2^(log2(RE) - se)), Upper.se.RE (The upper limit error bar for RE which
is 2^(log2(RE) + se)), Lower.se.log2FC (The lower limit error bar for
log2 RE), and Upper.se.log2FC (The upper limit error bar for log2 RE)

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

## Examples

``` r
data1 <- read.csv(system.file("extdata", "data_repeated_measure_1.csv", package = "rtpcr"))
REPEATED_DDCt(
  data1,
  numOfFactors = 1,
  numberOfrefGenes = 1,
  mainFactor.column = 1,
  block = NULL)
#> 
#> Relative Expression
#>     gene       contrast     ddCt      RE   log2FC  pvalue sig    LCL      UCL
#> 1 Target              1  0.00000 1.00000  0.00000 1.00000     0.0000  0.00000
#> 2 Target time2 vs time1  0.22333 0.85658 -0.22333 0.81660     0.0923  7.94918
#> 3 Target time3 vs time1 -2.23333 4.70219  2.23333 0.06847   . 0.5067 43.63677
#>        se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1 1.40507     0.37760     2.64830         0.00000         0.00000
#> 2 0.97527     0.43570     1.68405        -0.43907        -0.11360
#> 3 0.55411     3.20256     6.90403         1.52108         3.27911
#> The 1 level was used as calibrator.



data2 <- read.csv(system.file("extdata", "data_repeated_measure_2.csv", package = "rtpcr"))
REPEATED_DDCt(
  data2,
  numOfFactors = 2,
  numberOfrefGenes = 1,
  mainFactor.column = 2, 
  block = NULL,
  p.adj = "none")
#> 
#> Relative Expression
#>     gene       contrast     ddCt      RE  log2FC  pvalue sig     LCL     UCL
#> 1 Target              1  0.00000 1.00000 0.00000 1.00000     0.00000 0.00000
#> 2 Target time2 vs time1 -0.32833 1.25556 0.32833 0.55402     0.43805 3.59873
#> 3 Target time3 vs time1 -1.34667 2.54324 1.34667 0.03509   * 0.88731 7.28951
#>        se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1 0.70359     0.61404     1.62855         0.00000         0.00000
#> 2 0.53234     0.86814     1.81588         0.22702         0.47486
#> 3 0.70742     1.55752     4.15279         0.82472         2.19894
#> The 1 level was used as calibrator.
          
```
