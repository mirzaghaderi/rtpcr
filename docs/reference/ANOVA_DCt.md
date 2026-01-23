# Delta Ct ANOVA analysis

Performs Delta Ct (dCt) analysis of the data from a 1-, 2-, or 3-factor
experiment. Per-gene statistical grouping is also performed for all
treatment (T) combinations.

## Usage

``` r
ANOVA_DCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  block,
  alpha = 0.05,
  p.adj = "none",
  analyseAllTarget = TRUE
)
```

## Arguments

- x:

  The input data frame containing experimental design columns, target
  gene E/Ct column pairs, and reference gene E/Ct column pairs.
  Reference gene columns must be located at the end of the data frame.
  See "Input data structure" in vignettes for details about data
  structure.

- numOfFactors:

  Integer. Number of experimental factor columns (excluding `rep` and
  optional `block`).

- numberOfrefGenes:

  Integer. Number of reference genes. Each reference gene must be
  represented by two columns (E and Ct).

- block:

  Character. Block column name or `NULL`. When a qPCR experiment is done
  in multiple qPCR plates, variation resulting from the plates may
  interfere with the actual amount of gene expression. One solution is
  to conduct each plate as a randomized block so that at least one
  replicate of each treatment and control is present on a plate. Block
  effect is usually considered as random and its interaction with any
  main effect is not considered.

- alpha:

  statistical level for comparisons

- p.adj:

  Method for p-value adjustment. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- analyseAllTarget:

  Logical or character. If `TRUE` (default), all detected target genes
  are analysed. Alternatively, a character vector specifying the names
  (names of their Efficiency columns) of target genes to be analysed.

## Value

An object containing expression table, lm models, ANOVA table,
residuals, and raw data for each gene:

- dCt expression table for all treatment combinations along with the
  per-gene statistical grouping:

  `object$relativeExpression`

- ANOVA table for treatments:

  `object$perGene$gene_name$ANOVA_T`

- ANOVA table factorial:

  `object$perGene$gene_name$ANOVA_factorial`

- lm ANOVA for tratments:

  `object$perGene$gene_name$lm_T`

- lm ANOVA factorial:

  `object$perGene$gene_name$lm_factorial`

- Residuals:

  `resid(object$perGene$gene_name$lm_T)`

## Details

The function returns analysis of variance components and the expression
table which include these columns: gene (name of target genes), factor
columns, dCt (mean weighted delta Ct for each treatment combination), RE
(relative expression = 2^-dCt), log2FC (log(2) of relative expression),
LCL (95% lower confidence level), UCL (95% upper confidence level), se
(standard error of the mean calculated from the weighted delta Ct (wDCt)
values of each treatment combination), Lower.se.RE (The lower limit
error bar for RE which is 2^(log2(RE) - se)), Upper.se.RE (The upper
limit error bar for RE which is 2^(log2(RE) + se)), Lower.se.log2FC (The
lower limit error bar for log2 RE), Upper.se.log2FC (The upper limit
error bar for log2 RE), and sig (per-gene significance grouping
letters).

## Examples

``` r
data <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
res <- ANOVA_DCt(
  data,
  numOfFactors = 3,
  numberOfrefGenes = 1,
  block = NULL)
#> 
#> Relative Expression
#>    gene Type Conc SA      dCt      RE   log2FC     LCL     UCL      se
#> 1    PO    S    H A2 -2.37667 5.19335  2.37667 8.11969 3.32167 0.13094
#> 2    PO    S    H A1 -1.57000 2.96905  1.57000 4.64204 1.89900 0.05508
#> 3    PO    R    H A2 -0.79667 1.73708  0.79667 2.71589 1.11104 0.08373
#> 4    PO    S    L A2 -0.61667 1.53333  0.61667 2.39732 0.98072 0.08647
#> 5    PO    R    H A1  0.01667 0.98851 -0.01667 1.54552 0.63225 0.08413
#> 6    PO    S    L A1  0.33000 0.79554 -0.33000 1.24380 0.50883 0.21284
#> 7    PO    S    M A2  0.33000 0.79554 -0.33000 1.24380 0.50883 0.25710
#> 8    PO    R    M A1  0.67333 0.62706 -0.67333 0.98039 0.40107 0.43880
#> 9    PO    S    M A1  1.27000 0.41466 -1.27000 0.64831 0.26522 0.25403
#> 10   PO    R    M A2  1.66667 0.31498 -1.66667 0.49246 0.20146 0.28898
#> 11   PO    R    L A1  1.81000 0.28519 -1.81000 0.44589 0.18241 0.02082
#> 12   PO    R    L A2  3.96333 0.06411 -3.96333 0.10023 0.04100 0.82277
#>    Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC sig
#> 1      4.74277     5.68675         2.17046         2.60246   a
#> 2      2.85784     3.08458         1.51119         1.63109  ab
#> 3      1.63913     1.84088         0.75175         0.84427  bc
#> 4      1.44412     1.62805         0.58079         0.65476   c
#> 5      0.93252     1.04787        -0.01767        -0.01572  cd
#> 6      0.68642     0.92200        -0.38246        -0.28474   d
#> 7      0.66568     0.95072        -0.39437        -0.27613   d
#> 8      0.46261     0.84996        -0.91269        -0.49675  de
#> 9      0.34771     0.49450        -1.51452        -1.06496  ef
#> 10     0.25780     0.38484        -2.03630        -1.36413   f
#> 11     0.28111     0.28934        -1.83631        -1.78407   f
#> 12     0.03624     0.11340        -7.01032        -2.24070   g
```
