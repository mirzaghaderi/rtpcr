# Delta Delta Ct ANOVA analysis

Apply Delta Delta Ct (ddCt) analysis to each target gene and performs
per-gene statistical analysis.

## Usage

``` r
ANOVA_DDCt(
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

  The input data frame containing experimental design columns,
  replicates (integer), target gene E/Ct column pairs, and reference
  gene E/Ct column pairs. Reference gene columns must be located at the
  right end of the data frame. See "Input data structure" in vignettes
  for details about data structure.

- numOfFactors:

  Integer. Number of experimental factor columns (excluding `rep` and
  optional `block`).

- numberOfrefGenes:

  Integer. Number of reference genes.

- mainFactor.column:

  Integer. Column index of the factor for which the relative expression
  analysis is applied.

- block:

  Character. Block column name or `NULL`. When a qPCR experiment is done
  in multiple qPCR plates, variation resulting from the plates may
  interfere with the actual amount of gene expression. One solution is
  to conduct each plate as a randomized block so that at least one
  replicate of each treatment and control is present on a plate. Block
  effect is usually considered as random and its interaction with any
  main effect is not considered.

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

ddCt analysis of variance (ANOVA) is performed for the
`mainFactor.column` based on a full model factorial experiment by
default. However, if `ANCOVA_DDCt` function is used, analysis of
covariance is performed for the levels of the `mainFactor.column` and
the other factors are treated as covariates. if the interaction between
the main factor and the covariate is significant, ANCOVA is not
appropriate.

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
data1 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
ANOVA_DDCt(x = data1,
           numOfFactors = 2,
           numberOfrefGenes = 2,
           block = "block",
           mainFactor.column = 2,
           p.adj = "none")
#> 
#> Relative Expression
#>   gene contrast     ddCt      RE   log2FC  pvalue sig     LCL     UCL      se
#> 1   PO       L1  0.00000 1.00000  0.00000 1.00000     0.00000 0.00000 0.22166
#> 2   PO L2 vs L1 -0.74029 1.67052  0.74029 0.01386   * 1.03525 2.69562 0.16264
#> 3   PO L3 vs L1 -1.60816 3.04864  1.60816 0.00001 *** 1.95755 4.74786 0.24084
#> 4  NLM       L1  0.00000 1.00000  0.00000 1.00000     0.00000 0.00000 0.64932
#> 5  NLM L2 vs L1  1.12438 0.45870 -1.12438 0.00001 *** 0.33080 0.63606 0.17305
#> 6  NLM L3 vs L1 -0.91097 1.88031  0.91097 0.00021 *** 1.33699 2.64440 0.15207
#> 7 ref1       L1  0.00000 1.00000  0.00000 1.00000     0.00000 0.00000 0.70454
#> 8 ref1 L2 vs L1  0.52026 0.69725 -0.52026 0.25068     0.32030 1.51779 0.58376
#> 9 ref1 L3 vs L1  1.38263 0.38352 -1.38263 0.00571  ** 0.17618 0.83486 0.52174
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1     0.85758     1.16608         0.00000         0.00000
#> 2     1.49242     1.86987         0.66137         0.82864
#> 3     2.57991     3.60252         1.36091         1.90034
#> 4     0.63758     1.56843         0.00000         0.00000
#> 5     0.40685     0.51716        -1.26767        -0.99728
#> 6     1.69219     2.08933         0.81983         1.01223
#> 7     0.61364     1.62963         0.00000         0.00000
#> 8     0.46522     1.04500        -0.77973        -0.34713
#> 9     0.26713     0.55061        -1.98501        -0.96304
#> The L1 level was used as calibrator.
           
data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
ANOVA_DDCt(x = data2,
           numOfFactors = 1,
           numberOfrefGenes = 1,
           block = NULL,
           mainFactor.column = 1,
           p.adj = "none")
#> 
#> Relative Expression
#>      gene             contrast     ddCt      RE   log2FC  pvalue sig     LCL
#> 1 C2H2_26              control  0.00000 1.00000  0.00000 1.00000     0.00000
#> 2 C2H2_26 treatment vs control  1.19333 0.43729 -1.19333 0.04875   * 0.19262
#> 3 C2H2_01              control  0.00000 1.00000  0.00000 1.00000     0.00000
#> 4 C2H2_01 treatment vs control -2.00667 4.01853  2.00667 0.00141  ** 2.45984
#> 5 C2H2_12              control  0.00000 1.00000  0.00000 1.00000     0.00000
#> 6 C2H2_12 treatment vs control -0.72000 1.64718  0.72000 0.06239   . 0.95945
#>       UCL      se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1 0.00000 0.06009     0.95920     1.04253         0.00000         0.00000
#> 2 0.99273 0.42176     0.32644     0.58578        -1.59854        -0.89084
#> 3 0.00000 0.13017     0.91372     1.09442         0.00000         0.00000
#> 4 6.56487 0.21932     3.45180     4.67830         1.72367         2.33613
#> 5 0.00000 0.18502     0.87964     1.13683         0.00000         0.00000
#> 6 2.82787 0.21127     1.42280     1.90695         0.62192         0.83355
#> The control level was used as calibrator.
           
```
