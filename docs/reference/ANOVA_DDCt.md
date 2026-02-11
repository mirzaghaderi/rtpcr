# Delta Delta Ct ANOVA analysis with optional model specification

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
  analyseAllTarget = TRUE,
  model = NULL,
  set_missing_target_Ct_to_40 = FALSE,
  se.type = c("paired.group", "two.group", "single.group"),
  modelBased_se = TRUE
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

- model:

  Optional model formula. If provided, this overrides the automatic
  formula (uni - or multi-factorial CRD or RCBD based on `block` and
  `numOfFactors`). The formula uses `wDCt` as the response variable. For
  mixed models, random effects can be defined using `lmer` syntax (e.g.,
  `"wDCt ~ Treatment + (1 | id)"`). When using `model`, the `block` and
  `numOfFactors` arguments are ignored for model specification, but
  still used for data structure identification.

  for fixed effects only, the `"lm"` (ordinary least squares) is used.
  `"lmer"` is used for mixed effects models (requires the `lmerTest`
  package). If a custom formula is provided with random effects, the
  function will use
  [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html);
  otherwise it will use
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html). Note that `emmeans`
  supports both model types and will use appropriate degrees of freedom
  methods (Satterthwaite by default).

- set_missing_target_Ct_to_40:

  If `TRUE`, missing target gene Ct values become 40; if `FALSE`
  (default), they become NA.

- se.type:

  Character string specifying how standard error is calculated. One of
  `"paired.group"`, `"two.group"`, or `"single.group"`. `"paired.group"`
  computes SE from paired differences (used when a random `id` effect is
  present), `"two.group"` uses the unpaired two-group t-test standard
  error against the reference level, and `"single.group"` computes SE
  within each level using a one-group t-test.

- modelBased_se:

  Logical. If `TRUE` (default), standard errors are calculated from
  model-based residuals. If `FALSE`, standard errors are calculated
  directly from the observed `wDCt` values within each treatment group
  according to the selected `se.type`. For single factor data, both
  methods are the same. It is recommended to use modelBased_se = TRUE
  (default).

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
[`TTEST_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/TTEST_DDCt.md),
[`WILCOX_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/WILCOX_DDCt.md),
`ANOVA_DDCt()`, and
[`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md))
return the relative expression table which include fold change and
corresponding statistics. The output of `ANOVA_DDCt()`, and
[`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md)
also include lm models, residuals, raw data and ANOVA table for each
gene.

The expression table returned by
[`TTEST_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/TTEST_DDCt.md),
[`WILCOX_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/WILCOX_DDCt.md),
and `ANOVA_DDCt()` functions include these columns: gene (name of target
genes), contrast (calibrator level and contrasts for which the relative
expression is computed), ddCt (mean of weighted delta delta Ct values),
RE (relative expression or fold change = 2^-ddCt), log2FC (log(2) of
relative expression or fold change), pvalue, sig (per-gene
significance), LCL (95% lower confidence level), UCL (95% upper
confidence level), se (standard error of mean calculated from the
weighted delta Ct values of each of the main factor levels), Lower.se.RE
(The lower limit error bar for RE which is 2^(log2(RE) - se)),
Upper.se.RE (The upper limit error bar for RE which is 2^(log2(RE) +
se)), Lower.se.log2FC (The lower limit error bar for log2 RE), and
Upper.se.log2FC (The upper limit error bar for log2 RE)

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
           numberOfrefGenes = 3,
           block = "block",
           mainFactor.column = 2,
           p.adj = "none")
#> 
#> Relative Expression
#>   gene contrast     ddCt      RE   log2FC     LCL     UCL      se Lower.se.RE
#> 1   PO       L1  0.00000 1.00000  0.00000 0.00000 0.00000 0.00000     1.00000
#> 2   PO L2 vs L1 -0.94610 1.92666  0.94610 1.25860 2.94934 0.16301     1.72082
#> 3   PO L3 vs L1 -2.19198 4.56931  2.19198 3.08069 6.77724 0.20409     3.96657
#> 4  NLM       L1  0.00000 1.00000  0.00000 0.00000 0.00000 0.00000     1.00000
#> 5  NLM L2 vs L1  0.86568 0.54879 -0.86568 0.39830 0.75614 0.10971     0.50860
#> 6  NLM L3 vs L1 -1.44341 2.71964  1.44341 1.94670 3.79946 0.20250     2.36347
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1     1.00000         0.00000         0.00000 1.00000    
#> 2     2.15713         0.84502         1.05928 0.00116  **
#> 3     5.26364         1.90283         2.52506 0.00000 ***
#> 4     1.00000         0.00000         0.00000 1.00000    
#> 5     0.59215        -0.93407        -0.80229 0.00018 ***
#> 6     3.12947         1.25439         1.66093 0.00000 ***
#> 
#> The L1 level was used as calibrator.
#> Note: Using default model for statistical analysis: wDCt ~ block + Concentration * Type 
           
data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
ANOVA_DDCt(x = data2,
           numOfFactors = 1,
           numberOfrefGenes = 1,
           block = NULL,
           mainFactor.column = 1,
           p.adj = "none")
#> 
#> Relative Expression
#>      gene             contrast     ddCt      RE   log2FC     LCL     UCL
#> 1 C2H2_26              control  0.00000 1.00000  0.00000 0.00000 0.00000
#> 2 C2H2_26 treatment vs control  1.19333 0.43729 -1.19333 0.19262 0.99273
#> 3 C2H2_01              control  0.00000 1.00000  0.00000 0.00000 0.00000
#> 4 C2H2_01 treatment vs control -2.00667 4.01853  2.00667 2.45984 6.56487
#> 5 C2H2_12              control  0.00000 1.00000  0.00000 0.00000 0.00000
#> 6 C2H2_12 treatment vs control -0.72000 1.64718  0.72000 0.95945 2.82787
#>        se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1 0.00000     1.00000     1.00000         0.00000         0.00000 1.00000    
#> 2 0.42602     0.32548     0.58751        -1.60326        -0.88822 0.04875   *
#> 3 0.00000     1.00000     1.00000         0.00000         0.00000 1.00000    
#> 4 0.25504     3.36738     4.79558         1.68152         2.39469 0.00141  **
#> 5 0.00000     1.00000     1.00000         0.00000         0.00000 1.00000    
#> 6 0.28083     1.35582     2.00115         0.59264         0.87472 0.06239   .
#> 
#> The control level was used as calibrator.
#> Note: Using default model for statistical analysis: wDCt ~ Condition 
  
# Repeated measure analysis         
a <- ANOVA_DDCt(data_repeated_measure_1,
           numOfFactors = 1,
           numberOfrefGenes = 1,
           block = NULL,
           mainFactor.column = 1,
           p.adj = "none", model = wDCt ~ time + (1 | id))
#> 
#> Relative Expression
#>     gene       contrast     ddCt      RE   log2FC    LCL      UCL      se
#> 1 Target              1  0.00000 1.00000  0.00000 0.0000  0.00000 0.00000
#> 2 Target time2 vs time1  0.22333 0.85658 -0.22333 0.0923  7.94918 0.43682
#> 3 Target time3 vs time1 -2.23333 4.70219  2.23333 0.5067 43.63677 1.21754
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1     1.00000     1.00000         0.00000         0.00000 1.00000    
#> 2     0.63281     1.15949        -0.30231        -0.16499 0.81660    
#> 3     2.02201    10.93496         0.96037         5.19362 0.06847   .
#> 
#> The 1 level was used as calibrator.

a$perGene$Target$ANOVA_table
#> Type III Analysis of Variance Table with Satterthwaite's method
#>      Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#> time 11.073  5.5364     2     4  4.5382 0.09357 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# Repeated measure analysis: split-plot in time
a <- ANOVA_DDCt(data_repeated_measure_2,
           numOfFactors = 2, numberOfrefGenes = 1,
           mainFactor.column = 2, block = NULL,
           model = wDCt ~ treatment * time + (1 | id))
#> 
#> Relative Expression
#>     gene       contrast     ddCt      RE  log2FC     LCL     UCL      se
#> 1 Target              1  0.00000 1.00000 0.00000 0.00000 0.00000 0.00000
#> 2 Target time2 vs time1 -0.32833 1.25556 0.32833 0.43805 3.59873 0.38172
#> 3 Target time3 vs time1 -1.34667 2.54324 1.34667 0.88731 7.28951 0.58112
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1     1.00000     1.00000         0.00000         0.00000 1.00000    
#> 2     0.96367     1.63586         0.25200         0.42778 0.55402    
#> 3     1.70001     3.80471         0.90017         2.01463 0.03509   *
#> 
#> The 1 level was used as calibrator.
           
```
