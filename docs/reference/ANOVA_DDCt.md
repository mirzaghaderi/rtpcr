# Delta Delta Ct ANOVA analysis with optional model specification

Apply Delta Delta Ct (ddCt) analysis to each target gene and performs
per-gene statistical analysis.

## Usage

``` r
ANOVA_DDCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  specs,
  block,
  calibratorLevel = NULL,
  p.adj = "none",
  analyseAllTarget = TRUE,
  model = NULL,
  set_missing_target_Ct_to_40 = FALSE,
  se.type = c("single.group", "paired.group", "two.group"),
  modelBased_se = TRUE,
  ...
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

- specs:

  Example: "A", "A\|B" or A\|B\*C" if A, B and C are name of factor
  columns in the input data The first name (here A) is the factor for
  which the relative expression analysis is applied.

- block:

  Character. Block column name or `NULL`. When a qPCR experiment is done
  in multiple qPCR plates, variation resulting from the plates may
  interfere with the actual amount of gene expression. One solution is
  to conduct each plate as a randomized block so that at least one
  replicate of each treatment and control is present on a plate. Block
  effect is usually considered as random and its interaction with any
  main effect is not considered.

- calibratorLevel:

  NULL or one of the levels of the first selected factor in specs
  argument. If NULL the first level of that factor is used as
  calibrator. Optional character vector specifying the order of levels
  for the main factor. If `NULL`, the first observed level is used as
  the calibrator. If provided, the first element of the vector is used
  as the calibrator level.

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
  methods are the same. It is recommended to use `modelBased_se = TRUE`
  (default).

- ...:

  Additional arguments. Included for backward compatibility with
  deprecated `mainFactor.column`.

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
ANOVA_DDCt(data1,
           numOfFactors = 2,
           numberOfrefGenes = 3,
           block = "block",
           specs = "Concentration",
           p.adj = "none")
#> 
#> Relative Expression
#>   gene contrast     ddCt      RE     LCL     UCL   log2FC      se Lower.se.RE
#> 1   PO       L1  0.00000 1.00000 0.00000 0.00000  0.00000 0.09187     0.93830
#> 2   PO L2 vs L1 -0.94610 1.92666 1.35851 2.73242 -0.94610 0.13465     1.75497
#> 3   PO L3 vs L1 -2.19198 4.56931 3.30647 6.31447 -2.19198 0.18224     4.02710
#> 4  NLM       L1  0.00000 1.00000 0.00000 0.00000  0.00000 0.09011     0.93945
#> 5  NLM L2 vs L1  0.86568 0.54879 0.42174 0.71411  0.86568 0.06259     0.52549
#> 6  NLM L3 vs L1 -1.44341 2.71964 2.06638 3.57940 -1.44341 0.18135     2.39838
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC p.value sig
#> 1     1.06575        -0.09187         0.09187 1.00000    
#> 2     2.11515         0.81145         1.08076 0.00116  **
#> 3     5.18453         2.00974         2.37421 0.00000 ***
#> 4     1.06445        -0.09011         0.09011 1.00000    
#> 5     0.57312        -0.92826        -0.80309 0.00018 ***
#> 6     3.08392         1.26206         1.62477 0.00000 ***
#> 
#> The L1 level was used as calibrator.
#> Note: Using default model for statistical analysis: wDCt ~ block + Concentration * Type 
           
data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
ANOVA_DDCt(data2,
           numOfFactors = 1,
           numberOfrefGenes = 1,
           block = NULL,
           specs = "Condition",
           p.adj = "none",
           se.type = "single.group")
#> 
#> Relative Expression
#>      gene             contrast     ddCt      RE     LCL     UCL   log2FC
#> 1 C2H2_26              control  0.00000 1.00000 0.00000 0.00000  0.00000
#> 2 C2H2_26 treatment vs control  1.19333 0.43729 0.19262 0.99273  1.19333
#> 3 C2H2_01              control  0.00000 1.00000 0.00000 0.00000  0.00000
#> 4 C2H2_01 treatment vs control -2.00667 4.01853 2.45984 6.56487 -2.00667
#> 5 C2H2_12              control  0.00000 1.00000 0.00000 0.00000  0.00000
#> 6 C2H2_12 treatment vs control -0.72000 1.64718 0.95945 2.82787 -0.72000
#>        se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC p.value sig
#> 1 0.06009     0.95920     1.04253        -0.06009         0.06009 1.00000    
#> 2 0.42176     0.32644     0.58578        -1.61509        -0.77158 0.04875   *
#> 3 0.13017     0.91372     1.09442        -0.13017         0.13017 1.00000    
#> 4 0.21932     3.45180     4.67830         1.78735         2.22598 0.00141  **
#> 5 0.18502     0.87964     1.13683        -0.18502         0.18502 1.00000    
#> 6 0.21127     1.42280     1.90695         0.50873         0.93127 0.06239   .
#> 
#> The control level was used as calibrator.
#> Note: Using default model for statistical analysis: wDCt ~ Condition 
  
# Repeated measure analysis
a <- ANOVA_DDCt(data_repeated_measure_1,
           numOfFactors = 1,
           numberOfrefGenes = 1,
           block = NULL,
           specs = "time",
           p.adj = "none", model = wDCt ~ time + (1 | id))
#> 
#> Relative Expression
#>     gene contrast     ddCt      RE     LCL      UCL   log2FC      se
#> 1 Target       t1  0.00000 1.00000 0.00000  0.00000  0.00000 0.67643
#> 2 Target t2 vs t1  0.22333 0.85658 0.15102  4.85870  0.22333 0.29515
#> 3 Target t3 vs t1 -2.23333 4.70219 0.82899 26.67166 -2.23333 0.58051
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC p.value sig
#> 1     0.62571     1.59818        -0.67643         0.67643 1.00000    
#> 2     0.69810     1.05104        -0.51849         0.07182 0.81660    
#> 3     3.14447     7.03157         1.65282         2.81385 0.06847   .
#> 
#> The t1 level was used as calibrator.

a$perGene$Target$ANOVA_table
#> Type III Analysis of Variance Table with Satterthwaite's method
#>      Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#> time 11.073  5.5364     2     4  4.5382 0.09357 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# Repeated measure analysis: split-plot in time
a <- ANOVA_DDCt(data_repeated_measure_2,
           numOfFactors = 2, numberOfrefGenes = 1,
           specs = "time", block = NULL,
           model = wDCt ~ treatment * time + (1 | id))
#> 
#> Relative Expression
#>     gene contrast     ddCt      RE     LCL     UCL   log2FC      se Lower.se.RE
#> 1 Target       t1  0.00000 1.00000 0.00000 0.00000  0.00000 0.34003     0.79003
#> 2 Target t2 vs t1 -0.32833 1.25556 0.53676 2.93695 -0.32833 0.21216     1.08385
#> 3 Target t3 vs t1 -1.34667 2.54324 1.08725 5.94901 -1.34667 0.28931     2.08112
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC p.value sig
#> 1     1.26578        -0.34003         0.34003 1.00000    
#> 2     1.45447         0.11617         0.54050 0.55402    
#> 3     3.10798         1.05736         1.63598 0.03509   *
#> 
#> The t1 level was used as calibrator.
   
```
