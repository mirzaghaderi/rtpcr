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
  analyseAllTarget = TRUE,
  model = NULL
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
  formula (CRD or RCBD based on `block` and `numOfFactors`). The formula
  uses `wDCt` as the response variable. For mixed models, random effects
  can be defined using `lmer` syntax (e.g.,
  `"wDCt ~ Treatment + (1|Block)"`). When using `model`, the `block` and
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
\`TTEST_DDCt()\`, \`WILCOX_DDCt()\`, \`ANOVA_DDCt()\`, and
\`ANOVA_DCt()\`) return the relative expression table which include fold
change and corresponding statistics. The output of \`ANOVA_DDCt()\`, and
\`ANOVA_DCt()\` also include lm models, residuals, raw data and ANOVA
table for each gene.

The expression table returned by \`TTEST_DDCt()\`, \`WILCOX_DDCt()\`,
and \`ANOVA_DDCt()\` functions include these columns: gene (name of
target genes), contrast (calibrator level and contrasts for which the
relative expression is computed), ddCt (mean of weighted delta delta Ct
values), RE (relative expression or fold change = 2^-ddCt), log2FC
(log(2) of relative expression or fold change), pvalue, sig (per-gene
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
#> 1   PO       L1  0.00000 1.00000  0.00000 0.00000 0.00000 0.21159     0.86358
#> 2   PO L2 vs L1 -1.06105 2.08644  1.06105 1.29842 3.35272 0.14499     1.88695
#> 3   PO L3 vs L1 -2.30692 4.94825  2.30692 3.18964 7.67647 0.29402     4.03592
#> 4  NLM       L1  0.00000 1.00000  0.00000 0.00000 0.00000 0.96436     0.51251
#> 5  NLM L2 vs L1  0.75074 0.59430 -0.75074 0.42164 0.83767 0.36616     0.46108
#> 6  NLM L3 vs L1 -1.56512 2.95901  1.56512 2.06844 4.23303 0.17132     2.62768
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1     1.15797         0.00000         0.00000 1.00000    
#> 2     2.30703         0.95960         1.17322 0.00110  **
#> 3     6.06681         1.88158         2.82840 0.00000 ***
#> 4     1.95119         0.00000         0.00000 1.00000    
#> 5     0.76601        -0.96764        -0.58245 0.00124  **
#> 6     3.33212         1.38987         1.76246 0.00000 ***
#> 
#> Note: Using default model for analysis: wDCt ~ block + Concentration * Type 
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
#>      gene             contrast     ddCt      RE   log2FC     LCL     UCL
#> 1 C2H2_26              control  0.00000 1.00000  0.00000 0.00000 0.00000
#> 2 C2H2_26 treatment vs control  1.19333 0.43729 -1.19333 0.19262 0.99273
#> 3 C2H2_01              control  0.00000 1.00000  0.00000 0.00000 0.00000
#> 4 C2H2_01 treatment vs control -2.00667 4.01853  2.00667 2.45984 6.56487
#> 5 C2H2_12              control  0.00000 1.00000  0.00000 0.00000 0.00000
#> 6 C2H2_12 treatment vs control -0.72000 1.64718  0.72000 0.95945 2.82787
#>        se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1 0.06009     0.95920     1.04253         0.00000         0.00000 1.00000    
#> 2 0.42176     0.32644     0.58578        -1.59854        -0.89084 0.04875   *
#> 3 0.13017     0.91372     1.09442         0.00000         0.00000 1.00000    
#> 4 0.21932     3.45180     4.67830         1.72367         2.33613 0.00141  **
#> 5 0.18502     0.87964     1.13683         0.00000         0.00000 1.00000    
#> 6 0.21127     1.42280     1.90695         0.62192         0.83355 0.06239   .
#> 
#> Note: Using default model for analysis: wDCt ~ Condition 
#> The control level was used as calibrator.
  
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
#> 1 Target              1  0.00000 1.00000  0.00000 0.0000  0.00000 1.40507
#> 2 Target time2 vs time1  0.22333 0.85658 -0.22333 0.0923  7.94918 0.97527
#> 3 Target time3 vs time1 -2.23333 4.70219  2.23333 0.5067 43.63677 0.55411
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1     0.37760     2.64830         0.00000         0.00000 1.00000    
#> 2     0.43570     1.68405        -0.43907        -0.11360 0.81660    
#> 3     3.20256     6.90403         1.52108         3.27911 0.06847   .
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
#> 1 Target              1  0.00000 1.00000 0.00000 0.00000 0.00000 0.70359
#> 2 Target time2 vs time1 -0.32833 1.25556 0.32833 0.43805 3.59873 0.53234
#> 3 Target time3 vs time1 -1.34667 2.54324 1.34667 0.88731 7.28951 0.70742
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
#> 1     0.61404     1.62855         0.00000         0.00000 1.00000    
#> 2     0.86814     1.81588         0.22702         0.47486 0.55402    
#> 3     1.55752     4.15279         0.82472         2.19894 0.03509   *
#> 
#> The 1 level was used as calibrator.
           
```
