# Delta Ct ANOVA analysis with optional model specification

Performs Delta Ct (dCt) analysis of the data from a factorial experiment
with support for both fixed- and mixed-effect models. Per-gene
statistical analysis and grouping is also performed.

## Usage

``` r
ANOVA_DCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  block = NULL,
  alpha = 0.05,
  p.adj = "none",
  analyseAllTarget = TRUE,
  model = NULL,
  modelBased_se = TRUE,
  set_missing_target_Ct_to_40 = FALSE
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
  main effect is not considered. Note: This parameter is ignored if
  `model` is provided.

- alpha:

  Statistical level for comparisons (default: 0.05).

- p.adj:

  Method for p-value adjustment. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- analyseAllTarget:

  Logical or character. If `TRUE` (default), all detected target genes
  are analysed. Alternatively, a character vector specifying the names
  (names of their Efficiency columns) of target genes to be analysed.

- model:

  Optional model formula. If provided, this overrides the automatic
  formula (CRD or RCBD based on `block` and `numOfFactors`). The formula
  uses `wDCt` as the response variable. For mixed models, random effects
  can be defined using `lmer` syntax (e.g.,
  `"wDCt ~ Treatment + (1|Block)"`). When using `model`, the `block` and
  `numOfFactors` arguments are ignored for model specification, but
  still used for data structure identification.

- modelBased_se:

  Logical. If `TRUE` (default), standard errors are calculated from
  model-based residuals. If `FALSE`, standard errors are calculated
  directly from the observed `wDCt` values within each treatment group
  according to the selected `se.type`. For single factor data, both
  methods are the same. It is recommended to use `modelBased_se = TRUE`
  (default).

- set_missing_target_Ct_to_40:

  If `TRUE`, missing target gene Ct values become 40; if `FALSE`
  (default), they become NA.

## Value

An object containing expression tables, lm/lmer models, ANOVA tables,
residuals, and raw data for each gene:

- `relativeExpression`:

  dCt expression table for all treatment combinations along with
  per-gene statistical grouping

- `perGene`:

  Nested list containing detailed results for each target gene:

  - `ANOVA_table`: Full factorial ANOVA table

  - `lm`: lm/lmer model for factorial design

  - `Final_data`: Processed data with wDCt values

  - `resid(object$perGene$gene_name$lm)`: Residuals

## Details

The function performs ANOVA analysis on weighted delta Ct (wDCt) values
and returns variance components along with an expression table
containing:

- `gene`: Name of target genes

- Factor columns: Experimental design factors

- `dCt`: Mean weighted delta Ct for each treatment combination

- `RE`: Relative expression = fold change = 2^-dCt

- `log2FC`: log2 of RE

- `LCL`: 95% lower confidence level

- `UCL`: 95% upper confidence level

- `se`: Standard error

- `Lower.se.RE`: Lower limit error bar for RE

- `Upper.se.RE`: Upper limit error bar for RE

- `Lower.se.log2FC`: Lower limit error bar for log2FC

- `Upper.se.log2FC`: Upper limit error bar for log2FC

- `sig`: Per-gene significance grouping letters

## Examples

``` r
# Default usage with fixed effects
result <- ANOVA_DCt(data_2factorBlock3ref, numOfFactors = 2, numberOfrefGenes = 3, 
                    block = "block")
#> 
#> Relative Expression
#> 
#>    gene Type Concentration      dCt      RE   log2FC     LCL     UCL      se
#> 1    PO    R            L1  2.62828 0.16174 -2.62828 0.11704 0.22351 0.03803
#> 2    PO    S            L1  3.12306 0.11478 -3.12306 0.08306 0.15862 0.19479
#> 3    PO    R            L2  1.65834 0.31680 -1.65834 0.22089 0.46897 0.29573
#> 4    PO    S            L2  2.20079 0.21752 -2.20079 0.14694 0.31197 0.05659
#> 5    PO    R            L3  0.08550 0.94246 -0.08550 0.68199 1.30241 0.27565
#> 6    PO    S            L3  1.28189 0.41126 -1.28189 0.29760 0.56833 0.28106
#> 7   NLM    R            L1 -1.32600 2.50707  1.32600 1.92666 3.26233 0.15235
#> 8   NLM    S            L1  3.50923 0.08782 -3.50923 0.06749 0.11428 0.12116
#> 9   NLM    R            L2  1.00339 0.49883 -1.00339 0.38334 0.64910 0.09175
#> 10  NLM    S            L2  2.91119 0.13294 -2.91119 0.10216 0.17298 0.09931
#> 11  NLM    R            L3 -0.41248 1.33097  0.41248 0.97812 1.80481 0.40782
#> 12  NLM    S            L3 -0.29364 1.22573  0.29364 0.94196 1.59498 0.17875
#>    Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC sig
#> 1      0.15753     0.16606        -2.66631        -2.59026  de
#> 2      0.10028     0.13137        -3.31785        -2.92826   e
#> 3      0.25809     0.38888        -1.95407        -1.36261  bc
#> 4      0.20915     0.22622        -2.25738        -2.14421  cd
#> 5      0.77854     1.14089        -0.36115         0.19016   a
#> 6      0.33846     0.49971        -1.56295        -1.00083   b
#> 7      2.25583     2.78629         1.17366         1.47835   a
#> 8      0.08075     0.09552        -3.63039        -3.38807   e
#> 9      0.46809     0.53158        -1.09514        -0.91165   c
#> 10     0.12409     0.14241        -3.01049        -2.81188   d
#> 11     1.00323     1.76577         0.00466         0.82030   b
#> 12     1.08289     1.38741         0.11488         0.47239   b
#> 
#> Note: Using default model for statistical analysis: wDCt ~ block + Type * Concentration 

# Mixed model with random block effect
result <- ANOVA_DCt(data_2factorBlock3ref, numOfFactors = 2, numberOfrefGenes = 3,
                          block = "block")
#> 
#> Relative Expression
#> 
#>    gene Type Concentration      dCt      RE   log2FC     LCL     UCL      se
#> 1    PO    R            L1  2.62828 0.16174 -2.62828 0.11704 0.22351 0.03803
#> 2    PO    S            L1  3.12306 0.11478 -3.12306 0.08306 0.15862 0.19479
#> 3    PO    R            L2  1.65834 0.31680 -1.65834 0.22089 0.46897 0.29573
#> 4    PO    S            L2  2.20079 0.21752 -2.20079 0.14694 0.31197 0.05659
#> 5    PO    R            L3  0.08550 0.94246 -0.08550 0.68199 1.30241 0.27565
#> 6    PO    S            L3  1.28189 0.41126 -1.28189 0.29760 0.56833 0.28106
#> 7   NLM    R            L1 -1.32600 2.50707  1.32600 1.92666 3.26233 0.15235
#> 8   NLM    S            L1  3.50923 0.08782 -3.50923 0.06749 0.11428 0.12116
#> 9   NLM    R            L2  1.00339 0.49883 -1.00339 0.38334 0.64910 0.09175
#> 10  NLM    S            L2  2.91119 0.13294 -2.91119 0.10216 0.17298 0.09931
#> 11  NLM    R            L3 -0.41248 1.33097  0.41248 0.97812 1.80481 0.40782
#> 12  NLM    S            L3 -0.29364 1.22573  0.29364 0.94196 1.59498 0.17875
#>    Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC sig
#> 1      0.15753     0.16606        -2.66631        -2.59026  de
#> 2      0.10028     0.13137        -3.31785        -2.92826   e
#> 3      0.25809     0.38888        -1.95407        -1.36261  bc
#> 4      0.20915     0.22622        -2.25738        -2.14421  cd
#> 5      0.77854     1.14089        -0.36115         0.19016   a
#> 6      0.33846     0.49971        -1.56295        -1.00083   b
#> 7      2.25583     2.78629         1.17366         1.47835   a
#> 8      0.08075     0.09552        -3.63039        -3.38807   e
#> 9      0.46809     0.53158        -1.09514        -0.91165   c
#> 10     0.12409     0.14241        -3.01049        -2.81188   d
#> 11     1.00323     1.76577         0.00466         0.82030   b
#> 12     1.08289     1.38741         0.11488         0.47239   b
#> 
#> Note: Using default model for statistical analysis: wDCt ~ block + Type * Concentration 

# Custom mixed model formula with nested random effects
result_custom <- ANOVA_DCt(data_repeated_measure_2, numOfFactors = 2, numberOfrefGenes = 1,
                            block = NULL,
                            model = wDCt ~ treatment * time + (1 | id))
#> Using user defined formula.
#> 
#> Relative Expression
#> 
#>     gene treatment time       dCt       RE   log2FC       LCL       UCL      se
#> 1 Target untreated   t1 -13.34333 10393.06 13.34333  2420.745  44620.87 0.35860
#> 2 Target   treated   t1 -14.08667 17398.40 14.08667  4052.422  74697.09 0.67044
#> 3 Target untreated   t2 -14.22333 19127.14 14.22333  4455.080  82119.16 0.37483
#> 4 Target   treated   t2 -13.86333 14903.19 13.86333  3471.240  63984.34 0.29080
#> 5 Target untreated   t3 -13.80333 14296.09 13.80333  3329.836  61377.88 0.27656
#> 6 Target   treated   t3 -16.32000 81810.59 16.32000 19055.267 351240.05 0.58482
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC sig
#> 1    8105.743    13325.83        12.98473        13.70194   b
#> 2   10931.679    27690.55        13.41623        14.75711   b
#> 3   14750.765    24801.93        13.84850        14.59816  ab
#> 4   12182.565    18231.38        13.57253        14.15414   b
#> 5   11802.211    17316.95        13.52677        14.07990  ab
#> 6   54545.868   122703.57        15.73518        16.90482   a

```
