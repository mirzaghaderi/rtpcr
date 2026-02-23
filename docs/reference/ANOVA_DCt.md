# Delta Ct ANOVA analysis with optional model specification

Performs Delta Ct (dCt) analysis of the data from a 1-, 2-, or 3-factor
experiment with support for both fixed effects and mixed effects models.
Per-gene statistical grouping is performed for all treatment
combinations.

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

- `RE`: Relative expression = 2^-dCt

- `log2FC`: log2 of relative expression

- `LCL`: 95% lower confidence level

- `UCL`: 95% upper confidence level

- `se`: Standard error of the mean calculated from wDCt values

- `Lower.se.RE`: Lower limit error bar for RE (2^(log2(RE) - se))

- `Upper.se.RE`: Upper limit error bar for RE (2^(log2(RE) + se))

- `Lower.se.log2FC`: Lower limit error bar for log2 RE

- `Upper.se.log2FC`: Upper limit error bar for log2 RE

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
#> 1    PO    R            L1  2.62828 0.16174 -2.62828 0.11704 0.22351 0.06143
#> 2    PO    S            L1  3.12306 0.11478 -3.12306 0.08306 0.15862 0.21472
#> 3    PO    R            L2  1.65834 0.31680 -1.65834 0.22089 0.46897 0.25149
#> 4    PO    S            L2  2.20079 0.21752 -2.20079 0.14694 0.31197 0.05775
#> 5    PO    R            L3  0.08550 0.94246 -0.08550 0.68199 1.30241 0.26799
#> 6    PO    S            L3  1.28189 0.41126 -1.28189 0.29760 0.56833 0.30503
#> 7   NLM    R            L1 -1.32600 2.50707  1.32600 1.92666 3.26233 0.15025
#> 8   NLM    S            L1  3.50923 0.08782 -3.50923 0.06749 0.11428 0.11980
#> 9   NLM    R            L2  1.00339 0.49883 -1.00339 0.38334 0.64910 0.09128
#> 10  NLM    S            L2  2.91119 0.13294 -2.91119 0.10216 0.17298 0.10361
#> 11  NLM    R            L3 -0.41248 1.33097  0.41248 0.97812 1.80481 0.40799
#> 12  NLM    S            L3 -0.29364 1.22573  0.29364 0.94196 1.59498 0.17934
#>    Lower.se.RE Upper.se.RE Upper.se.log2FC Lower.se.log2FC sig
#> 1      0.15499     0.16877        -2.51872        -2.74261  de
#> 2      0.09891     0.13320        -2.69118        -3.62425   e
#> 3      0.26612     0.37713        -1.39305        -1.97415  bc
#> 4      0.20898     0.22640        -2.11443        -2.29068  cd
#> 5      0.78269     1.13484        -0.07100        -0.10295   a
#> 6      0.33288     0.50808        -1.03760        -1.58370   b
#> 7      2.25910     2.78226         1.47155         1.19485   a
#> 8      0.08083     0.09543        -3.22959        -3.81308   e
#> 9      0.46824     0.53141        -0.94187        -1.06893   c
#> 10     0.12372     0.14284        -2.70945        -3.12794   d
#> 11     1.00312     1.76598         0.54729         0.31087   b
#> 12     1.08244     1.38797         0.33250         0.25931   b
#> 
#> Note: Using default model for statistical analysis: wDCt ~ block + Type * Concentration 

# Mixed model with random block effect
result_mixed <- ANOVA_DCt(data_2factorBlock3ref, numOfFactors = 2, numberOfrefGenes = 3,
                          block = "block")
#> 
#> Relative Expression
#> 
#>    gene Type Concentration      dCt      RE   log2FC     LCL     UCL      se
#> 1    PO    R            L1  2.62828 0.16174 -2.62828 0.11704 0.22351 0.06143
#> 2    PO    S            L1  3.12306 0.11478 -3.12306 0.08306 0.15862 0.21472
#> 3    PO    R            L2  1.65834 0.31680 -1.65834 0.22089 0.46897 0.25149
#> 4    PO    S            L2  2.20079 0.21752 -2.20079 0.14694 0.31197 0.05775
#> 5    PO    R            L3  0.08550 0.94246 -0.08550 0.68199 1.30241 0.26799
#> 6    PO    S            L3  1.28189 0.41126 -1.28189 0.29760 0.56833 0.30503
#> 7   NLM    R            L1 -1.32600 2.50707  1.32600 1.92666 3.26233 0.15025
#> 8   NLM    S            L1  3.50923 0.08782 -3.50923 0.06749 0.11428 0.11980
#> 9   NLM    R            L2  1.00339 0.49883 -1.00339 0.38334 0.64910 0.09128
#> 10  NLM    S            L2  2.91119 0.13294 -2.91119 0.10216 0.17298 0.10361
#> 11  NLM    R            L3 -0.41248 1.33097  0.41248 0.97812 1.80481 0.40799
#> 12  NLM    S            L3 -0.29364 1.22573  0.29364 0.94196 1.59498 0.17934
#>    Lower.se.RE Upper.se.RE Upper.se.log2FC Lower.se.log2FC sig
#> 1      0.15499     0.16877        -2.51872        -2.74261  de
#> 2      0.09891     0.13320        -2.69118        -3.62425   e
#> 3      0.26612     0.37713        -1.39305        -1.97415  bc
#> 4      0.20898     0.22640        -2.11443        -2.29068  cd
#> 5      0.78269     1.13484        -0.07100        -0.10295   a
#> 6      0.33288     0.50808        -1.03760        -1.58370   b
#> 7      2.25910     2.78226         1.47155         1.19485   a
#> 8      0.08083     0.09543        -3.22959        -3.81308   e
#> 9      0.46824     0.53141        -0.94187        -1.06893   c
#> 10     0.12372     0.14284        -2.70945        -3.12794   d
#> 11     1.00312     1.76598         0.54729         0.31087   b
#> 12     1.08244     1.38797         0.33250         0.25931   b
#> 
#> Note: Using default model for statistical analysis: wDCt ~ block + Type * Concentration 

# Custom mixed model formula with nested random effects
result_custom <- ANOVA_DCt(data_repeated_measure_2, numOfFactors = 2, numberOfrefGenes = 1,
                            block = NULL,
                            model = wDCt ~ treatment * time + (1 | id))
#> Using user defined formula. Ignoring block and numOfFactors for model specification.
#> 
#> Relative Expression
#> 
#>     gene treatment time       dCt       RE   log2FC       LCL       UCL      se
#> 1 Target untreated   t1 -13.34333 10393.06 13.34333  2420.745  44620.87 0.60237
#> 2 Target   treated   t1 -14.08667 17398.40 14.08667  4052.422  74697.09 1.40507
#> 3 Target untreated   t2 -14.22333 19127.14 14.22333  4455.080  82119.16 0.65831
#> 4 Target   treated   t2 -13.86333 14903.19 13.86333  3471.240  63984.34 0.97527
#> 5 Target untreated   t3 -13.80333 14296.09 13.80333  3329.836  61377.88 0.78214
#> 6 Target   treated   t3 -16.32000 81810.59 16.32000 19055.267 351240.05 0.55411
#>   Lower.se.RE Upper.se.RE Upper.se.log2FC Lower.se.log2FC sig
#> 1    6845.630    15778.79        20.25790         8.78889   b
#> 2    6569.648    46076.16        37.30571         5.31914   b
#> 3   12119.302    30187.17        22.44780         9.01216  ab
#> 4    7580.446    29299.73        27.25537         7.05153   b
#> 5    8313.224    24584.72        23.73733         8.02668  ab
#> 6   55719.478   120119.09        23.96198        11.11521   a

```
