# Delta Delta Ct pairwise comparisons using a fitted model

Performs relative expression (fold change) analysis based on the Delta
Delta Ct (ddCt) methods using a fitted model object produced by
[`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md),
[`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
or `REPEATED_DDCt()`.

## Usage

``` r
Means_DDCt(model, specs, p.adj = "none")
```

## Arguments

- model:

  A fitted model object (typically an `lmer` or `lm` object) created by
  [`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md),
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  or `REPEATED_DDCt()`.

- specs:

  A character string or character vector specifying the predictors or
  combinations of predictors over which relative expression values are
  desired. This argument follows the specification syntax used by
  [`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
  (e.g., `"Factor"`, `"Factor1 | Factor2"`).

- p.adj:

  Character string specifying the method for adjusting p-values. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for available
  options.

## Value

A data frame containing estimated relative expression values, confidence
intervals, p-values, and significance levels derived from the fitted
model.

## Details

The `Means_DDCt` function performs pairwise comparisons of relative
expression values fo all combinations using estimated marginal means
derived from a fitted model. For ANOVA models, relative expression
values can be obtained for main effects, interactions, and sliced
(simple) effects. For ANCOVA models returned by the rtpcr package, only
simple effects are supported.

Internally, this function relies on the emmeans package to compute
marginal means and contrasts, which are then back-transformed to fold
change values using the ddCt framework.

## Author

Ghader Mirzaghaderi

## Examples

``` r
data <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))

# Obtain a fitted model from ANOVA_DDCt
res <- ANOVA_DDCt(
  data,
  numOfFactors = 3,
  numberOfrefGenes = 1,
  mainFactor.column = 1,
  block = NULL)
#> 
#> Relative Expression
#>   gene contrast     ddCt     RE  log2FC     LCL     UCL      se Lower.se.RE
#> 1   PO        R  0.00000 1.0000 0.00000 0.00000 0.00000 0.39385     0.76109
#> 2   PO   S vs R -1.66111 3.1626 1.66111 2.44335 4.09358 0.30640     2.55746
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC pvalue sig
#> 1     1.31390         0.00000         0.00000      1    
#> 2     3.91093         1.34327         2.05416      0    
#> 
#> The R level was used as calibrator.
#> Note: Using default model for statistical analysis: wDCt ~ Type * Conc * SA 

# Relative expression values for Type main effect
lm <- res$perGene$PO$lm
Means_DDCt(lm, specs = "Type")
#>  contrast     RE        SE df      LCL      UCL p.value sig
#>  S vs R   3.1626 0.1803631 24 2.443349 4.093578 <0.0001    
#> 
#> Results are averaged over the levels of: Conc, SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration main effect
Means_DDCt(lm, specs = "Conc")
#>  contrast       RE        SE df      LCL      UCL p.value sig
#>  M vs L   1.307369 0.2208988 24 0.953136 1.793254  0.0928 .  
#>  H vs L   5.869889 0.2208988 24 4.279436 8.051436 <0.0001    
#>  H vs M   4.489848 0.2208988 24 3.273318 6.158502 <0.0001    
#> 
#> Results are averaged over the levels of: Type, SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration sliced by Type
Means_DDCt(lm, specs = "Conc | Type")
#> Type = R:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   3.286761 0.3123981 24 2.102213  5.138776 <0.0001    
#>  H vs L   9.691142 0.3123981 24 6.198455 15.151878 <0.0001    
#>  H vs M   2.948538 0.3123981 24 1.885885  4.609972 <0.0001    
#> 
#> Type = S:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   0.520030 0.3123981 24 0.332611  0.813055  0.0059 ** 
#>  H vs L   3.555371 0.3123981 24 2.274015  5.558741 <0.0001    
#>  H vs M   6.836857 0.3123981 24 4.372854 10.689270 <0.0001    
#> 
#> Results are averaged over the levels of: SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration sliced by Type and SA
Means_DDCt(lm, specs = "Conc | Type * SA")
#> Type = R, SA = A1:
#>  contrast        RE        SE df       LCL      UCL p.value sig
#>  M vs L    2.198724 0.4417977 24  1.168649  4.13673  0.0167 *  
#>  H vs L    3.466148 0.4417977 24  1.842300  6.52130  0.0005    
#>  H vs M    1.576436 0.4417977 24  0.837895  2.96595  0.1502    
#> 
#> Type = S, SA = A1:
#>  contrast        RE        SE df       LCL      UCL p.value sig
#>  M vs L    0.521233 0.4417977 24  0.277042  0.98066  0.0438 *  
#>  H vs L    3.732132 0.4417977 24  1.983673  7.02173  0.0002    
#>  H vs M    7.160201 0.4417977 24  3.805733 13.47138 <0.0001    
#> 
#> Type = R, SA = A2:
#>  contrast        RE        SE df       LCL      UCL p.value sig
#>  M vs L    4.913213 0.4417977 24  2.611432  9.24384 <0.0001    
#>  H vs L   27.095850 0.4417977 24 14.401772 50.97880 <0.0001    
#>  H vs M    5.514895 0.4417977 24  2.931233 10.37586 <0.0001    
#> 
#> Type = S, SA = A2:
#>  contrast        RE        SE df       LCL      UCL p.value sig
#>  M vs L    0.518830 0.4417977 24  0.275764  0.97614  0.0425 *  
#>  H vs L    3.386981 0.4417977 24  1.800221  6.37235  0.0005    
#>  H vs M    6.528116 0.4417977 24  3.469773 12.28216 <0.0001    
#> 
#> Confidence level used: 0.95 
```
