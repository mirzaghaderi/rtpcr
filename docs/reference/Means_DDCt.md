# Delta Delta Ct pairwise comparisons using a fitted model

Performs ddCt expression analysis using a fitted model object produced
by
[`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md)
or
[`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md).

## Usage

``` r
Means_DDCt(model, specs, p.adj = "none")
```

## Arguments

- model:

  A fitted model object (typically an `lmer` or `lm` object) created by
  [`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md),
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md).

- specs:

  A character string or character vector specifying the predictors or
  combinations of predictors over which relative expression values are
  desired. This argument follows the specification syntax used by
  [`emmeans::emmeans()`](https://rdrr.io/pkg/emmeans/man/emmeans.html)
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
  specs = "Type",
  block = NULL)
#> 
#> Relative Expression
#>   gene contrast     ddCt      RE     LCL     UCL  log2FC      se Lower.se.RE
#> 1   PO        R  0.00000 1.00000 0.00000 0.00000 0.00000 0.13774     0.90894
#> 2   PO   S vs R -1.99444 3.98463 3.07843 5.15758 1.99444 0.06318     3.81389
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC p.value sig
#> 1     1.10018        -0.13774         0.13774       1    
#> 2     4.16301         1.93126         2.05763       0 ***
#> 
#> The R level was used as calibrator.
#> Note: Using default model for statistical analysis: wDCt ~ Type * Conc * SA 

# Relative expression values for Type main effect
lm <- res$perGene$PO$lm
Means_DDCt(lm, specs = "Type")
#>  contrast       RE        SE df      LCL      UCL p.value sig
#>  S vs R   3.984626 0.1803631 24 3.078427 5.157585 <0.0001 ***
#> 
#> Results are averaged over the levels of: Conc, SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration main effect
Means_DDCt(lm, specs = "Conc")
#>  contrast       RE        SE df       LCL      UCL p.value sig
#>  M vs L   1.307369 0.2208988 24 0.9531359 1.793254  0.0928 .  
#>  H vs L   4.150639 0.2208988 24 3.0260179 5.693225 <0.0001 ***
#>  H vs M   3.174802 0.2208988 24 2.3145855 4.354718 <0.0001 ***
#> 
#> Results are averaged over the levels of: Type, SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration sliced by Type
Means_DDCt(lm, specs = "Conc | Type")
#> Type = R:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   3.286761 0.3123981 24 2.102213  5.138776 <0.0001 ***
#>  H vs L   4.845571 0.3123981 24 3.099227  7.575939 <0.0001 ***
#>  H vs M   1.474269 0.3123981 24 0.942943  2.304986  0.0857 .  
#> 
#> Type = S:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   0.520030 0.3123981 24 0.332611  0.813055  0.0059 ** 
#>  H vs L   3.555371 0.3123981 24 2.274015  5.558741 <0.0001 ***
#>  H vs M   6.836857 0.3123981 24 4.372854 10.689270 <0.0001 ***
#> 
#> Results are averaged over the levels of: SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration sliced by Type and SA
Means_DDCt(lm, specs = "Conc | Type * SA")
#> Type = R, SA = A1:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   2.198724 0.4417977 24 1.168649  4.136734  0.0167 *  
#>  H vs L   3.466148 0.4417977 24 1.842300  6.521297  0.0005 ***
#>  H vs M   1.576436 0.4417977 24 0.837895  2.965946  0.1502    
#> 
#> Type = S, SA = A1:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   0.521233 0.4417977 24 0.277042  0.980660  0.0438 *  
#>  H vs L   3.732132 0.4417977 24 1.983673  7.021725  0.0002 ***
#>  H vs M   7.160201 0.4417977 24 3.805733 13.471379 <0.0001 ***
#> 
#> Type = R, SA = A2:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   4.913213 0.4417977 24 2.611432  9.243840 <0.0001 ***
#>  H vs L   6.773962 0.4417977 24 3.600443 12.744701 <0.0001 ***
#>  H vs M   1.378724 0.4417977 24 0.732808  2.593965  0.3047    
#> 
#> Type = S, SA = A2:
#>  contrast       RE        SE df      LCL       UCL p.value sig
#>  M vs L   0.518830 0.4417977 24 0.275764  0.976139  0.0425 *  
#>  H vs L   3.386981 0.4417977 24 1.800221  6.372350  0.0005 ***
#>  H vs M   6.528116 0.4417977 24 3.469773 12.282159 <0.0001 ***
#> 
#> Confidence level used: 0.95 
```
