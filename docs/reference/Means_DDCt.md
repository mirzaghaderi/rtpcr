# Pairwise comparisons of relative expression values (\\\Delta \Delta C_T\\) using a fitted model

Performs relative expression (fold change) analysis based on the
\\\Delta \Delta C_T\\ method using a fitted model object produced by
[`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
or
[`REPEATED_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/REPEATED_DDCt.md).

## Usage

``` r
Means_DDCt(model, specs, p.adj = "none")
```

## Arguments

- model:

  A fitted model object (typically an `lmer` or `lm` object) created by
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  or
  [`REPEATED_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/REPEATED_DDCt.md).

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
derived from a fitted model. For ANOVA models, FC values can be obtained
for main effects, interactions, and sliced (simple) effects. For ANCOVA
models returned by the rtpcr package, only simple effects are supported.

Internally, this function relies on the emmeans package to compute
marginal means and contrasts, which are then back-transformed to fold
change values using the \\\Delta \Delta C_T\\ framework.

## Author

Ghader Mirzaghaderi

## Examples

``` r
# Obtain a fitted model from ANOVA_DDCt
res <- ANOVA_DDCt(
  data_3factor,
  numOfFactors = 3,
  numberOfrefGenes = 1,
  mainFactor.column = 1,
  block = NULL)
#> NOTE: Results may be misleading due to involvement in interactions
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>              Df Sum Sq Mean Sq F value    Pr(>F)    
#> Type          1 24.834 24.8336 84.8207 2.392e-09 ***
#> Conc          2 45.454 22.7269 77.6252 3.319e-11 ***
#> SA            1  0.032  0.0324  0.1107 0.7422776    
#> Type:Conc     2 10.641  5.3203 18.1718 1.567e-05 ***
#> Type:SA       1  6.317  6.3168 21.5756 0.0001024 ***
#> Conc:SA       2  3.030  1.5150  5.1747 0.0135366 *  
#> Type:Conc:SA  2  3.694  1.8470  6.3086 0.0062852 ** 
#> Residuals    24  7.027  0.2928                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value    Pr(>F)    
#> SA         1  0.032  0.0324  0.0327    0.8577    
#> Conc       2 45.454 22.7269 22.9429 7.682e-07 ***
#> Type       1 24.834 24.8336 25.0696 2.105e-05 ***
#> Residuals 31 30.708  0.9906                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>   contrast     RE log2FC pvalue sig    LCL    UCL     se Lower.se.RE
#> 1        R 1.0000 0.0000      1     0.0000 0.0000 0.3939      0.7611
#> 2   S vs R 3.1626 1.6611      0 *** 2.4433 4.0936 0.3064      2.5575
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.3139          0.0000          0.0000
#> 2      3.9109          1.3433          2.0542
#> *** The R level was used as calibrator.
#> 
#> Relative Expression
#>   gene contrast     RE log2FC pvalue sig    LCL    UCL     se Lower.se.RE
#> 1 E_PO        R 1.0000 0.0000      1     0.0000 0.0000 0.3939      0.7611
#> 2 E_PO   S vs R 3.1626 1.6611      0 *** 2.4433 4.0936 0.3064      2.5575
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.3139          0.0000          0.0000
#> 2      3.9109          1.3433          2.0542

# Relative expression values for Type main effect
Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Type")
#>  contrast     RE        SE df      LCL      UCL p.value sig
#>  S vs R   3.1626 0.1803631 24 2.443349 4.093578 <0.0001 ***
#> 
#> Results are averaged over the levels of: Conc, SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration main effect
Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Conc")
#>  contrast        RE        SE df       LCL       UCL p.value sig
#>  L vs H   0.1703610 0.2208988 24 0.1242014 0.2336757 <0.0001 ***
#>  M vs H   0.2227247 0.2208988 24 0.1623772 0.3055004 <0.0001 ***
#>  M vs L   1.3073692 0.2208988 24 0.9531359 1.7932535  0.0928 .  
#> 
#> Results are averaged over the levels of: Type, SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration sliced by Type
Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Conc | Type")
#> Type = R:
#>  contrast       RE        SE df       LCL      UCL p.value sig
#>  L vs H   0.103187 0.3123981 24 0.0659984 0.161331 <0.0001 ***
#>  M vs H   0.339151 0.3123981 24 0.2169210 0.530255 <0.0001 ***
#>  M vs L   3.286761 0.3123981 24 2.1022126 5.138776 <0.0001 ***
#> 
#> Type = S:
#>  contrast       RE        SE df       LCL      UCL p.value sig
#>  L vs H   0.281265 0.3123981 24 0.1798969 0.439751 <0.0001 ***
#>  M vs H   0.146266 0.3123981 24 0.0935518 0.228684 <0.0001 ***
#>  M vs L   0.520030 0.3123981 24 0.3326112 0.813055  0.0059 ** 
#> 
#> Results are averaged over the levels of: SA 
#> Confidence level used: 0.95 

# Relative expression values for Concentration sliced by Type and SA
Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Conc | Type * SA")
#> Type = R, SA = A1:
#>  contrast       RE        SE df       LCL      UCL p.value sig
#>  L vs H   0.288505 0.4417977 24 0.1533437 0.542800  0.0005 ***
#>  M vs H   0.634342 0.4417977 24 0.3371606 1.193467  0.1502    
#>  M vs L   2.198724 0.4417977 24 1.1686485 4.136734  0.0167 *  
#> 
#> Type = S, SA = A1:
#>  contrast       RE        SE df       LCL      UCL p.value sig
#>  L vs H   0.267943 0.4417977 24 0.1424151 0.504115  0.0002 ***
#>  M vs H   0.139661 0.4417977 24 0.0742315 0.262761 <0.0001 ***
#>  M vs L   0.521233 0.4417977 24 0.2770416 0.980660  0.0438 *  
#> 
#> Type = R, SA = A2:
#>  contrast       RE        SE df       LCL      UCL p.value sig
#>  L vs H   0.036906 0.4417977 24 0.0196160 0.069436 <0.0001 ***
#>  M vs H   0.181327 0.4417977 24 0.0963776 0.341153 <0.0001 ***
#>  M vs L   4.913213 0.4417977 24 2.6114319 9.243840 <0.0001 ***
#> 
#> Type = S, SA = A2:
#>  contrast       RE        SE df       LCL      UCL p.value sig
#>  L vs H   0.295248 0.4417977 24 0.1569280 0.555487  0.0005 ***
#>  M vs H   0.153184 0.4417977 24 0.0814189 0.288203 <0.0001 ***
#>  M vs L   0.518830 0.4417977 24 0.2757643 0.976139  0.0425 *  
#> 
#> Confidence level used: 0.95 
```
