# Pairwise comparisons of relative expression values (Delta Ct or Delta Delta Ct) using a fitted model

Performs relative expression (fold change) analysis based on the
\\\Delta C_T\\ or \\\Delta \Delta C_T\\ methods using a fitted model
object produced by
[`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md),
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
  [`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md),
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
derived from a fitted model. For ANOVA models, relative expression
values can be obtained for main effects, interactions, and sliced
(simple) effects. For ANCOVA models returned by the rtpcr package, only
simple effects are supported.

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




data <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
res <- ANOVA_DCt(
  data,
  numOfFactors = 3,
  numberOfrefGenes = 1,
  block = NULL)
#> NULL
#> 
#> Relative expression (DCt method)
#>    Type Conc SA     RE  log2FC    LCL    UCL     se Lower.se.RE Upper.se.RE
#> 1     S    H A2 5.1934  2.3767 8.1197 3.3217 0.1309      4.7428      5.6867
#> 2     S    H A1 2.9690  1.5700 4.6420 1.8990 0.0551      2.8578      3.0846
#> 3     R    H A2 1.7371  0.7967 2.7159 1.1110 0.0837      1.6391      1.8409
#> 4     S    L A2 1.5333  0.6167 2.3973 0.9807 0.0865      1.4441      1.6280
#> 5     R    H A1 0.9885 -0.0167 1.5455 0.6323 0.0841      0.9325      1.0479
#> 6     S    L A1 0.7955 -0.3300 1.2438 0.5088 0.2128      0.6864      0.9220
#> 7     S    M A2 0.7955 -0.3300 1.2438 0.5088 0.2571      0.6657      0.9507
#> 8     R    M A1 0.6271 -0.6733 0.9804 0.4011 0.4388      0.4626      0.8500
#> 9     S    M A1 0.4147 -1.2700 0.6483 0.2652 0.2540      0.3477      0.4945
#> 10    R    M A2 0.3150 -1.6667 0.4925 0.2015 0.2890      0.2578      0.3848
#> 11    R    L A1 0.2852 -1.8100 0.4459 0.1824 0.0208      0.2811      0.2893
#> 12    R    L A2 0.0641 -3.9633 0.1002 0.0410 0.8228      0.0362      0.1134
#>    Lower.se.log2FC Upper.se.log2FC sig
#> 1           2.1705          2.6025   a
#> 2           1.5112          1.6311  ab
#> 3           0.7517          0.8443  bc
#> 4           0.5808          0.6548   c
#> 5          -0.0177         -0.0157  cd
#> 6          -0.3825         -0.2847   d
#> 7          -0.3944         -0.2761   d
#> 8          -0.9127         -0.4968  de
#> 9          -1.5145         -1.0650  ef
#> 10         -2.0363         -1.3641   f
#> 11         -1.8363         -1.7841   f
#> 12         -7.0103         -2.2407   g
#> 
#> Combined Expression Table (all genes)
#>    gene Type Conc SA     RE  log2FC    LCL    UCL     se Lower.se.RE
#> 1    PO    S    H A2 5.1934  2.3767 8.1197 3.3217 0.1309      4.7428
#> 2    PO    S    H A1 2.9690  1.5700 4.6420 1.8990 0.0551      2.8578
#> 3    PO    R    H A2 1.7371  0.7967 2.7159 1.1110 0.0837      1.6391
#> 4    PO    S    L A2 1.5333  0.6167 2.3973 0.9807 0.0865      1.4441
#> 5    PO    R    H A1 0.9885 -0.0167 1.5455 0.6323 0.0841      0.9325
#> 6    PO    S    L A1 0.7955 -0.3300 1.2438 0.5088 0.2128      0.6864
#> 7    PO    S    M A2 0.7955 -0.3300 1.2438 0.5088 0.2571      0.6657
#> 8    PO    R    M A1 0.6271 -0.6733 0.9804 0.4011 0.4388      0.4626
#> 9    PO    S    M A1 0.4147 -1.2700 0.6483 0.2652 0.2540      0.3477
#> 10   PO    R    M A2 0.3150 -1.6667 0.4925 0.2015 0.2890      0.2578
#> 11   PO    R    L A1 0.2852 -1.8100 0.4459 0.1824 0.0208      0.2811
#> 12   PO    R    L A2 0.0641 -3.9633 0.1002 0.0410 0.8228      0.0362
#>    Upper.se.RE Lower.se.log2FC Upper.se.log2FC sig
#> 1       5.6867          2.1705          2.6025   a
#> 2       3.0846          1.5112          1.6311  ab
#> 3       1.8409          0.7517          0.8443  bc
#> 4       1.6280          0.5808          0.6548   c
#> 5       1.0479         -0.0177         -0.0157  cd
#> 6       0.9220         -0.3825         -0.2847   d
#> 7       0.9507         -0.3944         -0.2761   d
#> 8       0.8500         -0.9127         -0.4968  de
#> 9       0.4945         -1.5145         -1.0650  ef
#> 10      0.3848         -2.0363         -1.3641   f
#> 11      0.2893         -1.8363         -1.7841   f
#> 12      0.1134         -7.0103         -2.2407   g

# lm <- res$perGene$PO$lm_factorial
# Means_DDCt(lm, specs = "Type * Conc * SA", p.adj = "none")
```
