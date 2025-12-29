# ΔΔCt ANOVA analysis

Apply ΔΔCt analysis to each target gene in the input data frame. Target
and reference genes must be provided as paired efficiency (E) and Ct
columns located after the experimental design columns. columns.

## Usage

``` r
ANOVA_DDCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  mainFactor.column,
  analysisType = "anova",
  mainFactor.level.order = NULL,
  block = NULL,
  p.adj = "none",
  plot = FALSE,
  plotType = "RE",
  analyseAllTarget = TRUE
)
```

## Arguments

- x:

  A data frame containing experimental design columns, target gene E/Ct
  column pairs, and reference gene E/Ct column pairs. Reference gene
  columns must be located at the end of the data frame.

- numOfFactors:

  Integer. Number of experimental factor columns (excluding `rep` and
  optional `block`).

- numberOfrefGenes:

  Integer. Number of reference genes. Each reference gene must be
  represented by two columns (E and Ct).

- mainFactor.column:

  Column index or name of the factor for which relative expression is
  calculated. When `analysisType = "ancova"`, remaining factors are
  treated as covariates.

- analysisType:

  Character string specifying the analysis type; one of `"anova"`
  (default) or `"ancova"`.

- mainFactor.level.order:

  Optional character vector specifying the order of levels for the main
  factor. If `NULL`, the first observed level is used as the calibrator.
  If provided, the first element of the vector is used as the calibrator
  level.

- block:

  Character or `NULL`. Name of the blocking factor column.

- p.adj:

  Method for p-value adjustment. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- plot:

  Logical; if `FALSE`, per gene-plots are not generated.

- plotType:

  Plot scale to use: `"RE"` for relative expression or `"log2FC"` for
  log2 fold change.

- analyseAllTarget:

  Logical or character. If `TRUE` (default), all target genes are
  analysed. Alternatively, a character vector specifying the names
  (names of their Efficiency columns) of target genes to be analysed.

## Value

An object containing expression table, lm models, residuals, raw data
and ANOVA table for each gene.

- ΔΔCt combined expression table:

  `object$combinedFoldChange`

- ANOVA table:

  `object$perGene$gene_name$ANOVA_table`

- lm ANOVA:

  `object$perGene$gene_name$lm_ANOVA`

- lm ANCOVA:

  `object$perGene$gene_name$lm_ANCOVA`

- Residuals:

  `resid(object$perGene$gene_name$lm_ANOVA)`

- log2FC_Plot:

  `object$perGene$gene_name$log2FC_Plot`

- RE_Plot:

  `object$perGene$gene_name$RE_Plot`

## Details

ΔΔCt analysis is performed for the `mainFactor.column` based on a full
model factorial experiment by default. However, if `ancova`, the
`analysisType` argument, analysis of covariance is performed for the
levels of the `mainFactor.column` and the other factors are treated as
covariates. if the interaction between the main factor and the covariate
is significant, ANCOVA is not appropriate. ANCOVA is basically used when
a factor is affected by uncontrolled quantitative covariate(s). For
example, suppose that wDCt of a target gene in a plant is affected by
temperature. The gene may also be affected by drought. Since we already
know that temperature affects the target gene, we are interested to know
if the gene expression is also altered by the drought levels. We can
design an experiment to understand the gene behavior at both temperature
and drought levels at the same time. The drought is another factor (the
covariate) that may affect the expression of our gene under the levels
of the first factor i.e. temperature. The data of such an experiment can
be analyzed by ANCOVA or using ANOVA based on a factorial experiment.
ANCOVA is done even there is only one factor (without covariate or
factor variable).

## Examples

``` r
data1 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
ANOVA_DDCt(x = data1,
           numOfFactors = 2,
           numberOfrefGenes = 1,
           block = "block",
           mainFactor.column = 2,
           plot = FALSE,
           p.adj = "none")
#> NOTE: Results may be misleading due to involvement in interactions
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>                 Df Sum Sq Mean Sq F value    Pr(>F)    
#> block            1  3.956   3.956  2.1292 0.1665850    
#> factor2          2  0.866   0.433  0.2331 0.7950523    
#> factor1          1 50.978  50.978 27.4353 0.0001256 ***
#> factor2:factor1  2  1.901   0.950  0.5115 0.6103927    
#> Residuals       14 26.014   1.858                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value    Pr(>F)    
#> block      1  3.956   3.956  2.2677    0.1516    
#> factor1    1 51.704  51.704 29.6359 5.409e-05 ***
#> factor2    2  0.140   0.070  0.0402    0.9607    
#> Residuals 16 27.914   1.745                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>   contrast     RE log2FC pvalue sig    LCL    UCL     se Lower.se.RE
#> 1       L1 1.0000 0.0000 1.0000     0.0000 0.0000 0.9290      0.5252
#> 2 L2 vs L1 1.0601 0.0842 0.9138     0.2652 4.2370 0.5943      0.7022
#> 3 L3 vs L1 1.0919 0.1268 0.8610     0.3005 3.9677 0.6948      0.6746
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.9040          0.0000          0.0000
#> 2      1.6005          0.0558          0.1271
#> 3      1.7674          0.0784          0.2053
#> *** The L1 level was used as calibrator.
#> NOTE: Results may be misleading due to involvement in interactions
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>                 Df  Sum Sq Mean Sq F value    Pr(>F)    
#> block            1  4.3883  4.3883  3.5437 0.0780855 .  
#> factor2          2  3.9386  1.9693  1.5903 0.2344519    
#> factor1          1 25.5443 25.5443 20.6282 0.0003334 ***
#> factor2:factor1  2  6.1005  3.0503  2.4632 0.1167857    
#> Residuals       16 19.8130  1.2383                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df  Sum Sq Mean Sq F value    Pr(>F)    
#> block      1  4.3883  4.3883  3.0482 0.0978731 .  
#> factor1    1 24.6349 24.6349 17.1118 0.0006196 ***
#> factor2    2  4.8480  2.4240  1.6838 0.2136332    
#> Residuals 18 25.9136  1.4396                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>   contrast     RE  log2FC pvalue sig   LCL    UCL     se Lower.se.RE
#> 1       L1 1.0000  0.0000 1.0000     0.000 0.0000 0.6237      0.6490
#> 2 L2 vs L1 0.7935 -0.3338 0.5733     0.281 2.2407 0.6037      0.5221
#> 3 L3 vs L1 0.4520 -1.1457 0.0659   . 0.160 1.2763 0.5520      0.3083
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.5408          0.0000          0.0000
#> 2      1.2058         -0.5072         -0.2196
#> 3      0.6626         -1.6798         -0.7815
#> *** The L1 level was used as calibrator.
#> NOTE: Results may be misleading due to involvement in interactions
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>                 Df Sum Sq Mean Sq F value Pr(>F)
#> block            1  6.914  6.9138  2.2121 0.1564
#> factor2          2  9.373  4.6864  1.4994 0.2530
#> factor1          1  0.136  0.1356  0.0434 0.8377
#> factor2:factor1  2  4.004  2.0021  0.6406 0.5400
#> Residuals       16 50.007  3.1254               
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value Pr(>F)
#> block      1  6.914  6.9138  2.3041 0.1464
#> factor1    1  0.024  0.0245  0.0082 0.9290
#> factor2    2  9.484  4.7419  1.5803 0.2332
#> Residuals 18 54.011  3.0006               
#> 
#> Expression table
#>   contrast     RE log2FC pvalue sig    LCL     UCL     se Lower.se.RE
#> 1       L1 1.0000 0.0000 1.0000     0.0000  0.0000 1.0167      0.4942
#> 2 L2 vs L1 2.4017 1.2640 0.1894     0.4616 12.4961 0.2760      1.9834
#> 3 L3 vs L1 2.5785 1.3665 0.1578     0.4956 13.4163 0.3309      2.0500
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      2.0233          0.0000          0.0000
#> 2      2.9081          1.0439          1.5306
#> 3      3.2433          1.0864          1.7189
#> *** The L1 level was used as calibrator.
#> 
#> Relative Expression
#>     gene contrast     RE  log2FC pvalue sig    LCL     UCL     se Lower.se.RE
#> 1   E.PO       L1 1.0000  0.0000 1.0000     0.0000  0.0000 0.9290      0.5252
#> 2   E.PO L2 vs L1 1.0601  0.0842 0.9138     0.2652  4.2370 0.5943      0.7022
#> 3   E.PO L3 vs L1 1.0919  0.1268 0.8610     0.3005  3.9677 0.6948      0.6746
#> 4 E.PAO5       L1 1.0000  0.0000 1.0000     0.0000  0.0000 0.6237      0.6490
#> 5 E.PAO5 L2 vs L1 0.7935 -0.3338 0.5733     0.2810  2.2407 0.6037      0.5221
#> 6 E.PAO5 L3 vs L1 0.4520 -1.1457 0.0659   . 0.1600  1.2763 0.5520      0.3083
#> 7 E.ref1       L1 1.0000  0.0000 1.0000     0.0000  0.0000 1.0167      0.4942
#> 8 E.ref1 L2 vs L1 2.4017  1.2640 0.1894     0.4616 12.4961 0.2760      1.9834
#> 9 E.ref1 L3 vs L1 2.5785  1.3665 0.1578     0.4956 13.4163 0.3309      2.0500
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.9040          0.0000          0.0000
#> 2      1.6005          0.0558          0.1271
#> 3      1.7674          0.0784          0.2053
#> 4      1.5408          0.0000          0.0000
#> 5      1.2058         -0.5072         -0.2196
#> 6      0.6626         -1.6798         -0.7815
#> 7      2.0233          0.0000          0.0000
#> 8      2.9081          1.0439          1.5306
#> 9      3.2433          1.0864          1.7189
           
data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
ANOVA_DDCt(
           x = data2,
           numOfFactors = 1,
           numberOfrefGenes = 1,
           block = NULL,
           mainFactor.column = 1,
           plot = FALSE,
           p.adj = "none"
           )
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value  Pr(>F)  
#> Condition  1 2.1361 2.13607  7.8465 0.04875 *
#> Residuals  4 1.0889 0.27223                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value  Pr(>F)  
#> Condition  1 2.1361 2.13607  7.8465 0.04875 *
#> Residuals  4 1.0889 0.27223                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>               contrast     RE  log2FC pvalue sig    LCL    UCL     se
#> 1              control 1.0000  0.0000 1.0000     0.0000 0.0000 0.0601
#> 2 treatment vs control 0.4373 -1.1933 0.0488   * 0.1926 0.9927 0.4218
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      0.9592      1.0425          0.0000          0.0000
#> 2      0.3264      0.5858         -1.5985         -0.8908
#> *** The control level was used as calibrator.
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value  Pr(>F)   
#> Condition  1 6.0401  6.0401  61.907 0.00141 **
#> Residuals  4 0.3903  0.0976                   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value  Pr(>F)   
#> Condition  1 6.0401  6.0401  61.907 0.00141 **
#> Residuals  4 0.3903  0.0976                   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>               contrast     RE log2FC pvalue sig    LCL    UCL     se
#> 1              control 1.0000 0.0000 1.0000     0.0000 0.0000 0.1302
#> 2 treatment vs control 4.0185 2.0067 0.0014  ** 2.4598 6.5649 0.2193
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      0.9137      1.0944          0.0000          0.0000
#> 2      3.4518      4.6783          1.7237          2.3361
#> *** The control level was used as calibrator.
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value  Pr(>F)  
#> Condition  1 0.7776  0.7776  6.5731 0.06239 .
#> Residuals  4 0.4732  0.1183                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value  Pr(>F)  
#> Condition  1 0.7776  0.7776  6.5731 0.06239 .
#> Residuals  4 0.4732  0.1183                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>               contrast     RE log2FC pvalue sig    LCL    UCL     se
#> 1              control 1.0000   0.00 1.0000     0.0000 0.0000 0.1850
#> 2 treatment vs control 1.6472   0.72 0.0624   . 0.9595 2.8279 0.2113
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      0.8796      1.1368          0.0000          0.0000
#> 2      1.4228      1.9069          0.6219          0.8335
#> *** The control level was used as calibrator.
#> 
#> Relative Expression
#>      gene             contrast     RE  log2FC pvalue sig    LCL    UCL     se
#> 1 C2H2_26              control 1.0000  0.0000 1.0000     0.0000 0.0000 0.0601
#> 2 C2H2_26 treatment vs control 0.4373 -1.1933 0.0488   * 0.1926 0.9927 0.4218
#> 3 C2H2_01              control 1.0000  0.0000 1.0000     0.0000 0.0000 0.1302
#> 4 C2H2_01 treatment vs control 4.0185  2.0067 0.0014  ** 2.4598 6.5649 0.2193
#> 5 C2H2_12              control 1.0000  0.0000 1.0000     0.0000 0.0000 0.1850
#> 6 C2H2_12 treatment vs control 1.6472  0.7200 0.0624   . 0.9595 2.8279 0.2113
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      0.9592      1.0425          0.0000          0.0000
#> 2      0.3264      0.5858         -1.5985         -0.8908
#> 3      0.9137      1.0944          0.0000          0.0000
#> 4      3.4518      4.6783          1.7237          2.3361
#> 5      0.8796      1.1368          0.0000          0.0000
#> 6      1.4228      1.9069          0.6219          0.8335
```
