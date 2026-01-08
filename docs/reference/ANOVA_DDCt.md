# Delta Delta Ct ANOVA analysis

Apply \\\Delta \Delta C_T\\ analysis to each target gene in the input
data frame. Target and reference genes must be provided as paired
efficiency (E) and Ct columns located after the experimental design
columns. columns.

## Usage

``` r
ANOVA_DDCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  mainFactor.column,
  block,
  analysisType = "anova",
  mainFactor.level.order = NULL,
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

- block:

  Character or `NULL`. Name of the blocking factor column. When a qPCR
  experiment is done in multiple qPCR plates, variation resulting from
  the plates may interfere with the actual amount of gene expression.
  One solution is to conduct each plate as a randomized block so that at
  least one replicate of each treatment and control is present on a
  plate. Block effect is usually considered as random and its
  interaction with any main effect is not considered.

- analysisType:

  Character string specifying the analysis type; one of `"anova"`
  (default) or `"ancova"`.

- mainFactor.level.order:

  Optional character vector specifying the order of levels for the main
  factor. If `NULL`, the first observed level is used as the calibrator.
  If provided, the first element of the vector is used as the calibrator
  level.

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

- \\\Delta \Delta C_T\\ combined expression table:

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

\\\Delta \Delta C_T\\ analysis is performed for the `mainFactor.column`
based on a full model factorial experiment by default. However, if
`ancova`, the `analysisType` argument, analysis of covariance is
performed for the levels of the `mainFactor.column` and the other
factors are treated as covariates. if the interaction between the main
factor and the covariate is significant, ANCOVA is not appropriate.
ANCOVA is basically used when a factor is affected by uncontrolled
quantitative covariate(s). For example, suppose that wDCt of a target
gene in a plant is affected by temperature. The gene may also be
affected by drought. Since we already know that temperature affects the
target gene, we are interested to know if the gene expression is also
altered by the drought levels. We can design an experiment to understand
the gene behavior at both temperature and drought levels at the same
time. The drought is another factor (the covariate) that may affect the
expression of our gene under the levels of the first factor i.e.
temperature. The data of such an experiment can be analyzed by ANCOVA or
using ANOVA based on a factorial experiment. ANCOVA is done even there
is only one factor (without covariate or factor variable).

## Examples

``` r
data1 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
ANOVA_DDCt(x = data1,
           numOfFactors = 2,
           numberOfrefGenes = 2,
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
#> block            1  1.836   1.836  1.2206    0.2866    
#> factor2          2  1.725   0.862  0.5732    0.5756    
#> factor1          1 58.636  58.636 38.9781 1.574e-05 ***
#> factor2:factor1  2  0.566   0.283  0.1882    0.8304    
#> Residuals       15 22.565   1.504                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value    Pr(>F)    
#> block      1  1.836   1.836  1.3495    0.2614    
#> factor1    1 58.636  58.636 43.0936 4.816e-06 ***
#> factor2    2  1.725   0.862  0.6338    0.5427    
#> Residuals 17 23.131   1.361                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>   contrast     RE  log2FC pvalue sig   LCL    UCL     se Lower.se.RE
#> 1       L1 1.0000  0.0000 1.0000     0.000 0.0000 0.8318      0.5618
#> 2 L2 vs L1 0.6525 -0.6159 0.3672     0.198 2.1504 0.6218      0.4240
#> 3 L3 vs L1 0.6818 -0.5525 0.3818     0.226 2.0568 0.7245      0.4126
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.7799          0.0000          0.0000
#> 2      1.0041         -0.9478         -0.4003
#> 3      1.1266         -0.9130         -0.3344
#> *** The L1 level was used as calibrator.
#> NOTE: Results may be misleading due to involvement in interactions
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>                 Df  Sum Sq Mean Sq F value    Pr(>F)    
#> block            1  2.2118  2.2118  2.1255  0.163088    
#> factor2          2 13.6670  6.8335  6.5669  0.007707 ** 
#> factor1          1 30.9361 30.9361 29.7291 4.301e-05 ***
#> factor2:factor1  2  3.1292  1.5646  1.5035  0.250470    
#> Residuals       17 17.6902  1.0406                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df  Sum Sq Mean Sq F value    Pr(>F)    
#> block      1  2.2118  2.2118  2.0185  0.171592    
#> factor1    1 30.9361 30.9361 28.2326 3.967e-05 ***
#> factor2    2 13.6670  6.8335  6.2363  0.008275 ** 
#> Residuals 19 20.8194  1.0958                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>   contrast     RE  log2FC pvalue sig    LCL    UCL     se Lower.se.RE
#> 1       L1 1.0000  0.0000 1.0000     0.0000 0.0000 0.5922      0.6633
#> 2 L2 vs L1 0.5053 -0.9847 0.0704   . 0.2040 1.2516 0.5838      0.3372
#> 3 L3 vs L1 0.2780 -1.8471 0.0021  ** 0.1122 0.6884 0.5217      0.1936
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.5076          0.0000          0.0000
#> 2      0.7573         -1.4759         -0.6570
#> 3      0.3991         -2.6518         -1.2866
#> *** The L1 level was used as calibrator.
#> 
#> Relative Expression
#>     gene contrast     RE  log2FC pvalue sig    LCL    UCL     se Lower.se.RE
#> 1   E.PO       L1 1.0000  0.0000 1.0000     0.0000 0.0000 0.8318      0.5618
#> 2   E.PO L2 vs L1 0.6525 -0.6159 0.3672     0.1980 2.1504 0.6218      0.4240
#> 3   E.PO L3 vs L1 0.6818 -0.5525 0.3818     0.2260 2.0568 0.7245      0.4126
#> 4 E.PAO5       L1 1.0000  0.0000 1.0000     0.0000 0.0000 0.5922      0.6633
#> 5 E.PAO5 L2 vs L1 0.5053 -0.9847 0.0704   . 0.2040 1.2516 0.5838      0.3372
#> 6 E.PAO5 L3 vs L1 0.2780 -1.8471 0.0021  ** 0.1122 0.6884 0.5217      0.1936
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      1.7799          0.0000          0.0000
#> 2      1.0041         -0.9478         -0.4003
#> 3      1.1266         -0.9130         -0.3344
#> 4      1.5076          0.0000          0.0000
#> 5      0.7573         -1.4759         -0.6570
#> 6      0.3991         -2.6518         -1.2866
           
data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
ANOVA_DDCt(
           x = data2,
           numOfFactors = 1,
           numberOfrefGenes = 1,
           block = NULL,
           mainFactor.column = 1,
           plot = FALSE,
           p.adj = "none")
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
