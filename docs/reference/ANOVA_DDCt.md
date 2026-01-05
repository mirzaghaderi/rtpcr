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
#>                 Df Sum Sq Mean Sq F value Pr(>F)
#> block            1  10.99  10.992  0.2980 0.5932
#> factor2          2  63.50  31.748  0.8608 0.4427
#> factor1          1   1.23   1.225  0.0332 0.8578
#> factor2:factor1  2  35.01  17.504  0.4746 0.6312
#> Residuals       15 553.24  36.883               
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value Pr(>F)
#> block      1  10.99  10.992  0.3177 0.5804
#> factor1    1   1.23   1.225  0.0354 0.8530
#> factor2    2  63.50  31.748  0.9175 0.4184
#> Residuals 17 588.25  34.603               
#> 
#> Expression table
#>   contrast      RE log2FC pvalue sig    LCL      UCL     se Lower.se.RE
#> 1       L1  1.0000 0.0000 1.0000     0.0000    0.000 3.1603      0.1119
#> 2 L2 vs L1 11.3694 3.5071 0.3018     0.0310 4171.614 0.5943      7.5307
#> 3 L3 vs L1 11.7105 3.5497 0.2606     0.0495 2772.692 0.6948      7.2347
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      8.9402           0.000          0.0000
#> 2     17.1649           2.323          5.2948
#> 3     18.9551           2.193          5.7458
#> *** The L1 level was used as calibrator.
#> NOTE: Results may be misleading due to involvement in interactions
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>                 Df Sum Sq Mean Sq F value Pr(>F)
#> block            1   8.84   8.838  0.2761 0.6060
#> factor2          2  40.41  20.203  0.6312 0.5440
#> factor1          1   0.07   0.071  0.0022 0.9631
#> factor2:factor1  2  46.33  23.165  0.7237 0.4993
#> Residuals       17 544.14  32.008               
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value Pr(>F)
#> block      1   8.84  8.8383  0.2844 0.6000
#> factor1    1   0.07  0.0705  0.0023 0.9625
#> factor2    2  40.41 20.2031  0.6501 0.5332
#> Residuals 19 590.47 31.0775               
#> 
#> Expression table
#>   contrast     RE log2FC pvalue sig    LCL       UCL     se Lower.se.RE
#> 1       L1 1.0000 0.0000 1.0000     0.0000    0.0000 3.1677      0.1113
#> 2 L2 vs L1 8.3810 3.0671 0.2934     0.0548 1281.8250 0.6037      5.5151
#> 3 L3 vs L1 4.7739 2.2552 0.4363     0.0312  730.1455 0.5520      3.2562
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      8.9858          0.0000          0.0000
#> 2     12.7360          2.0183          4.6609
#> 3      6.9992          1.5382          3.3064
#> *** The L1 level was used as calibrator.
#> NOTE: Results may be misleading due to involvement in interactions
#> ANOVA table 
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>                 Df Sum Sq Mean Sq F value Pr(>F)
#> block            1   7.89   7.889  0.2450 0.6269
#> factor2          2 114.25  57.124  1.7742 0.1996
#> factor1          1  31.92  31.924  0.9915 0.3333
#> factor2:factor1  2  34.84  17.418  0.5410 0.5919
#> Residuals       17 547.35  32.197               
#> 
#> ANCOVA table
#> Analysis of Variance Table
#> 
#> Response: wDCt
#>           Df Sum Sq Mean Sq F value Pr(>F)
#> block      1   7.89   7.889  0.2575 0.6177
#> factor1    1  31.92  31.924  1.0419 0.3202
#> factor2    2 114.25  57.124  1.8643 0.1823
#> Residuals 19 582.18  30.641               
#> 
#> Expression table
#>   contrast      RE log2FC pvalue sig    LCL      UCL     se Lower.se.RE
#> 1       L1  1.0000 0.0000 1.0000     0.0000    0.000 3.3047      0.1012
#> 2 L2 vs L1 23.8555 4.5763 0.1252     0.1537 3702.933 0.2760     19.7011
#> 3 L3 vs L1 25.6120 4.6788 0.1175     0.1650 3975.588 0.3309     20.3623
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      9.8816          0.0000          0.0000
#> 2     28.8860          3.7793          5.5413
#> 3     32.2152          3.7197          5.8850
#> *** The L1 level was used as calibrator.
#> 
#> Relative Expression
#>     gene contrast      RE log2FC pvalue sig    LCL       UCL     se Lower.se.RE
#> 1   E.PO       L1  1.0000 0.0000 1.0000     0.0000    0.0000 3.1603      0.1119
#> 2   E.PO L2 vs L1 11.3694 3.5071 0.3018     0.0310 4171.6138 0.5943      7.5307
#> 3   E.PO L3 vs L1 11.7105 3.5497 0.2606     0.0495 2772.6918 0.6948      7.2347
#> 4 E.PAO5       L1  1.0000 0.0000 1.0000     0.0000    0.0000 3.1677      0.1113
#> 5 E.PAO5 L2 vs L1  8.3810 3.0671 0.2934     0.0548 1281.8250 0.6037      5.5151
#> 6 E.PAO5 L3 vs L1  4.7739 2.2552 0.4363     0.0312  730.1455 0.5520      3.2562
#> 7 E.ref1       L1  1.0000 0.0000 1.0000     0.0000    0.0000 3.3047      0.1012
#> 8 E.ref1 L2 vs L1 23.8555 4.5763 0.1252     0.1537 3702.9329 0.2760     19.7011
#> 9 E.ref1 L3 vs L1 25.6120 4.6788 0.1175     0.1650 3975.5884 0.3309     20.3623
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      8.9402          0.0000          0.0000
#> 2     17.1649          2.3230          5.2948
#> 3     18.9551          2.1930          5.7458
#> 4      8.9858          0.0000          0.0000
#> 5     12.7360          2.0183          4.6609
#> 6      6.9992          1.5382          3.3064
#> 7      9.8816          0.0000          0.0000
#> 8     28.8860          3.7793          5.5413
#> 9     32.2152          3.7197          5.8850
           
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
