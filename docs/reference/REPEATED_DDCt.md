# \\\Delta\Delta C_T\\ ANOVA analysis on repeated measure data

`REPEATED_DDCt` function performs \\\Delta \Delta C_T\\ method analysis
of observations repeatedly taken over different time courses. Data may
be obtained over time from a uni- or multi-factorial experiment. Target
genes must be provided as paired efficiency (E) and Ct columns followed
by the columns and the reference gene columns.

## Usage

``` r
REPEATED_DDCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  repeatedFactor,
  calibratorLevel,
  block = NULL,
  p.adj = "none",
  plot = FALSE,
  analyseAllTarget = TRUE
)
```

## Arguments

- x:

  input data frame in which the first column is id, followed by the
  factor column(s) which include at least time. The first level of time
  in data frame is used as calibrator or reference level. Additional
  factor(s) may also be present. Other columns are efficiency and Ct
  values of target and reference genes. **NOTE:** In the "id" column, a
  unique number is assigned to each individual from which samples have
  been taken over time, for example see `data_repeated_measure_1`, all
  the three number 1 indicate one individual which has been sampled over
  three different time courses. See example data sets or refer
  [`vignette`](https://mirzaghaderi.github.io/rtpcr/doc/vignette.md),
  section "data structure and column arrangement" for details.

- numOfFactors:

  Integer. Number of experimental factor columns (excluding optional
  `block`).

- numberOfrefGenes:

  Integer. Number of reference genes. Each reference gene must be
  represented by two columns (E and Ct).

- repeatedFactor:

  Character string specifying the factor for which fold changes are
  analysed (commonly `"time"`).

- calibratorLevel:

  A level of `repeatedFactor` to be used as the calibrator (reference
  level) which is the reference level or sample that all others are
  compared to. Examples are untreated or time 0.

- block:

  Optional blocking factor column name. If supplied, block effects are
  treated as random effects.

- p.adj:

  Method for p-value adjustment.

- plot:

  Logical; if `FALSE`, plots are not produced.

- analyseAllTarget:

  Logical or character. If `TRUE` (default), all detected target genes
  are analysed. Alternatively, a character vector specifying the names
  (names of their Efficiency columns) of target genes to be analysed.

## Value

An object containing expression table, lm models, residuals, raw data
and ANOVA table for each gene.

- \\\Delta\Delta C_T\\ combined expression table:

  `object$Relative_Expression_table`

- ANOVA table:

  `object$perGene$gene_name$ANOVA_table`

- lm ANOVA:

  `object$perGene$gene_name$lm`

- Residuals:

  `resid(object$perGene$gene_name$lm)`

- log2FC_Plot:

  `object$perGene$gene_name$log2FC_Plot`

- RE_Plot:

  `object$perGene$gene_name$RE_Plot`

## Details

Column layout requirements for `x`:

- Target gene columns: E/Ct pairs located between design and reference
  columns

- Reference gene columns: columns located at the end

## Examples

``` r
data <- read.csv(system.file("extdata", "data_repeated_measure_1.csv", package = "rtpcr"))
REPEATED_DDCt(
  data,
  numOfFactors = 1,
  numberOfrefGenes = 1,
  repeatedFactor = "time",
  calibratorLevel = "1",
  block = NULL)
#> Type III Analysis of Variance Table with Satterthwaite's method
#>       Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#> Time_ 11.073  5.5364     2     4  4.5382 0.09357 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Expression table
#>           contrast     RE  log2FC pvalue sig    LCL     UCL     se Lower.se.RE
#> 1           Time_1 1.0000  0.0000 1.0000     0.0000  0.0000 1.4051      0.3776
#> 2 Time_2 vs Time_1 0.8566 -0.2233 0.8166     0.0923  7.9492 0.9753      0.4357
#> 3 Time_3 vs Time_1 4.7022  2.2333 0.0685   . 0.5067 43.6368 0.5541      3.2026
#>   Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      2.6483          0.0000          0.0000
#> 2      1.6840         -0.4391         -0.1136
#> 3      6.9040          1.5211          3.2791
#> The level 1  of the selected factor was used as calibrator.
#> 
#> Combined Relative Expression Table (all genes)
#>     gene         contrast     RE  log2FC pvalue sig    LCL     UCL     se
#> 1 Target           Time_1 1.0000  0.0000 1.0000     0.0000  0.0000 1.4051
#> 2 Target Time_2 vs Time_1 0.8566 -0.2233 0.8166     0.0923  7.9492 0.9753
#> 3 Target Time_3 vs Time_1 4.7022  2.2333 0.0685   . 0.5067 43.6368 0.5541
#>   Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC
#> 1      0.3776      2.6483          0.0000          0.0000
#> 2      0.4357      1.6840         -0.4391         -0.1136
#> 3      3.2026      6.9040          1.5211          3.2791
```
