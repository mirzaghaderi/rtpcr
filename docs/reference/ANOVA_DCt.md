# \\\Delta C_T\\ ANOVA analysis

Performs \\\Delta C_T\\ analysis for target genes by applying \\\Delta
C_T\\ method to each target gene. Target genes must be provided as
paired efficiency (E) and Ct columns followed by the the reference
gene(s) columns. See "Input data structure and column arrangement" in
vignettes for details about data structure.

## Usage

``` r
ANOVA_DCt(
  x,
  numOfFactors,
  numberOfrefGenes,
  block,
  alpha = 0.05,
  p.adj = "none",
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

- block:

  Character or `NULL`. Name of the blocking factor column. When a qPCR
  experiment is done in multiple qPCR plates, variation resulting from
  the plates may interfere with the actual amount of gene expression.
  One solution is to conduct each plate as a randomized block so that at
  least one replicate of each treatment and control is present on a
  plate. Block effect is usually considered as random and its
  interaction with any main effect is not considered.

- alpha:

  statistical level for comparisons

- p.adj:

  Method for p-value adjustment. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- analyseAllTarget:

  Logical or character. If `TRUE` (default), all detected target genes
  are analysed. Alternatively, a character vector specifying the names
  (names of their Efficiency columns) of target genes to be analysed.

## Value

An object containing expression table, lm models, ANOVA tables,
residuals, raw data and ANOVA table for each gene.

- \\\Delta C_T\\ combined expression table:

  `object$combinedResults`

- ANOVA table for treatments:

  `object$perGene$gene_name$ANOVA_T`

- ANOVA table factorial:

  `object$perGene$gene_name$ANOVA_factorial`

- lm ANOVA for tratments:

  `object$perGene$gene_name$lm_T`

- lm ANOVA factorial:

  `object$perGene$gene_name$lm_factorial`

- Residuals:

  `resid(object$perGene$gene_name$lm_T)`

## Examples

``` r
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
```
