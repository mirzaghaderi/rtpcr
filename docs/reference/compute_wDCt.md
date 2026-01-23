# Cleaning data and weighted delta Ct (wDCt) calculation

The `compute_wDCt` function cleans the data and computes wDCt. This
function is automatically applied to the expression analysis functions
like `ANOVA_DDCt`, `TTEST_DDCt`, etc. So it should not be applied in
advance of expression analysis functions.

## Usage

``` r
compute_wDCt(x, numOfFactors, numberOfrefGenes, block)
```

## Arguments

- x:

  A data frame containing experimental design columns, replicates
  (integer), target gene E/Ct column pairs, and reference gene E/Ct
  column pairs. Reference gene columns must be located at the end of the
  data frame.

- numOfFactors:

  Integer. Number of experimental factor columns (excluding `rep` and
  optional `block`).

- numberOfrefGenes:

  Integer. Number of reference genes.

- block:

  Character or `NULL`. Name of the blocking factor column. When a qPCR
  experiment is done in multiple qPCR plates, each plate is considered
  as a random block so that at least one replicate of each treatment and
  control is present on a plate.

## Value

The original data frame along with the weighted delta Ct column.

## Details

The `compute_wDCt` function computes weighted delta Ct (wDCt) for the
input data. Missing data can be denoted by NA in the input data frame.
Values such as '0' and 'undetermined' (for any E and Ct) are
automatically converted to NA. For target genes, NA for E or Ct
measurements cause returning NA for the corresponding delta Ct for that
replicate (row). If there are more than one reference gene, NA in the
place of the E or the Ct value cause skipping that gene and remaining
references are geometrically averaged. The `compute_wDCt` function is
automatically applied to the expression analysis functions.

## Examples

``` r
           
data <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
data
#>    Type Concentration  block Rep   PO        Ct_PO NLM       Ct_NLM ref1
#> 1     R            L1 plate1   1 1.88         33.3 1.9        28.77 1.75
#> 2     R            L1 plate1   2 1.88        33.39 1.9         28.2 1.75
#> 3     R            L1 plate2   3 1.88        33.34 1.9        28.48 1.75
#> 4     R            L1 plate2   4 1.88        33.34 1.9        28.64 1.75
#> 5     R            L2 plate1   1 1.88        32.73 1.9        31.17 1.75
#> 6     R            L2 plate1   2 1.88 undetermined 1.9        31.05 1.75
#> 7     R            L2 plate2   3 1.88         32.6 1.9        31.42 1.75
#> 8     R            L2 plate2   4 1.88        32.15 1.9         31.4 1.75
#> 9     R            L3 plate1   1 1.88        31.48 1.9        30.19 1.75
#> 10    R            L3 plate1   2 1.88        31.27 1.9 undetermined 1.75
#> 11    R            L3 plate2   3 1.88        31.32 1.9        30.22 1.75
#> 12    R            L3 plate2   4 1.88        31.32 1.9         30.4 1.75
#> 13    S            L1 plate1   1 1.88        32.85 1.9        32.26 1.75
#> 14    S            L1 plate1   2 1.88        32.17 1.9        32.85 1.75
#> 15    S            L1 plate2   3 1.88        32.99 1.9        32.55 1.75
#> 16    S            L1 plate2   4 1.88        32.89 1.9        32.75 1.75
#> 17    S            L2 plate1   1 1.88        32.41 1.9        32.97 1.75
#> 18    S            L2 plate1   2 1.88        32.49 1.9        32.83 1.75
#> 19    S            L2 plate2   3 1.88 undetermined 1.9         32.6 1.75
#> 20    S            L2 plate2   4 1.88         32.3 1.9        32.33 1.75
#> 21    S            L3 plate1   1 1.88        31.03 1.9         29.4 1.75
#> 22    S            L3 plate1   2 1.88        31.73 1.9        29.76 1.75
#> 23    S            L3 plate2   3 1.88        31.83 1.9        29.22 1.75
#> 24    S            L3 plate2   4 1.88        31.83 1.9        29.15 1.75
#>    Ct_ref1 ref2 Ct_ref2 ref3      Ct_ref3
#> 1    31.53    2   30.81    2        27.01
#> 2    31.57    2   30.71    2        27.17
#> 3    31.50    2   30.03    2        27.53
#> 4    30.50    2   30.03    2         27.9
#> 5    31.30    2   29.91    2         27.6
#> 6    32.55    2   30.05    2         27.3
#> 7    31.92    2   29.99    2        28.46
#> 8    30.81    2   30.55    2         28.4
#> 9    33.30    2   30.10    2         27.7
#> 10   33.37    2   30.17    2        28.14
#> 11   33.35    2   30.51    2        29.61
#> 12   29.35    2   30.51    2        29.51
#> 13   26.94    2   30.73    2         27.7
#> 14   27.69    2   30.68    2        27.68
#> 15   27.34    2   30.00    2 undetermined
#> 16   27.36    2   30.71    2        27.11
#> 17   28.70    2   30.58    2         28.7
#> 18   28.66    2   30.03    2        28.92
#> 19   28.71    2   31.06    2        28.43
#> 20   28.72    2   31.06    2        27.43
#> 21   30.61    2   31.14    2        27.42
#> 22   30.20    2   30.24    2        27.81
#> 23   30.49    2   30.11    2        28.29
#> 24   29.34    2   30.11    2        27.24
compute_wDCt(x = data,
             numOfFactors = 2,
             numberOfrefGenes = 3,
             block = "block")
#>    Type Concentration  block Rep   PO Ct_PO NLM Ct_NLM ref1 Ct_ref1 ref2
#> 1     R            L1 plate1   1 1.88 33.30 1.9  28.77 1.75   31.53    2
#> 2     R            L1 plate1   2 1.88 33.39 1.9  28.20 1.75   31.57    2
#> 3     R            L1 plate2   3 1.88 33.34 1.9  28.48 1.75   31.50    2
#> 4     R            L1 plate2   4 1.88 33.34 1.9  28.64 1.75   30.50    2
#> 5     R            L2 plate1   1 1.88 32.73 1.9  31.17 1.75   31.30    2
#> 6     R            L2 plate1   2   NA    NA 1.9  31.05 1.75   32.55    2
#> 7     R            L2 plate2   3 1.88 32.60 1.9  31.42 1.75   31.92    2
#> 8     R            L2 plate2   4 1.88 32.15 1.9  31.40 1.75   30.81    2
#> 9     R            L3 plate1   1 1.88 31.48 1.9  30.19 1.75   33.30    2
#> 10    R            L3 plate1   2 1.88 31.27  NA     NA 1.75   33.37    2
#> 11    R            L3 plate2   3 1.88 31.32 1.9  30.22 1.75   33.35    2
#> 12    R            L3 plate2   4 1.88 31.32 1.9  30.40 1.75   29.35    2
#> 13    S            L1 plate1   1 1.88 32.85 1.9  32.26 1.75   26.94    2
#> 14    S            L1 plate1   2 1.88 32.17 1.9  32.85 1.75   27.69    2
#> 15    S            L1 plate2   3 1.88 32.99 1.9  32.55 1.75   27.34    2
#> 16    S            L1 plate2   4 1.88 32.89 1.9  32.75 1.75   27.36    2
#> 17    S            L2 plate1   1 1.88 32.41 1.9  32.97 1.75   28.70    2
#> 18    S            L2 plate1   2 1.88 32.49 1.9  32.83 1.75   28.66    2
#> 19    S            L2 plate2   3   NA    NA 1.9  32.60 1.75   28.71    2
#> 20    S            L2 plate2   4 1.88 32.30 1.9  32.33 1.75   28.72    2
#> 21    S            L3 plate1   1 1.88 31.03 1.9  29.40 1.75   30.61    2
#> 22    S            L3 plate1   2 1.88 31.73 1.9  29.76 1.75   30.20    2
#> 23    S            L3 plate2   3 1.88 31.83 1.9  29.22 1.75   30.49    2
#> 24    S            L3 plate2   4 1.88 31.83 1.9  29.15 1.75   29.34    2
#>    Ct_ref2 ref3 Ct_ref3        wDCt
#> 1    30.81    2   27.01 -1.16612574
#> 2    30.71    2   27.17 -1.73033264
#> 3    30.03    2   27.53 -1.36500151
#> 4    30.03    2   27.90 -1.04254847
#> 5    29.91    2   27.60  1.19827424
#> 6    30.05    2   27.30  0.78209047
#> 7    29.99    2   28.46  0.93696586
#> 8    30.55    2   28.40  1.09624592
#> 9    30.10    2   27.70 -0.38014430
#> 10   30.17    2   28.14          NA
#> 11   30.51    2   29.61 -1.13474523
#> 12   30.51    2   29.51  0.27745548
#> 13   30.73    2   27.70  3.28655511
#> 14   30.68    2   27.68  3.60944806
#> 15   30.00   NA      NA  4.26074112
#> 16   30.71    2   27.11  3.79970050
#> 17   30.58    2   28.70  3.09910490
#> 18   30.03    2   28.92  3.07813001
#> 19   31.06    2   28.43  2.69725429
#> 20   31.06    2   27.43  2.77025159
#> 21   31.14    2   27.42 -0.54668179
#> 22   30.24    2   27.81  0.05099101
#> 23   30.11    2   28.29 -0.65484402
#> 24   30.11    2   27.24 -0.02400916
```
