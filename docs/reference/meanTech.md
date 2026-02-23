# Computing the mean of technical replicates

Computes the arithmetic mean of technical replicates for each sample or
group. This is often performed before ANOVA or other statistical
analyses to simplify comparisons between experimental groups.

## Usage

``` r
meanTech(
  x,
  groups,
  numOfFactors,
  numberOfrefGenes,
  block,
  set_missing_target_Ct_to_40 = FALSE
)
```

## Arguments

- x:

  A raw data frame containing technical replicates.

- groups:

  An integer vector or character vector specifying the column(s) to
  group by before calculating the mean of technical replicates.

- numOfFactors:

  Integer. Number of experimental factor columns

- numberOfrefGenes:

  Integer. Number of reference genes.

- block:

  Character. Block column name or `NULL`.

- set_missing_target_Ct_to_40:

  If `TRUE`, missing target gene Ct values become 40; if `FALSE`
  (default), they become NA.

## Value

A data frame with the mean of technical replicates for each group.

## Details

The `meanTech` function calculates the mean of technical replicates
based on one or more grouping columns. This reduces the dataset to a
single representative value per group, facilitating downstream analysis
such as fold change calculation or ANOVA.

## Author

Ghader Mirzaghaderi

## Examples

``` r
# Example input data frame with technical replicates
data1 <- read.csv(system.file("extdata", "data_withTechRep.csv", package = "rtpcr"))

# Calculate mean of technical replicates using first four columns as groups
meanTech(data1,
         groups = 1:2,
         numOfFactors = 1,
         numberOfrefGenes = 1,
         block = NULL)
#>   Condition biolrep target Ct_target ref   Ct_ref
#> 1      Ctrl       1      2  28.89400   2 30.29300
#> 2      Ctrl       2      2  32.33733   2 30.56600
#> 3      Ctrl       3      2  33.23000   2 31.59467
#> 4     Treat       1      2  32.98800   2 31.52700
#> 5     Treat       2      2  32.60000   2 31.77900
#> 6     Treat       3      2  30.68500   2 30.95233

# Another example using different dataset and grouping columns
data2 <- read.csv(system.file("extdata", "data_Lee_etal2020qPCR.csv", package = "rtpcr"))
meanTech(data2, groups = 1:3,
         numOfFactors = 2,
         numberOfrefGenes = 1,
         block = NULL)
#>    factor1  DS biolRep  APOE Ct_APOE GAPDH Ct_GAPDH
#> 1    DSWHi D12       1 2.120  24.710 2.190   14.675
#> 2    DSWHi D12       2 1.985  27.785 1.935   16.385
#> 3    DSWHi D12       3 2.030  26.395 2.025   14.290
#> 4    DSWHi D15       1 2.100  27.030 2.045   15.130
#> 5    DSWHi D15       2 2.005  29.780 2.075   15.870
#> 6    DSWHi D15       3 2.005  28.205 2.155   13.810
#> 7    DSWHi D18       1 2.110  28.865 2.060   15.095
#> 8    DSWHi D18       2 1.935  31.380 1.885   16.110
#> 9    DSWHi D18       3 2.105  30.385 2.075   14.005
#> 10   DSWHi  D7       1 2.080  24.175 2.100   15.260
#> 11   DSWHi  D7       2 1.970  26.480 2.025   16.715
#> 12   DSWHi  D7       3 2.195  24.830 2.175   14.725
#> 13    DSWi D12       1 1.945  25.330 2.080   15.840
#> 14    DSWi D12       2 1.920  24.160 2.030   14.725
#> 15    DSWi D12       3 2.105  23.830 2.335   14.745
#> 16    DSWi D15       1 1.905  28.615 1.895   16.180
#> 17    DSWi D15       2 1.935  28.280 1.990   14.755
#> 18    DSWi D15       3 2.125  27.005 2.395   14.735
#> 19    DSWi D18       1 1.930  29.730 2.040   16.585
#> 20    DSWi D18       2 2.085  28.655 1.970   14.490
#> 21    DSWi D18       3 2.090  28.680 2.140   14.060
#> 22    DSWi  D7       1 2.035  25.835 2.000   16.625
#> 23    DSWi  D7       2 2.085  24.145 2.080   14.895
#> 24    DSWi  D7       3 2.305  23.330 2.145   14.470
#> 25     DSi D12       1 1.970  26.460 1.960   16.315
#> 26     DSi D12       2 2.095  25.305 1.945   16.275
#> 27     DSi D12       3 1.965  24.650 2.065   14.835
#> 28     DSi D15       1 1.915  29.240 1.940   16.200
#> 29     DSi D15       2 1.975  27.275 1.975   16.390
#> 30     DSi D15       3 1.935  26.940 2.045   14.855
#> 31     DSi D18       1 2.005  30.980 1.965   16.230
#> 32     DSi D18       2 1.960  30.215 1.970   16.305
#> 33     DSi D18       3 2.110  29.890 2.075   14.735
#> 34     DSi  D7       1 2.000  25.090 2.020   17.065
#> 35     DSi  D7       2 2.045  24.260 1.915   16.280
#> 36     DSi  D7       3 1.910  23.680 2.035   15.205
```
