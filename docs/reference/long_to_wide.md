# Converts a 4-column qPCR long data format to wide format

Converts a 4-column (Condition, gene, Efficiency, Ct) qPCR long data
format to wide format

## Usage

``` r
long_to_wide(x)
```

## Arguments

- x:

  a 4-column (Condition, gene, Efficiency, Ct) qPCR long data

## Value

A wide qPCR data frame

## Details

Converts a 4-column (Condition, gene, Efficiency, Ct) qPCR long data
format to wide format

## Author

Ghader Mirzaghaderi

## Examples

``` r
df <- read.table(header = TRUE, text = "
Condition  Gene  E  Ct
control  C2H2-26  1.8  31.26
control  C2H2-26  1.8  31.01
control  C2H2-26  1.8  30.97
treatment  C2H2-26  1.8  32.65
treatment  C2H2-26  1.8  32.03
treatment  C2H2-26  1.8  32.4
control  C2H2-01  1.75  31.06
control  C2H2-01  1.75  30.41
control  C2H2-01  1.75  30.97
treatment  C2H2-01  1.75  28.85
treatment  C2H2-01  1.75  28.93
treatment  C2H2-01  1.75  28.9
control  C2H2-12  2  28.5
control  C2H2-12  2  28.4
control  C2H2-12  2  28.8
treatment  C2H2-12  2  27.9
treatment  C2H2-12  2  28
treatment  C2H2-12  2  27.9
control  ref  1.9  28.87
control  ref  1.9  28.42
control  ref  1.9  28.53
treatment  ref  1.9  28.31
treatment  ref  1.9  29.14
treatment  ref  1.9  28.63")

long_to_wide(df)
#>   Condition Rep E.C2H2-26 Ct.C2H2-26 E.C2H2-01 Ct.C2H2-01 E.C2H2-12 Ct.C2H2-12
#> 1   control   1       1.8      31.26      1.75      31.06         2       28.5
#> 2   control   2       1.8      31.01      1.75      30.41         2       28.4
#> 3   control   3       1.8      30.97      1.75      30.97         2       28.8
#> 4 treatment   1       1.8      32.65      1.75      28.85         2       27.9
#> 5 treatment   2       1.8      32.03      1.75      28.93         2       28.0
#> 6 treatment   3       1.8      32.40      1.75      28.90         2       27.9
#>   E.ref Ct.ref
#> 1   1.9  28.87
#> 2   1.9  28.42
#> 3   1.9  28.53
#> 4   1.9  28.31
#> 5   1.9  29.14
#> 6   1.9  28.63
```
