# Sample data in (one factor with one reference gene)

A sample qPCR data set with one experimental factor (condition) and one
reference gene. Each line belongs to a separate individual (non-repeated
measure experiment).

## Usage

``` r
data_Yuan2006PMCBioinf
```

## Format

A data frame with 24 observations and 6 variables:

- condition:

  Experimental factor with two levels (control, treatment)

- rep:

  Biological replicates

- target:

  Mean amplification efficiency of target gene

- Ct_target:

  Ct values of target gene. Each is the mean of technical replicates

- ref:

  Mean amplification efficiency of reference gene

- Ct_ref:

  Ct values of reference gene. Each is the mean of technical replicates

## Source

Yuan2006PMCBioinf
