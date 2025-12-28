# Sample data in (two factor with blocking factor and 3 reference genes)

A sample qPCR data set with blocking factor and 3 reference genes. Each
line belongs to a separate individual (non-repeated measure experiment).

## Usage

``` r
data_2factorBlock3ref
```

## Format

A data frame with 18 observations and 8 variables:

- factor1:

  First experimental factor

- factor2:

  Second experimental factor

- block:

  blocking factor

- Rep:

  Biological replicates

- E_PO:

  Mean amplification efficiency of PO gene

- Ct_PO:

  Ct values of PO gene. Each is the mean of technical replicates

- E_GAPDH:

  Mean amplification efficiency of GAPDH gene

- Ct_GAPDH:

  Ct values of GAPDH gene. Each is the mean of technical replicates

- ref2E:

  Mean amplification efficiency of ref2 gene

- ref2Ct:

  Ct values of ref2 gene. Each is the mean of technical replicates

- ref3E:

  Mean amplification efficiency of ref3 gene

- ref3Ct:

  Ct values of GAPDH gene. Each is the mean of technical replicates

## Source

Not applicable
