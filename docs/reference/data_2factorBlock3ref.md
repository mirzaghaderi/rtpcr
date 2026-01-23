# Sample data in (two factor with blocking factor and 3 reference genes)

A sample qPCR data set with blocking factor and 3 reference genes. Each
line belongs to a separate individual (non-repeated measure experiment).

## Usage

``` r
data_2factorBlock3ref
```

## Format

A data frame with 18 observations and 8 variables:

- Type:

  First experimental factor

- Concentration:

  Second experimental factor

- block:

  blocking factor

- Rep:

  Biological replicates

- PO:

  Mean amplification efficiency of PO gene

- Ct_PO:

  Ct values of PO gene. Each is the mean of technical replicates

- NLM:

  Mean amplification efficiency of NLM gene

- Ct_NLM:

  Ct values of NLM gene. Each is the mean of technical replicates

- ref1:

  Mean amplification efficiency of ref1 gene

- Ct_ref1:

  Ct values of ref1 gene. Each is the mean of technical replicates

- ref2:

  Mean amplification efficiency of ref2 gene

- Ct_ref2:

  Ct values of ref2 gene. Each is the mean of technical replicates

- ref3:

  Mean amplification efficiency of ref3 gene

- Ct_ref3:

  Ct values of GAPDH gene. Each is the mean of technical replicates

## Source

Not applicable
