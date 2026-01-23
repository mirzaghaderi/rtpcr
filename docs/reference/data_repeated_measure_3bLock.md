# Sample data in (repeated-measures design with two reference genes)

A sample qPCR data set with repeated measurements over time, a blocking
factor, and two reference genes. The experiment includes two treatments
(untreated and treated) with repeated measures on the same individuals
(id).

## Usage

``` r
data_repeated_measure_3bLock
```

## Format

A data frame with 18 observations and 10 variables:

- treatment:

  Experimental treatment (untreated, treated)

- time:

  Time point of measurement

- blk:

  Blocking factor

- id:

  Subject/individual identifier for repeated measures

- Target:

  Mean amplification efficiency of target gene

- Ct_Target:

  Ct values of target gene

- Ref1:

  Mean amplification efficiency of reference gene 1

- Ct_Ref1:

  Ct values of reference gene 1

- Ref2:

  Mean amplification efficiency of reference gene 2

- Ct_Ref2:

  Ct values of reference gene 2

## Source

Not applicable

## Details

This dataset is suitable for repeated-measures mixed-model analysis and
normalization using multiple reference genes.
