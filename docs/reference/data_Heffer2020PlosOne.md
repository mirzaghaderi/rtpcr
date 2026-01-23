# Sample data in (multi-group design with 6 target genes)

A sample qPCR data set with four experimental groups (Uninjected,
SNConly, DMSORPE, SNCRPE) and six target genes normalized to one
reference gene (Gapdh). Each line belongs to a separate individual
(non-repeated measure experiment).

## Usage

``` r
data_Heffer2020PlosOne
```

## Format

A data frame with 18 observations and 16 variables:

- Treatment:

  Experimental treatment group

- rep:

  Biological replicates

- Fn1, Ct.Fn1:

  Efficiency and Ct values of Fn1 gene

- Col1a1, Ct.Col1a1:

  Efficiency and Ct values of Col1a1 gene

- Acta2, Ct.Acta2:

  Efficiency and Ct values of Acta2 gene

- TgfB, Ct.TgfB:

  Efficiency and Ct values of TgfB gene

- Tnfa, Ct.Tnfa:

  Efficiency and Ct values of Tnfa gene

- Mcp1, Ct.Mcp1:

  Efficiency and Ct values of Mcp1 gene

- Gapdh, Ct.Gapdh:

  Efficiency and Ct values of Gapdh reference gene

## Source

Not applicable

## Details

This dataset is suitable for t-test or one-way ANOVA based expression
analysis of multiple genes.
