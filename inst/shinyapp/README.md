This is the shiny version of the rtpcr package, a web application developed using R/Shiny for interactive analysis of qPCR data

rtpcr uses dCt and ddCt methods, including t-tests and ANOVA,
repeated-measures models, and publication-ready visualizations. The
package implements a general calculation method adopted from Ganger et
al. (2017) and Taylor et al. (2019), covering both the Livak and Pfaffl
methods. The rtpcr package gets efficiency (E) and Ct values of genes and
performs different expression analysis.


# Input data structure

For relative expression analysis (using `TTEST_DDCt`, `WILCOX_DDCt`,
`ANOVA_DCt`, and `ANOVA_DDCt` functions), the amplification efficiency
(`E`) and `Ct` or `Cq` values (the mean of technical replicates) is used
for the input table. If the `E` values are not available you should use
‘2’ instead representing the complete primer amplification efficiency.
The input data table should include the following columns from left to
wright:

1.  Experimental condition columns (and one block if available [NOTE
    1](#note-1))
2.  Replicates information (biological replicates or subjects; see [NOTE
    2](#note-2), and [NOTE 3](#note-3))  
3.  Target genes efficiency and Ct values (a pair column for each gene).
4.  Reference genes efficiency and Ct values (a pair column for each
    gene) [NOTE 4](#note-4).

The package supports **one or more target or reference gene(s)**,
supplied as efficiency–Ct column pairs. Reference gene columns must
always appear last. Two sample input data sets are presented below.

<figure>
<img src="readme/sampleData1.png" class="center"
style="width:100.0%"
alt="Figure 1: A sample input data with one experimetal factor, replicate column and E/Ct information of target and reference genes" />
<figcaption aria-hidden="true">Figure 1: A sample input data with one
experimetal factor, replicate column and E/Ct information of target and
reference genes</figcaption>
</figure>

<br>

If there is no blocking factor, the block column should be omitted.
However, a column for biological replicates (which may be named “Rep”,
“id” or similar) is always required.

<br>

<figure>
<img src="readme/dataStructure1.png" class="center"
style="width:100.0%"
alt="Figure 2: A sample input data with two experimetal factors, blocking factor, replicate column and E/Ct information of target and reference genes" />
<figcaption aria-hidden="true">Figure 2: A sample input data with two
experimetal factors, blocking factor, replicate column and E/Ct
information of target and reference genes</figcaption>
</figure>

#### NOTE 1

When a qPCR experiment is done in multiple qPCR plates, variation
resulting from the plates may interfere with the actual amount of gene
expression. One solution is to conduct each plate as a randomized block
so that at least one replicate of each treatment and control is present
on a plate. Block effect is usually considered as random and its
interaction with any main effect is not considered.

#### NOTE 2

For `TTEST_DDCt` and `WILCOX_DDCt` (independent groups), `ANOVA_DCt`,
and `ANOVA_DDCt` each row is from a separate and unique biological
replicate. For example, a data frame with 12 rows has come from an
experiment with 12 individuals. The repeated measure models are intended
for experiments with repeated observations (e.g. time-course data). In
repeated measure experiments the Replicate column contains identifiers
for each individual (id or subject). For example, all rows with a `1` at
Rep column correspond to a single individual, all rows with a `2`
correspond to another individual, and so on, which have been sampled at
specific time points.

#### NOTE 3

Your data table may also include a column of technical replicates (For
example, using one target and one reference genes, if you want to have 4
biological and 3 technical replicates under Control and Treatment
conditions, there would be a table of 24 rows containing both biological
replicates and technical replicate columns in the data). In this case,
the `meanTech` function should be applied first to calculate the mean of
the technical replicates. The resulting collapsed table is then used as
the input for expression analysis. To use the `meanTech` function
correctly, the technical replicate column must appear immediately after
the biological replicate column (see [Mean of technical
replicates](#mean-of-technical-replicates) for an example).

#### NOTE 4

Complete amplification efficiency (E) in the rtpcr package is denoted by
2. This means that 2 indicates 100%, and 1.85 and 1.70 indicate 0.85%
and 0.70% amplification efficiencies.

# Handling missing data

Missing Ct values for target genes is Handled using the
`set_missing_target_Ct_to_40` function. If `TRUE`, missing target gene
Ct values become 40; if `FALSE` (default), they become NA. missing Ct
values of reference genes are always converted to NA. If there are more
than one reference gene, NA in the place of the E or the Ct value of
cause skipping that gene and remaining references are geometrically
averaged in that replicate.



### Example of data format for amplification Efficiency

The `efficiency` function calculates the amplification efficiency (E),
slope, and R² statistics for genes, and performs pairwise comparisons of
slopes. It takes a data frame in which the first column contains the
dilution ratios, followed by the Ct value columns for each gene.

``` r
# dilutions Gene1   Gene2   Gene3
# 1.00  25.58   24.25   22.61
# 1.00  25.54   24.13   22.68
# 1.00  25.50   24.04   22.63
# 0.50  26.71   25.56   23.67
# 0.50  26.73   25.43   23.65
# 0.50  26.87   26.01   23.70
# 0.20  28.17   27.37   25.11
# 0.20  28.07   26.94   25.12
# 0.20  28.11   27.14   25.11
# 0.10  29.20   28.05   26.17
# 0.10  29.49   28.89   26.15
# 0.10  29.07   28.32   26.15
# 0.05  30.17   29.50   27.12
# 0.05  30.14   29.93   27.14
# 0.05  30.12   29.71   27.16
# 0.02  31.35   30.69   28.52
# 0.02  31.35   30.54   28.57
# 0.02  31.35   30.04   28.53
# 0.01  32.55   31.12   29.49
# 0.01  32.45   31.29   29.48
# 0.01  32.28   31.15   29.26
```

For further details please visit: https://github.com/mirzaghaderi/rtpcr

# Contact

Email: gh.mirzaghaderi at uok.ac.ir

# Citation

``` md
To cite the package ‘rtpcr’ in publications, please use:

  Ghader Mirzaghaderi (2025). rtpcr: a package for statistical analysis and graphical presentation of qPCR data in R. PeerJ 13:e20185. https://doi.org/10.7717/peerj.20185

A BibTeX entry for LaTeX users is

  @Article{,
    title = {rtpcr: A package for statistical analysis and graphical presentation of qPCR data in R},
    author = {Ghader Mirzaghaderi},
    journal = {PeerJ},
    volume = {13},
    pages = {e20185},
    year = {2025},
    doi = {10.7717/peerj.20185},
  }
```



