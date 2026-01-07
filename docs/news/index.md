# Changelog

## rtpcr 2.1.1

### New features

- Introduced a unified input data structure for all expression analysis
  functions.

- Removed restrictions on the number of target or reference genes.
  Weighted Ct values are now calculated from efficiency and Ct values
  using the geometric mean of the specified reference genes.

- Added support for analyzing all genes or selected subsets. The
  returned object includes the final expression table for all genes as
  well as per-gene results.

- Added publication-ready graphical outputs via the `plotFacto()`
  function for one- to three-factor experimental designs.

- Improved handling of missing values for more robust analyses.
