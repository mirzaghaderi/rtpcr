# Changelog

## rtpcr 2.1.5

### New Features

- [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  now supports complex experimental designs. Users can perform ddCt
  analysis across interactions of multiple factors (e.g.,
  `specs = "FactorA | FactorB"`).

### Deprecated

- `mainFactor.column` (Integer) is now deprecated in favor of `specs`.
- `specs` now accepts one or more column names as a character (e.g. “A”
  or “A \| B”, or “A \| B\* C”).
- `mainFactor.level.order` whic already accepted a vector, is now
  deprecated in favor of `calibratorLevel` which accepts one of the
  levels of the main factor (e.g. `A`), Default (NULL) is the first
  level.
- This allows for ddCt calculations to be performed separately within
  each level of a secondary factor or across factor combinations (e.g.,
  `B`, or `B*C`).

## rtpcr 2.1.4

CRAN release: 2026-02-12

### New features

- In
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  function the `se.type` argument was added to control how standard
  errors (SE) are calculated for relative expression (RE = fold change)
  estimates. If set to `"paired.sample"` the se is computed from paired
  differences between factor levels, matching samples by `id`. This is
  automatically used when an `id` random effect is detected in a
  user-provided mixed model. `"two.sample"` computes SE using an
  unpaired samples against the reference level, and `"single.sample"`
  computes SE within each factor level. When a random `id` effect is
  detected in the model, paired standard errors are used automatically
  (with a warning if another `se.type` is requested).

- The default behavior in
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  function for standard error calculation has been updated. Standard
  errors are now calculated from model-based residuals
  (`modelBased_se = TRUE`) by default. Setting `modelBased_se = FALSE`
  restores the previous behavior i.e. direct computing from observed
  wDCt values. For single factor data, both methods are the same. It is
  recommended let it use `modelBased_se = TRUE` (default).

- A function called
  [`plotSingleGene()`](https://mirzaghaderi.github.io/rtpcr/reference/plotSingleGene.md)
  was added that creates a bar plot of relative gene expression (fold
  change) values from single gene analysis showing all pairwise
  significances.

## rtpcr 2.1.3

CRAN release: 2026-02-03

### New features

- In version 2.1.3, optional custom model formula can be supplied by
  user to
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  and
  [`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md)
  functions via `model` argument. If provided, this overrides the
  default formula (single- or multi-factorial CRD or RCBD design based
  on the availability of block and the number of factors). The formula
  uses `wDCt` as the response variable. For mixed models, random effects
  can be supplied (e.g., `wDCt ~ Treatment * (1| id)`).

- Handling missing Ct values for target genes using the
  `set_missing_target_Ct_to_40` function. If `TRUE`, missing target gene
  Ct values become 40; if `FALSE` (default), they become NA. missing Ct
  values of reference genes are always converted to NA.

- The
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  and
  [`ANOVA_DCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DCt.md)
  functions detect singular fits when a user-defined mixed model (lmer)
  is used, and genes that have singular fits after the analysis will be
  reported.

- Because, currently a model (including repeated measure and ancova
  models) can be supplied to
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  by user, `REPEATED_DDCt()` and `ANCOVA_DDCt()`functions were removed.
  Refer to the vignette for a sample of most common models and
  description

- For CRD, RCBD, and factorial experiments under CRD or RCBD designs no
  model is required as one of these models is appropriately selected
  based on the input arguments.

## rtpcr 2.1.2

CRAN release: 2026-01-23

### New features

- The non-parametric
  [`WILCOX_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/WILCOX_DDCt.md)
  function was added to analyze the gene expression (ddCt method) using
  wilcox.test.

- ANCOVA analysis which already was performed via the
  `analysisType = "ancova"` argument of the
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  function, now is performed using the `ANCOVA_DDCt()` new function. The
  `analysisType` argument was removed.

- A function named
  [`compute_wDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/compute_wDCt.md)
  was added which cleans the data and computes weighted delta Ct (wDCt)
  values. Handling missing data by this function has been described in
  the function details.

- Two arguments of the `REPEATED_DDCt()` function were updated in order
  to be similar to those of the
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md)
  function: The `repeatedFactor` and `calibratorLevel` arguments were
  changed to `mainFactor.column` and `mainFactor.level.order`,
  respectively.

- The font size and legend position controls were added to the
  [`efficiency()`](https://mirzaghaderi.github.io/rtpcr/reference/efficiency.md)
  function.

- It is now possible to return the lm formula that is used for variance
  analysis in
  [`ANOVA_DDCt()`](https://mirzaghaderi.github.io/rtpcr/reference/ANOVA_DDCt.md),
  and `REPEATED_DDCt()` functions using
  `object$perGene$gene_name$lm_formula` command.

## rtpcr 2.1.1

CRAN release: 2026-01-08

### New features

- The rtpcr package now accepts any number of references and any number
  of target genes in a single analysis. The package computes the
  geometric mean of reference genes considering efficiency values.

- Missing data can be denoted by NA in the input data frame. Values such
  as “0” and “undetermined” (for any E or Ct) are automatically
  converted to NA before being passed to downstream analyses. For target
  genes, NA values for E or Ct measurements result in NA for the
  corresponding ΔCt for that replicate, which is then propagated through
  downstream statistical analyses. When more than one reference gene is
  used, NA values in either the E or Ct field for any reference gene
  cause that reference gene to be skipped, and the remaining reference
  genes are geometrically averaged.

- The
  [`long_to_wide()`](https://mirzaghaderi.github.io/rtpcr/reference/long_to_wide.md)
  function was added to convert a 4-column qPCR data with long format to
  wide.
