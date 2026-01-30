# rtpcr 2.1.3

## New features

- In current version, optional custom model formula can be supplied to `ANOVA_DDCt()` and `ANOVA_DCt()` functions. If provided, this overrides the automatic formula (factorial CRD or RCBD design based on the availability of block and the number of factors). The formula uses `wDCt` as the response variable.  For mixed models, random effects can be supplied (e.g., `wDCt ~ Treatment * (1|Block)`). 

- The functions detect singular fits when a user-defined mixed model (lmer) is used, and store a per-gene flag (singular = TRUE/FALSE). Genes that have singular fits after the analysis will be reported.

- Because, currently a model can be supplied `ANOVA_DDCt()` by user, the `REPEATED_DDCt()` function was removed. 

- Refer to the vignette for a sample of most common models and description

- For CRD, RCBD, and factorial experiments under CRD or RCBD designs no model is required as as one of these models is appropriately selected based on the arguments.



# rtpcr 2.1.2

## New features

- The non-parametric `WILCOX_DDCt()` function was added to analyze the gene expression (ddCt method) using wilcox.test. 

- ANCOVA analysis which already was performed via the `analysisType = "ancova"` argument of the `ANOVA_DDCt()` function, now is performed using the `ANCOVA_DDCt()` new function. The `analysisType` argument was removed.

- A function named `compute_wDCt()` was added which cleans the data and computes weighted delta Ct (wDCt) values. Handling missing data by this function has been described in the function details.

- Two arguments of the `REPEATED_DDCt()` function were updated in order to be similar to those of the `ANOVA_DDCt()` function: The `repeatedFactor` and `calibratorLevel` arguments were changed to `mainFactor.column` and `mainFactor.level.order`, respectively.   

- The font size and legend position controls were added to the `efficiency()` function.

- It is now possible to return the lm formula that is used for variance analysis in `ANOVA_DDCt()`, and `REPEATED_DDCt()` functions using `object$perGene$gene_name$lm_formula` command. 



# rtpcr 2.1.1

## New features

- The rtpcr package now accepts any number of references and any number of target genes in a single analysis. The package computes the geometric mean of reference genes considering efficiency values. 

- Missing data can be denoted by NA in the input data frame. Values such as "0" and "undetermined" (for any E or Ct) are automatically converted to NA before being passed to downstream analyses. For target genes, NA values for E or Ct measurements result in NA for the corresponding ΔCt for that replicate, which is then propagated through downstream statistical analyses. When more than one reference gene is used, NA values in either the E or Ct field for any reference gene cause that reference gene to be skipped, and the remaining reference genes are geometrically averaged.

- The `long_to_wide()` function was added to convert a 4-column qPCR data with long format to wide. 




