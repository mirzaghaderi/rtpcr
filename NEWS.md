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

- Missing data can be denoted by NA in the input data frame. Values such as '0' and 'undetermined' (for any E and Ct) are automatically converted to NA before passing to the downstream analysis. For target genes, NA for E or Ct measurements cause returning NA for the corresponding ΔCt for that replicate which is passed along to downstream statistical analyses. If there are more than one reference genes, NA in the place of the E or the Ct value of any reference gene cause skipping that reference gene and remaining genes are averaged in that replicate.

- The `long_to_wide()` function was added to convert a 4-column qPCR data with long format to wide. 

- The last gene name in the plot output of the `TTEST_DDCt()` function was fixed.




