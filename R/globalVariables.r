#' Internal global variables
#'
#' These objects are declared to avoid R CMD check notes about
#' non-standard evaluation (e.g., in ggplot2).
#'
#' @keywords internal
#' @name globalVariables
NULL

utils::globalVariables(c(
  "ymin", "ymax", "lower", "upper", "sig", "d",
  "contrast", "RE", "log2FC", "Lower.se", "Upper.se", 
  "Lower.se.RE", "Upper.se.RE", "factors", "as.formula", "tail", "Gene",
  "se", "pvalue", "ave", "reshape", "lm", "formula", "anova", "wDCt", 
  "na.exclude", "residuals", "block", "bwDCt"
))
