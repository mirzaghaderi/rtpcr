#' @title Pairwise comparisons of relative expression values (\eqn{\Delta \Delta C_T}) using a fitted model
#'
#' @description
#' Performs relative expression (fold change) analysis based on the
#' \eqn{\Delta \Delta C_T} method using a fitted model object produced by
#' \code{ANOVA_DDCt()} or \code{REPEATED_DDCt()}.
#'
#' @details
#' The \code{Means_DDCt} function performs pairwise comparisons of relative expression values fo all combinations using
#' estimated marginal means derived from a fitted model.
#' For ANOVA models, FC values can be obtained for main effects,
#' interactions, and sliced (simple) effects.
#' For ANCOVA models returned by the \pkg{rtpcr} package, only simple
#' effects are supported.
#'
#' Internally, this function relies on the \pkg{emmeans} package to
#' compute marginal means and contrasts, which are then back-transformed
#' to fold change values using the \eqn{\Delta \Delta C_T} framework.
#'
#' @author
#' Ghader Mirzaghaderi
#'
#' @export
#'
#' @import emmeans
#'
#' @param model
#' A fitted model object (typically an \code{lmer} or \code{lm} object)
#' created by \code{ANOVA_DDCt()} or \code{REPEATED_DDCt()}.
#'
#' @param specs
#' A character string or character vector specifying the predictors or
#' combinations of predictors over which relative expression values are desired.
#' This argument follows the specification syntax used by
#' \code{emmeans::emmeans()} (e.g., \code{"Factor"},
#' \code{"Factor1 | Factor2"}).
#'
#' @param p.adj
#' Character string specifying the method for adjusting p-values.
#' See \code{\link[stats]{p.adjust}} for available options.
#'
#' @return
#' A data frame containing estimated relative expression values, confidence
#' intervals, p-values, and significance levels derived from the fitted
#' model.
#'
#' @examples
#'
#' # Obtain a fitted model from ANOVA_DDCt
#' res <- ANOVA_DDCt(
#'   data_3factor,
#'   NumOfFactors = 3,
#'   numberOfrefGenes = 1,
#'   mainFactor.column = 1,
#'   block = NULL)
#'
#' # Relative expression values for Type main effect
#' Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Type")
#'
#' # Relative expression values for Concentration main effect
#' Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Conc")
#'
#' # Relative expression values for Concentration sliced by Type
#' Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Conc | Type")
#'
#' # Relative expression values for Concentration sliced by Type and SA
#' Means_DDCt(res$perGene$E_PO$lm_ANOVA, specs = "Conc | Type * SA")




Means_DDCt <- function(model, specs, p.adj = "none"){
  
  if(any(c("lmerModLmerTest", "lm") %in% class(model)) == FALSE){
    stop(deparse(substitute(model)), " is not an accepted model object created by ANOVA_DDCt, qpcrANOVARE, or REPEATED_DDCt functions!")
  }
  
  f <- stats::as.formula(paste("pairwise ~", specs))
  base::suppressMessages(Pv <- emmeans(model, f, adjust = p.adj))
  base::suppressMessages(emm2 <- stats::confint(emmeans(model, f, adjust = p.adj)))
emm2$contrasts$p.value <- as.data.frame(Pv$contrasts)$p.value
emm2 <- as.data.frame(emm2$contrasts)
emm2$estimate <- 2^emm2$estimate
emm2$lower.CL <- 2^emm2$lower.CL
emm2$upper.CL <- 2^emm2$upper.CL
emm2$contrast <- as.character(emm2$contrast)
emm2$contrast <- sapply(strsplit(emm2$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))
colnames(emm2)[which(names(emm2) == "lower.CL")] <- "LCL"
colnames(emm2)[which(names(emm2) == "upper.CL")] <- "UCL"
#emm2$se = (emm2$UCL - emm2$LCL)/(2*stats::qt(0.975, emm2$df))
emm2$sig <- .convert_to_character(emm2$p.value)
colnames(emm2)[which(names(emm2) == "estimate")] <- "RE"
colnames(emm2)[which(names(emm2) == "FC")] <- "RE"
return(emm2)
}


