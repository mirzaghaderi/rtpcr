#' @title Fold change (\eqn{\Delta \Delta C_T} method) analysis using a model
#' 
#' @description Fold change (\eqn{\Delta \Delta C_T} method) analysis using a model produced by the
#' \code{qpcrANCOVA} or \code{qpcrREPEATED}.
#' 
#' @details The \code{qpcrMeans} performs fold change (\eqn{\Delta \Delta C_T} method) analysis using a model produced by the
#' \code{qpcrANCOVA} or \code{qpcrREPEATED}. The values can be returned for any effects in the model including simple effects,
#' interactions and slicing if an ANOVA model is used, but ANCOVA models returned by rtpcr package only include simple effects.
#' 
#' @author Ghader Mirzaghaderi
#' @export qpcrMeans
#' @import emmeans
#' @param model an `lmer` fitted model object
#' @param specs A character vector specifying the names of the predictors over which FC values are desired
#' @param p.adj Method for adjusting p values
#' @return Table of FC values, significance and confidence limits.
#' 
#' 
#' @examples
#' 
#' # Returning fold change values from a fitted model.
#' # Firstly, result of `qpcrANCOVA` or `qpcrREPEATED` is 
#' # acquired which includes a model object:
#' res <- qpcrANCOVA(data_3factor, numberOfrefGenes = 1, mainFactor.column = 1)
#' 
#' # Returning fold change values of Conc levels from a fitted model:
#' qpcrMeans(res$lm_ANOVA, specs = "Conc")
#' 
#' # Returning fold change values of Conc levels sliced by Type:
#' qpcrMeans(res$lm_ANOVA, specs = "Conc | Type")
#' 
#' # Returning fold change values of Conc levels sliced by Type*SA:
#' qpcrMeans(res$lm_ANOVA, specs = "Conc | (Type*SA)")
#' 
#' # Returning fold change values of Conc
#' qpcrMeans(res$lm_ANOVA, specs = "Conc * Type")


qpcrMeans <- function(model, specs, p.adj = "none"){
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
colnames(emm2)[which(names(emm2) == "estimate")] <- "FC"
#emm2$se = (emm2$UCL - emm2$LCL)/(2*stats::qt(0.975, emm2$df))
emm2$sig <- .convert_to_character(emm2$p.value)
return(emm2)
}



