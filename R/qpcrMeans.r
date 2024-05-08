#' @title Fold change (\eqn{\Delta \Delta C_T} method) analysis using model.
#' 
#' @description Fold change (\eqn{\Delta \Delta C_T} method) analysis using a model produced by the
#' \code{qpcrANCOVA} or \code{qpcrREPEATED}.
#' 
#' @details The \code{qpcrMeans} performs fold change (\eqn{\Delta \Delta C_T} method) analysis using a model produced by the
#' \code{qpcrANCOVA} or \code{qpcrREPEATED}.
#' 
#' @author Ghader Mirzaghaderi
#' @export qpcrMeans
#' @import emmeans
#' @param model an `lmer` fitted model object
#' @param specs A character vector specifying the names of the predictors over which FC values are desired
#' @return Table of FC values, significance and confidence limits.
#' 
#' 
#' @examples
#' 
#' res <- qpcrANCOVA(data_3factor, numberOfrefGenes = 1, mainFactor.column = 1)
#' FCslice(res$lm_ANOVA, specs = "Conc | Type")
#' 
#' FCslice(res$lm_ANOVA, specs = "Conc | (Type*SA)")
#' 
#' FCslice(res$lm_ANOVA, specs = "Conc * Type")


qpcrMeans <- function(model, specs){
  f <- as.formula(paste("pairwise ~", specs))
Pv <- emmeans(ml, f)
emm2 <- confint(emmeans(ml, f))
emm2$contrasts$p.value <- as.data.frame(Pv$contrasts)$p.value
emm2 <- as.data.frame(emm2$contrasts)
emm2$estimate <- 2^emm2$estimate
emm2$lower.CL <- 2^emm2$lower.CL
emm2$upper.CL <- 2^emm2$upper.CL
emm2$contrast <- as.character(emm2$contrast)
emm2$contrast <- sapply(strsplit(emm2$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))
colnames(emm2)[which(names(emm2) == "estimate")] <- "FC"
return(emm2)
}



