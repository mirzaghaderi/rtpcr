#' @title Relative expression analysis using the \eqn{\Delta\Delta C_T} method with ANOVA and ANCOVA
#'
#' @description
#' The \code{ANOVA_DDCt} function performs relative expression (RE) analysis based on
#' the \eqn{\Delta\Delta C_T} method using analysis of variance (ANOVA) or
#' analysis of covariance (ANCOVA). It supports uni- and multi-factorial qPCR
#' experimental designs.
#'
#' Bar plots of relative expression (RE) or log2 fold change (log2FC), together
#' with standard errors and confidence intervals, are optionally returned.
#'
#' @details
#' The function calculates weighted Delta Ct (wDCt) values using as many specified
#' reference genes and then performs statistical analysis on the resulting
#' relative expression values.
#'
#' For multi-factorial experiments, relative expression is calculated for the
#' levels of the factor specified by \code{mainFactor.column}.
#'
#' If \code{analysisType = "anova"}, a full factorial ANOVA model is fitted by
#' default.
#'
#' If \code{analysisType = "ancova"}, relative expression is calculated for the
#' levels of \code{mainFactor.column}, while the remaining factor(s), if any, are
#' treated as covariates. In such cases, the ANCOVA table should be examined
#' carefully: a significant interaction between the main factor and a covariate
#' indicates that ANCOVA assumptions are violated and the model may be inappropriate.
#'
#' ANCOVA is typically used when gene expression is influenced by one or more
#' uncontrolled quantitative variables. For example, gene expression may depend
#' on temperature while the primary interest is in treatment or stress effects.
#'
#' The function also supports single-factor experiments, in which case ANOVA
#' reduces to a one-way analysis.
#'
#' @author Ghader Mirzaghaderi
#'
#' @export
#'
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lmerTest
#' @import emmeans
#'
#' @param x
#' A data frame containing experimental conditions, biological replicates,
#' amplification efficiency (E), and Ct values for target and reference genes.
#' Each Ct value should represent the mean of technical replicates.
#'
#' \strong{NOTE:} Each row corresponds to a different biological individual,
#' reflecting a non-repeated-measures experimental design.
#' See the package vignette for details on data structure and column arrangement.
#'
#' @param numberOfrefGenes
#' Integer specifying the number of reference genes used for normalization
#' (must be \eqn{\ge 1}).
#'
#' @param analysisType
#' Character string specifying the analysis type; one of \code{"anova"} (default)
#' or \code{"ancova"}.
#'
#' @param mainFactor.column
#' Column index or name of the factor for which relative expression is calculated.
#' When \code{analysisType = "ancova"}, remaining factors are treated as covariates.
#'
#' @param mainFactor.level.order
#' Optional character vector specifying the order of levels for the main factor.
#' If \code{NULL}, the first observed level is used as the calibrator.
#' If provided, the first element of the vector is used as the calibrator level.
#'
#' @param x.axis.labels.rename
#' Optional character vector used to relabel the x-axis in bar plots.
#'
#' @param block
#' Optional column name specifying a blocking factor.
#' Blocking is commonly used to account for variation between qPCR plates.
#' Block effects are treated as random, and interactions with main effects
#' are not considered.
#'
#' @param p.adj
#' Method for p-value adjustment.
#'
#' @param plot
#' Logical; if \code{FALSE}, plots are not generated.
#'
#' @param plotType
#' Plot scale to use: \code{"RE"} for relative expression or
#' \code{"log2FC"} for log2 fold change.
#'
#' @return
#' A list containing the following components:
#' \describe{
#'   \item{Final_data}{Input data frame augmented with weighted Delta Ct (wDCt) values.}
#'   \item{lm_ANOVA}{Linear model object for ANOVA analysis (if applicable).}
#'   \item{lm_ANCOVA}{Linear model object for ANCOVA analysis (if applicable).}
#'   \item{ANOVA_table}{ANOVA table.}
#'   \item{ANCOVA_table}{ANCOVA table.}
#'   \item{Expression_Table}{Table of RE values, log2FC, p-values, significance codes,
#'   confidence intervals, standard errors, and lower/upper SE limits.}
#'   \item{RE_Plot}{Bar plot of relative expression values for main factor levels.}
#'   \item{log2FC_Plot}{Bar plot of log2 fold change values for main factor levels.}
#' }
#'
#' @references
#' Livak, K. J. and Schmittgen, T. D. (2001).
#' Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR
#' and the Double Delta CT Method.
#' \emph{Methods}, 25(4), 402–408.
#' doi:10.1006/meth.2001.1262
#'
#' Ganger, M. T., Dietz, G. D., and Ewing, S. J. (2017).
#' A common base method for analysis of qPCR data and the application of simple
#' blocking in qPCR experiments.
#' \emph{BMC Bioinformatics}, 18, 1–11.
#'
#' Yuan, J. S., Reed, A., Chen, F., and Stewart, N. (2006).
#' Statistical Analysis of Real-Time PCR Data.
#' \emph{BMC Bioinformatics}, 7, 85.
#'
#' @examples
#' ANOVA_DDCt(data_1factor,
#'   numberOfrefGenes = 1,
#'   mainFactor.column = 1,
#'   block = NULL
#' )
#'
#' ANOVA_DDCt(data_2factor,
#'   numberOfrefGenes = 1,
#'   mainFactor.column = 2,
#'   analysisType = "ancova",
#'   block = NULL
#' )
#'
#' df <- meanTech(Lee_etal2020qPCR, groups = 1:3)
#'
#' ANOVA_DDCt(df,
#'   numberOfrefGenes = 1,
#'   analysisType = "ancova",
#'   mainFactor.column = 2,
#'   plotType = "log2FC",
#'   block = NULL
#' )


ANOVA_DDCt <- function(
    x, 
    numberOfrefGenes = 1, 
    mainFactor.column = 1, 
    analysisType = "anova",
    mainFactor.level.order = NULL, 
    block = NULL, 
    x.axis.labels.rename = "none",
    p.adj = "none",  
    plot = TRUE,
    plotType = "RE") {
  
  
  x <- x[, c(mainFactor.column, (1:ncol(x))[-mainFactor.column])]

  
  # basic checks
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame")
  }
  if (!is.numeric(numberOfrefGenes) || length(numberOfrefGenes) != 1) {
    stop("`numberOfrefGenes` must be a single numeric value")
  }
  if (!is.numeric(mainFactor.column) || length(mainFactor.column) != 1) {
    stop("`mainFactor.column` must be a single numeric value")
  }
  
  
  if (missing(numberOfrefGenes)) {
    stop("argument 'numberOfrefGenes' is missing, with no default")
  }
  if (missing(mainFactor.column)) {
    stop("argument 'mainFactor.column' is missing, with no default")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
  
  

  
  
  if (is.null(mainFactor.level.order)) {
    mainFactor.level.order <- unique(x[,1])
    calibrartor <- x[,1][1]
    on.exit(cat(structure(paste("*** The", calibrartor, "level was used as calibrator.\n"))))
  } else if (any(is.na(match(unique(x[,1]), mainFactor.level.order))) == TRUE){
    stop("The `mainFactor.level.order` doesn't match main factor levels.")
  } else {
    x <- x[order(match(x[,1], mainFactor.level.order)), ]
    #x <- x[order(match(as.character(x[,1]), as.character(mainFactor.level.order))), ]
    calibrartor <- x[,1][1]
    on.exit(cat(structure(paste("*** The", calibrartor, "level was used as calibrator.\n"))))
  }
  
  
  x <- .compute_wDCt(x, numberOfrefGenes, block)
  x[,1] <- factor(
    x[,1],
    levels = mainFactor.level.order
  )

  # get names of factor columns
  factors <- names(x)[vapply(x, is.factor, logical(1))]
  
  
  # Check if there is on target gene
  # if (is.null(block)) {
  #   if(ncol(x) - (2*numberOfrefGenes) - length(factors) -1 > 2){
  #     stop("Only one target gene is allowed!")
  #   }
  # } elseif(block = "block") {
  #   if(ncol(x) - (2*numberOfrefGenes) - length(factors) -2 > 2){
  #     stop("Only one target gene is allowed!")
  # }
  
  
  # Check if there is block
  if (is.null(block)) {
    
    # ANOVA based on factorial design
    formula_ANOVA <- paste("wDCt ~", paste(factors, collapse = " * "), "+ (1 | rep)")
    base::suppressMessages(lmf <- lmerTest::lmer(formula_ANOVA, data = x))
    ANOVA <- stats::anova(lmf)
    # ANCOVA 
    formula_ANCOVA <- paste("wDCt ~", paste(rev(factors), collapse = " + "), "+ (1 | rep)")
    base::suppressMessages(lmc <- lmerTest::lmer(formula_ANCOVA, data = x))
    ANCOVA <- stats::anova(lmc)
    
  } else {
    # If ANOVA based on factorial design was desired with blocking factor:
    formula_ANOVA <- paste("wDCt ~ block +", paste(factors, collapse = " * "), "+ (1 | rep)")
    lmfb <- lmerTest::lmer(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmfb)
    # ANCOVA 
    formula_ANCOVA <- paste("wDCt ~ block +", paste(rev(factors), collapse = " + "), "+ (1 | rep)")
    lmcb <- lmerTest::lmer(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmcb)
  }
  
  
  
  
  
  # Type of analysis: ancova or anova
  if (is.null(block)) {
    if(analysisType == "ancova") {
      lm <- lmc
    } 
    else{
      lm <- lmf
    }
  } else {
    if(analysisType == "ancova") {
      lm <- lmcb
    } 
    else{
      lm <- lmfb
    } 
  }
  
  
  
  
  pp1 <- emmeans(lm, colnames(x)[1], data = x, adjust = p.adj, mode = "satterthwaite")
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  pp3 <- pp2[1:length(mainFactor.level.order)-1,]
  ci <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(unique(x[,1]))-1,]
  pp <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
  
  

  
  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(factor = x[,1], bwDCt = bwDCt), x[,1]),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))  
  
  
  sig <- .convert_to_character(pp$p.value)
  contrast <- pp$contrast
  post_hoc_test <- data.frame(contrast, 
                              RE = 1/(2^-(pp$estimate)),
                              log2FC = log2(1/(2^-(pp$estimate))),
                              pvalue = pp$p.value,
                              sig = sig,
                              LCL = 1/(2^-pp$lower.CL),
                              UCL = 1/(2^-pp$upper.CL),
                              se = se$se[-1])
  
  reference <- data.frame(contrast = mainFactor.level.order[1],
                          RE = 1,
                          log2FC = 0,
                          pvalue = 1, 
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])
  
  tableC <- rbind(reference, post_hoc_test)
  
 
  FINALDATA <- x
  
  tableC$contrast <- as.character(tableC$contrast)
  tableC$contrast <- sapply(strsplit(tableC$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))
  
  if(any(x.axis.labels.rename == "none")){
    tableC
  }else{
    tableC$contrast <- x.axis.labels.rename
  }
  
  
  
  
  tableC$contrast <- factor(tableC$contrast, levels = unique(tableC$contrast))
  contrast <- tableC$contrast
  LCL <- tableC$LCL
  UCL <- tableC$UCL
  REp <- as.numeric(tableC$RE)
  FCp <- as.numeric(tableC$log2FC)
  significance <- tableC$sig
  se <- tableC$se
  
  tableC <- data.frame(tableC, 
                       Lower.se.RE = 2^(log2(tableC$RE) - tableC$se), 
                       Upper.se.RE = 2^(log2(tableC$RE) + tableC$se))  
  ##################################################
  a <- data.frame(tableC, d = 0)

  for (i in 1:length(tableC$RE)) {
    if (tableC$RE[i] < 1) {
      a$Lower.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] - 0.2
    } else {
      a$Lower.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] + 0.2
    }
  }
  pfc1 <- ggplot(a, aes(contrast,RE)) + 
    geom_col() +
    geom_errorbar(aes(ymin = tableC$Lower.se.RE, ymax=tableC$Upper.se.RE), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = tableC$Upper.se.RE + 0.2)) +
    ylab("Relative Expression (DDCt)")
  pfc2 <- ggplot(a, aes(contrast,log2FC)) +
    geom_col() +
    geom_errorbar(aes(ymin = Upper.se, ymax=Lower.se), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = d)) +
    ylab("log2FC")

  tableC <- data.frame(tableC, Lower.se.log2FC = a$Lower.se, Upper.se.log2FC = a$Upper.se)
  ##################################################  
  
  if (is.null(block)) {
    lm_ANOVA <- lmf
    lm_ANCOVA <- lmc
  } else {
    lm_ANOVA <- lmfb
    lm_ANCOVA <- lmcb
  }
  
  #round tableC to 4 decimal places
  tableC <- tableC %>%
    mutate_if(is.numeric, ~ round(., 4))
  
  outlist2 <- structure(list(Final_data = x,
                             lm_ANOVA = lm_ANOVA,
                             lm_ANCOVA = lm_ANCOVA,
                             ANOVA_table = ANOVA,
                             ANCOVA_table = ANCOVA,
                             Fold_Change  = tableC,
                             RE_Plot = pfc1,
                             log2FC_Plot = pfc2), class = "XX")
  
  print.XX <- function(outlist2){
    cat("ANOVA table", "\n")
    print(outlist2$ANOVA_table)
    cat("\n", sep = '', "ANCOVA table", "\n")
    print(outlist2$ANCOVA_table)
    cat("\n", sep = '', "Expression table", "\n")
    print(outlist2$Fold_Change)
    
    
    if (plot == TRUE){
      if(plotType == "RE"){
        cat("RE_Plot\n")
        print(outlist2$RE_Plot)
      }else{
        cat("log2FC_Plot\n")
        print(outlist2$log2FC_Plot)
      }
    }
    invisible(outlist2)
  }
  print.XX(outlist2)
}
