#' Delta Delta Ct ANCOVA analysis
#'
#' Apply Delta Delta Ct (ddCt) analysis to each target gene
#' and performs per-gene statistical analysis.
#'
#' @param x The input data frame containing experimental design columns, 
#' replicates (integer), target gene E/Ct column pairs, and reference gene E/Ct 
#' column pairs. Reference gene columns must be located at the end of the data frame.  
#' See "Input data structure" in vignettes for details about data structure.
#' @param numOfFactors Integer. Number of experimental factor columns
#'   (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes.
#' @param block Character or \code{NULL}. Name of the blocking factor column.
#' When a qPCR experiment is done in multiple qPCR plates, 
#' variation resulting from the plates may interfere with the actual amount of 
#' gene expression. One solution is to conduct each plate as a randomized block 
#' so that at least one replicate of each treatment and control is present 
#' on a plate. Block effect is usually considered as random and its interaction 
#' with any main effect is not considered.
#' @param mainFactor.column
#' Integer. Column index of the factor for which the relative expression analysis is applied.
#' The remaining factors are treated as covariate(s).
#' @param mainFactor.level.order
#' Optional character vector specifying the order of levels for the main factor.
#' If \code{NULL}, the first observed level is used as the calibrator.
#' If provided, the first element of the vector is used as the calibrator level.
#' @param analyseAllTarget Logical or character.
#' If \code{TRUE} (default), all target genes are analysed.
#' Alternatively, a character vector specifying the names (names of their Efficiency columns) of target genes
#' to be analysed.
#' @param p.adj
#' Method for p-value adjustment. See \code{\link[stats]{p.adjust}}.
#' 
#' @importFrom stats setNames
#'
#' @details
#' ddCt analysis of covariance (ANCOVA) is performed for the levels of the \code{mainFactor.column} and the other factors are 
#' treated as covariates. if the interaction between the main factor and the covariate is significant, ANCOVA is not appropriate. 
#' ANCOVA is basically used when a factor is affected by uncontrolled quantitative covariate(s). 
#' For example, suppose that wDCt of a target gene in a plant is affected by temperature. The gene may 
#' also be affected by drought. Since we already know that temperature affects the target gene, we are 
#' interested to know if the gene expression is also altered by the drought levels. We can design an 
#' experiment to understand the gene behavior at both temperature and drought levels at the same time. 
#' The drought is another factor (the covariate) that may affect the expression of our gene under the 
#' levels of the first factor i.e. temperature. The data of such an experiment can be analyzed by ANCOVA 
#' or using ANOVA based on a factorial experiment. ANCOVA is done  
#' even there is only one factor (without covariate or factor  variable).
#' 
#' All the functions for relative expression analysis (including `TTEST_DDCt()`, 
#' `WILCOX_DDCt()`, `ANOVA_DDCt()`, `ANCOVA_DDCt()`, `REPEATED_DDCt()`, and `ANOVA_DCt()`) return the 
#' relative expression table which include fold change and corresponding 
#' statistics. The output of `ANOVA_DDCt()`, `ANCOVA_DDCt()`, `ANCOVA_DDCt()`, `REPEATED_DDCt()`, 
#' and `ANOVA_DCt()` also include lm models, residuals, raw data and ANOVA table 
#' for each gene. 
#' 
#' The expression table returned by `TTEST_DDCt()`, 
#' `WILCOX_DDCt()`, `ANOVA_DDCt()`, `ANCOVA_DDCt()`, and `REPEATED_DDCt()` functions 
#' include these columns: gene (name of target genes), 
#' contrast (calibrator level and contrasts for which the relative expression is computed), 
#' ddCt (mean of weighted delta delta Ct values), RE (relative expression or 
#' fold change = 2^-ddCt), log2FC (log(2) of relative expression or fold change), 
#' pvalue, sig (per-gene significance), LCL (95\% lower confidence level), UCL (95\% upper confidence level),
#' se (standard error of mean calculated from the weighted delta Ct values of each of the main factor levels),
#' Lower.se.RE (The lower limit error bar for RE which is 2^(log2(RE) - se)), 
#' Upper.se.RE (The upper limit error bar for RE which is 2^(log2(RE) + se)),
#' Lower.se.log2FC (The lower limit error bar for log2 RE), and 
#' Upper.se.log2FC (The upper limit error bar for log2 RE)
#'
#' @import emmeans
#' @return
#' An object containing expression table, lm model, residuals, raw data and ANOVA table for each gene:
#' \describe{  
#' \item{ddCt expression table along with per-gene statistical comparison outputs}{\code{object$relativeExpression}}
#' \item{ANOVA table}{\code{object$perGene$gene_name$ANOVA_table}}
#' \item{lm ANOVA}{\code{object$perGene$gene_name$lm}}
#' \item{lm_formula}{\code{object$perGene$gene_name$lm_formula}}
#' \item{Residuals}{\code{resid(object$perGene$gene_name$lm)}}
#' }
#' @export
#' 
#' @references
#' LivakKJ, Schmittgen TD (2001).
#' Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR
#' and the Double Delta CT Method.
#' \emph{Methods}, 25(4), 402–408.
#' doi:10.1006/meth.2001.1262
#'
#' Ganger MT, Dietz GD, and Ewing SJ (2017).
#' A common base method for analysis of qPCR data and the application of
#' simple blocking in qPCR experiments.
#' \emph{BMC Bioinformatics}, 18, 1–11.
#' 
#' Taylor SC, Nadeau K, Abbasi M, Lachance C, Nguyen M, Fenrich, J. (2019). 
#' The ultimate qPCR experiment: producing publication quality, reproducible 
#' data the first time. \emph{Trends in Biotechnology}, 37, 761-774. 
#' 
#' Yuan JS, Reed A, Chen F, Stewart N (2006).
#' Statistical Analysis of Real-Time PCR Data.
#' \emph{BMC Bioinformatics}, 7, 85.
#' 
#' @examples
#' data1 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
#' ANCOVA_DDCt(x = data1,
#'            numOfFactors = 2,
#'            numberOfrefGenes = 2,
#'            block = "block",
#'            mainFactor.column = 2,
#'            p.adj = "none")
#'            
#' data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
#' ANCOVA_DDCt(x = data2,
#'            numOfFactors = 1,
#'            numberOfrefGenes = 1,
#'            block = NULL,
#'            mainFactor.column = 1,
#'            p.adj = "none")
#'            



ANCOVA_DDCt <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    mainFactor.column,
    block,
    mainFactor.level.order = NULL,
    p.adj = "none",
    analyseAllTarget = TRUE
) {
  

  # Basic checks
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!is.numeric(numOfFactors) || numOfFactors < 1) stop("numOfFactors must be a positive integer")
  if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1) stop("numberOfrefGenes must be a positive integer")
  if (!(isTRUE(analyseAllTarget) || is.character(analyseAllTarget))) {
    stop("analyseAllTarget must be TRUE or a character vector of target names")
  }
  
  n <- ncol(x)
  nDesign <- if (is.null(block)) numOfFactors + 1 else numOfFactors + 2
  designCols <- seq_len(nDesign)
  nRefCols <- 2 * numberOfrefGenes
  refCols <- (n - nRefCols + 1):n
  targetCols <- setdiff(seq_len(n), c(designCols, refCols))
  
  if (length(targetCols) == 0 || length(targetCols) %% 2 != 0) {
    stop("Target genes must be supplied as E/Ct column pairs")
  }
  
  targetPairs <- split(targetCols, ceiling(seq_along(targetCols)/2))
  targetNames <- vapply(targetPairs, function(tc) colnames(x)[tc[1]], character(1))
  
  # Subset target genes if requested
  if (!isTRUE(analyseAllTarget)) {
    keep <- targetNames %in% analyseAllTarget
    if (!any(keep)) stop("None of the specified target genes were found in the data.")
    targetPairs <- targetPairs[keep]
    targetNames <- targetNames[keep]
  }
  

  # Analyse each target gene
  perGene <- lapply(seq_along(targetPairs), function(i) {
    
    tc <- targetPairs[[i]]
    gene_name <- targetNames[i]
    
    gene_df <- x[, c(designCols, tc, refCols), drop = FALSE]
    
    res <- .ANOVA(
      x = gene_df,
      numOfFactors = numOfFactors,
      numberOfrefGenes = numberOfrefGenes,
      mainFactor.column = mainFactor.column,
      analysisType = "ancova",
      mainFactor.level.order = mainFactor.level.order,
      block = block,
      p.adj = p.adj
    )
    
    res$Fold_Change$gene <- gene_name
    res
  })
  

  # Combine fold-change tables
  relativeExpression <- do.call(rbind, lapply(perGene, function(g) g$Fold_Change))
  rownames(relativeExpression) <- NULL
  
  relativeExpression <- relativeExpression[, c(ncol(relativeExpression), 1:(ncol(relativeExpression) - 1))]
  
  # Print automatically
  cat("\nRelative Expression\n")
  print(relativeExpression)
  
  # Return list with all results
  res_list <- list(
    perGene = setNames(perGene, targetNames),
    relativeExpression = relativeExpression
  )
  on.exit(cat(paste("The", res_list$relativeExpression$contrast[1], "level was used as calibrator.\n")))
  invisible(res_list)
}


