#' ΔΔCt ANOVA analysis
#'
#' Apply ΔΔCt analysis to each target gene
#' in the input data frame. Target and reference genes must be provided as paired
#' efficiency (E) and Ct columns located after the experimental design columns.
#' columns.
#'
#' @param x A data frame containing experimental design columns, target gene
#'   E/Ct column pairs, and reference gene E/Ct column pairs. Reference gene
#'   columns must be located at the end of the data frame. 
#' @param numOfFactors Integer. Number of experimental factor columns
#'   (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes. Each reference gene
#'   must be represented by two columns (E and Ct).
#' @param block Character or \code{NULL}. Name of the blocking factor column.
#' @param mainFactor.column
#' Column index or name of the factor for which relative expression is calculated.
#' When \code{analysisType = "ancova"}, remaining factors are treated as covariates.
#' @param mainFactor.level.order
#' Optional character vector specifying the order of levels for the main factor.
#' If \code{NULL}, the first observed level is used as the calibrator.
#' If provided, the first element of the vector is used as the calibrator level.
#' @param analyseAllTarget Logical or character.
#' If \code{TRUE} (default), all target genes are analysed.
#' Alternatively, a character vector specifying the names (names of their Efficiency columns) of target genes
#' to be analysed.
#' @param analysisType
#' Character string specifying the analysis type; one of \code{"anova"} (default)
#' or \code{"ancova"}.
#' @param p.adj
#' Method for p-value adjustment. See \code{\link[stats]{p.adjust}}.
#' @param plot
#' Logical; if \code{FALSE}, per gene-plots are not generated.
#' @param plotType
#' Plot scale to use: \code{"RE"} for relative expression or
#' \code{"log2FC"} for log2 fold change.
#' 
#' @importFrom stats setNames
#'
#' @details
#' ΔΔCt analysis is performed for 
#' the \code{mainFactor.column} based on a full model factorial 
#' experiment by default. However, if \code{ancova}, the \code{analysisType} argument,
#' analysis of covariance is performed for the levels of the \code{mainFactor.column} and the other factors are 
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
#' @import emmeans
#' @return
#' An object containing expression table, lm models, residuals, raw data and ANOVA table for each gene.
#' \describe{  
#' \item{ΔΔCt combined expression table}{\code{object$combinedFoldChange}}
#' \item{ANOVA table}{\code{object$perGene$gene_name$ANOVA_table}}
#' \item{lm ANOVA}{\code{object$perGene$gene_name$lm_ANOVA}}
#' \item{lm ANCOVA}{\code{object$perGene$gene_name$lm_ANCOVA}}
#' \item{Residuals}{\code{resid(object$perGene$gene_name$lm_ANOVA)}}
#' \item{log2FC_Plot}{\code{object$perGene$gene_name$log2FC_Plot}}
#' \item{RE_Plot}{\code{object$perGene$gene_name$RE_Plot}}
#' }
#' @export
#' 
#' @examples
#' data1 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
#' ANOVA_DDCt(x = data1,
#'            numOfFactors = 2,
#'            numberOfrefGenes = 1,
#'            block = "block",
#'            mainFactor.column = 2,
#'            plot = FALSE,
#'            p.adj = "none")
#'            
#' data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
#' ANOVA_DDCt(
#'            x = data2,
#'            numOfFactors = 1,
#'            numberOfrefGenes = 1,
#'            block = NULL,
#'            mainFactor.column = 1,
#'            plot = FALSE,
#'            p.adj = "none"
#'            )



ANOVA_DDCt <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    mainFactor.column,
    analysisType = "anova",
    mainFactor.level.order = NULL,
    block = NULL,
    p.adj = "none",
    plot = FALSE,
    plotType = "RE",
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
    
    res <- .ANOVA_DDCt_uniTarget(
      x = gene_df,
      numberOfrefGenes = numberOfrefGenes,
      mainFactor.column = mainFactor.column,
      analysisType = analysisType,
      mainFactor.level.order = mainFactor.level.order,
      block = block,
      p.adj = p.adj,
      plot = plot,
      plotType = plotType
    )
    
    res$Fold_Change$gene <- gene_name
    res
  })
  

  # Combine fold-change tables
  combinedFoldChange <- do.call(rbind, lapply(perGene, function(g) g$Fold_Change))
  rownames(combinedFoldChange) <- NULL
  
  combinedFoldChange <- combinedFoldChange[, c(ncol(combinedFoldChange), 1:(ncol(combinedFoldChange) - 1))]
  
  # Print automatically
  cat("\nRelative Expression\n")
  print(combinedFoldChange)
  
  # Return list with all results
  res_list <- list(
    perGene = setNames(perGene, targetNames),
    combinedFoldChange = combinedFoldChange
  )
  
  invisible(res_list)
}


