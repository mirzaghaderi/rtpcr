#' Multi-target Delta-Delta Ct ANOVA analysis
#'
#' Performs Delta-Delta Ct (DDCt) analysis for multiple target genes by
#' applying DDCt method to each target gene present
#' in the input data frame. Target genes must be provided as paired
#' efficiency (E) and Ct columns located between the experimental design
#' columns and the reference gene columns.
#'
#' @param x A data frame containing experimental design columns, target gene
#'   E/Ct column pairs, and reference gene E/Ct column pairs. Reference gene
#'   columns must be located at the end of the data frame.
#' @param NumOfFactors Integer. Number of experimental factor columns
#'   (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes. Each reference gene
#'   must be represented by two columns (E and Ct).
#' @param block Character or \code{NULL}. Name of the blocking factor column.
#' @param mainFactor.column
#' Column index or name of the factor for which relative expression is calculated.
#' When \code{analysisType = "ancova"}, remaining factors are treated as covariates.
#'   If \code{NULL}, no blocking factor is used.
#' @param mainFactor.level.order
#' Optional character vector specifying the order of levels for the main factor.
#' If \code{NULL}, the first observed level is used as the calibrator.
#' If provided, the first element of the vector is used as the calibrator level.
#' @param analyseAllTarget Logical or character.
#' If \code{TRUE} (default), all detected target genes are analysed.
#' Alternatively, a character vector specifying the names (names of their Efficiency columns) of target genes
#' to be analysed.
#' @param analysisType
#' Character string specifying the analysis type; one of \code{"anova"} (default)
#' or \code{"ancova"}.
#' @param x.axis.labels.rename
#' Optional character vector used to relabel the x-axis in bar plots.
#' @param p.adj
#' Method for p-value adjustment.
#' @param plot
#' Logical; if \code{FALSE}, plots are not generated.
#' @param plotType
#' Plot scale to use: \code{"RE"} for relative expression or
#' \code{"log2FC"} for log2 fold change.
#' 
#' @importFrom stats setNames
#'
#' @details
#' Column layout requirements for \code{x}:
#' \itemize{
#'   \item Target gene columns: E/Ct pairs located between design and reference columns
#'   \item Reference gene columns: \code{2 * numberOfrefGenes} columns located at the end
#' }
#'
#' @import emmeans
#' @return
#' A data frame combining DDCt method expression tables for all target genes.
#' @export
#' 
#' @examples
#' data1 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
#' ANOVA_DDCt(x = data1,
#'            NumOfFactors = 2,
#'            numberOfrefGenes = 1,
#'            block = "block",
#'            mainFactor.column = 2,
#'            plot = FALSE,
#'            p.adj = "none")
#'            
#' data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
#' ANOVA_DDCt(
#'            x = data2,
#'            NumOfFactors = 1,
#'            numberOfrefGenes = 1,
#'            block = NULL,
#'            mainFactor.column = 1,
#'            plot = FALSE,
#'            p.adj = "none"
#'            )



ANOVA_DDCt <- function(
    x,
    NumOfFactors,
    numberOfrefGenes,
    mainFactor.column,
    analysisType = "anova",
    mainFactor.level.order = NULL,
    block = NULL,
    x.axis.labels.rename = "none",
    p.adj = "none",
    plot = FALSE,
    plotType = "RE",
    analyseAllTarget = TRUE
) {
  

  # Basic checks
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!is.numeric(NumOfFactors) || NumOfFactors < 1) stop("NumOfFactors must be a positive integer")
  if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1) stop("numberOfrefGenes must be a positive integer")
  if (!(isTRUE(analyseAllTarget) || is.character(analyseAllTarget))) {
    stop("analyseAllTarget must be TRUE or a character vector of target names")
  }
  
  n <- ncol(x)
  nDesign <- if (is.null(block)) NumOfFactors + 1 else NumOfFactors + 2
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
      x.axis.labels.rename = x.axis.labels.rename,
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


