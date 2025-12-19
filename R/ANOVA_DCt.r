#' Multi-target Delta-Delta Ct ANOVA analysis
#'
#' Performs Delta-Delta Ct (DCt) analysis for multiple target genes by
#' applying DCt method to each target gene. Target genes must be provided as paired
#' efficiency (E) and Ct columns followed by the columns and the reference gene columns.
#'
#' @param x A data frame containing experimental design columns, target gene
#'   E/Ct column pairs, and reference gene E/Ct column pairs. Reference gene
#'   columns must be located at the end of the data frame.
#' @param NumOfFactors Integer. Number of experimental factor columns
#'   (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes. Each reference gene
#'   must be represented by two columns (E and Ct).
#' @param analyseAllTarget Logical or character.
#'   If \code{TRUE} (default), all detected target genes are analysed.
#'   Alternatively, a character vector specifying the names (names of their Efficiency columns) of target genes
#'   to be analysed.
#' @param block Character or \code{NULL}. Name of the blocking factor column.
#' @param adjust
#' Method for p-value adjustment.
#' @param alpha
#' statistical level for comparisons
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
#'
#' @return
#' A data frame combining DCt method expression tables for all target genes.
#' The output includes all columns.
#' @export
#' 
#' @examples
#' res <- ANOVA_DCt(
#'   data_3factor,
#'   NumOfFactors = 3,
#'   numberOfrefGenes = 1,
#'   block = NULL
#' )



ANOVA_DCt <- function(
    x,
    NumOfFactors,
    numberOfrefGenes,
    block = NULL,
    alpha = 0.05,
    adjust = "none",
    analyseAllTarget = TRUE
) {
  
  # -----------------------------
  # Basic checks
  # -----------------------------
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!is.numeric(NumOfFactors) || NumOfFactors < 1) stop("NumOfFactors must be a positive integer")
  if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1) stop("numberOfrefGenes must be a positive integer")
  if (!(isTRUE(analyseAllTarget) || is.character(analyseAllTarget))) {
    stop("analyseAllTarget must be TRUE or a character vector of target gene names")
  }
  
  n <- ncol(x)
  
  # -----------------------------
  # Design columns
  # -----------------------------
  nDesign <- if (is.null(block)) NumOfFactors + 1 else NumOfFactors + 2
  if (nDesign >= n) stop("Not enough columns for target and reference genes")
  designCols <- seq_len(nDesign)
  
  # -----------------------------
  # Reference gene columns
  # -----------------------------
  nRefCols <- 2 * numberOfrefGenes
  refCols  <- (n - nRefCols + 1):n
  
  # -----------------------------
  # Target gene columns
  # -----------------------------
  targetCols <- setdiff(seq_len(n), c(designCols, refCols))
  if (length(targetCols) == 0 || length(targetCols) %% 2 != 0) {
    stop("Target genes must be supplied as E/Ct column pairs")
  }
  
  targetPairs <- split(targetCols, ceiling(seq_along(targetCols) / 2))
  targetNames <- vapply(targetPairs, function(tc) colnames(x)[tc[1]], character(1))
  
  # -----------------------------
  # Subset target genes if requested
  # -----------------------------
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
    
    res <- .ANOVA_DCt_uniTarget(
      x = gene_df,
      numberOfrefGenes = numberOfrefGenes,
      block = block,
      alpha = alpha,
      adjust = adjust
    )
    
    # Add gene name to results
    res$Results$gene <- gene_name
    res
  })
  

  # Combine results for all genes
  combinedResults <- do.call(rbind, lapply(perGene, function(g) g$Results))
  rownames(combinedResults) <- NULL
  

  # Print combined results automatically
  cat("\nCombined Expression Table (all genes)\n")
  print(combinedResults)
  

  # Return full structured object invisibly
  invisible(list(
    perGene = setNames(perGene, targetNames),  # All individual gene outputs with models
    combinedResults = combinedResults          # Combined table
  ))
}
