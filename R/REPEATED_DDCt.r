#' \eqn{\Delta\Delta C_T} ANOVA analysis on repeated measure data
#'
#' @description \code{REPEATED_DDCt} function performs \eqn{\Delta \Delta C_T} method
#' analysis of observations repeatedly taken over different time courses. 
#' Data may be obtained over time from a uni- or multi-factorial experiment. Target genes must be provided as paired
#' efficiency (E) and Ct columns followed by the columns and the reference gene columns.
#'
#' @param x input data frame in which the first column is id, 
#' followed by the factor column(s) which include at least time. 
#' The first level of time in data frame is used as calibrator or reference level.
#' Additional factor(s) may also be present. Other columns are efficiency and Ct values of target and reference genes.
#'  \strong{NOTE:} In the "id" column, a unique number is assigned to each individual from which samples have been taken over time, 
#' for example see \code{data_repeated_measure_1}, 
#' all the three number 1 indicate one individual which has been sampled over three different time courses.
#' See example data sets or refer \href{../doc/vignette.html}{\code{vignette}}, section "data structure and column arrangement" for details.
#' 
#' @param NumOfFactors Integer. Number of experimental factor columns
#'   (excluding optional \code{block}).
#' @param factor
#' Character string specifying the factor for which fold changes are analysed
#' (commonly \code{"time"}).
#' @param numberOfrefGenes Integer. Number of reference genes. Each reference gene
#'   must be represented by two columns (E and Ct).
#' @param calibratorLevel
#' A level of \code{factor} to be used as the calibrator (reference level) which is the reference level or sample that all others are compared to. Examples are untreated 
#' or time 0.
#' @param block
#' Optional blocking factor column name. If supplied, block effects are treated
#' as random effects.
#' @param x.axis.labels.rename
#' Optional character vector used to replace x-axis labels in the bar plot.
#' @param analyseAllTarget Logical or character.
#'   If \code{TRUE} (default), all detected target genes are analysed.
#'   Alternatively, a character vector specifying the names (names of their Efficiency columns) of target genes
#'   to be analysed.
#' @param p.adj
#' Method for p-value adjustment.
#' @param plot
#' Logical; if \code{FALSE}, plots are not produced.
#' 
#' @importFrom stats setNames
#'
#' @details
#' Column layout requirements for \code{x}:
#' \itemize{
#'   \item Target gene columns: E/Ct pairs located between design and reference columns
#'   \item Reference gene columns: columns located at the end
#' }
#'
#'
#' @return
#' An object containing expression table, lm models, residuals, raw data and ANOVA table for each gene.
#' \describe{ 
#' \item{\eqn{\Delta\Delta C_T} combined expression table}{\code{object$Relative_Expression_table}}
#' \item{ANOVA table}{\code{object$perGene[["gene_name"]]$ANOVA_table}}
#' \item{lm ANOVA}{\code{object$perGene[["gene_name"]]$lm}}
#' \item{Residuals}{\code{resid(object$perGene[["gene_name"]]$lm)}}
#' \item{log2FC_Plot}{\code{object$perGene[["gene_name"]]$log2FC_Plot}}
#' \item{RE_Plot}{\code{object$perGene[["gene_name"]]$RE_Plot}}
#' }
#' @export
#' 
#'
#' @examples
#' data <- read.csv(system.file("extdata", "data_repeated_measure_1.csv", package = "rtpcr"))
#' REPEATED_DDCt(
#'   data,
#'   NumOfFactors = 1,
#'   numberOfrefGenes = 1,
#'   factor = "time",
#'   calibratorLevel = "1",
#'   block = NULL)




REPEATED_DDCt <- function(
    x,
    NumOfFactors,
    numberOfrefGenes,
    factor, 
    calibratorLevel,
    block = NULL,
    x.axis.labels.rename = "none",
    p.adj = "none",
    plot = FALSE,
    analyseAllTarget = TRUE
) {
  

  # Basic checks
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!is.numeric(NumOfFactors) || NumOfFactors < 1) stop("NumOfFactors must be a positive integer")
  if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1) stop("numberOfrefGenes must be a positive integer")
  if (!(isTRUE(analyseAllTarget) || is.character(analyseAllTarget))) {
    stop("analyseAllTarget must be TRUE or a character vector of target gene names")
  }
  if (!factor %in% colnames(x)) stop("The specified time factor column was not found in x")
  if (!is.null(block) && !block %in% colnames(x)) stop("The specified block column was not found in x")
  
  n <- ncol(x)
  

  # Design columns
  # id + NumOfFactors + time (+ block if present)
  nDesign <- if (is.null(block)) NumOfFactors + 1 else NumOfFactors + 2
  if (nDesign >= n) stop("Not enough columns for target and reference genes")
  designCols <- seq_len(nDesign)
  

  # Reference gene columns
  nRefCols <- 2 * numberOfrefGenes
  if (nRefCols >= n) stop("Not enough columns for reference genes")
  refCols <- (n - nRefCols + 1):n
  

  # Target gene columns
  targetCols <- setdiff(seq_len(n), c(designCols, refCols))
  if (length(targetCols) == 0 || length(targetCols) %% 2 != 0) {
    stop("Target gene columns must be provided as E/Ct pairs")
  }
  
  targetPairs <- split(targetCols, ceiling(seq_along(targetCols) / 2))
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
    
    # Skip genes with no usable data
    if (all(is.na(gene_df[, tc]))) {
      warning("Skipping target gene ", gene_name, " (all values are NA)")
      return(NULL)
    }
    
    res <- .REPEATED_DDCt_uniTarget(
      x = gene_df,
      numberOfrefGenes = numberOfrefGenes,
      factor = factor,
      block = block,
      calibratorLevel = calibratorLevel,
      x.axis.labels.rename = x.axis.labels.rename,
      p.adj = p.adj,
      plot = plot
    )
    
    # Add gene name to Relative Expression table
    res$Relative_Expression_table$gene <- gene_name
    res
  })
  
  # Remove skipped genes
  perGene <- Filter(Negate(is.null), perGene)
  if (length(perGene) == 0) stop("No target genes could be analysed")
  

  # Combine Relative Expression tables
  combinedFoldChange <- do.call(rbind, lapply(perGene, function(g) g$Relative_Expression_table))
  rownames(combinedFoldChange) <- NULL
  
  combinedFoldChange <- combinedFoldChange[, c(ncol(combinedFoldChange), 1:(ncol(combinedFoldChange) - 1))]

  # Print combined Fold Change / RE table automatically
  cat("\nCombined Relative Expression Table (all genes)\n")
  print(combinedFoldChange)
  

  # Return full structured object invisibly
  invisible(list(
    perGene = setNames(perGene, targetNames),   # All individual gene outputs with models
    combinedFoldChange = combinedFoldChange     # Combined table
  ))
}
