#' ΔΔCt ANOVA analysis on repeated measure data
#'
#' @description \code{REPEATED_DDCt} function performs ΔΔCt method
#' analysis of observations repeatedly taken over different time courses. 
#' Data may be obtained over time from a uni- or multi-factorial experiment. Target genes must be provided as paired
#' efficiency (E) and Ct columns followed by the E/Ct column pairs of reference genes.
#'
#' @param x input data frame in which the first column is \code{id}, 
#' followed by the factor column(s) which include at least time. 
#' The first level of time in data frame is used as calibrator or reference level.
#' Additional factor(s) may also be present. Other columns are efficiency and Ct values of target and reference genes.
#' In the \code{id} column, a unique number is assigned to each individual from which samples have been taken over time, 
#' for example see \code{data_repeated_measure_1}, 
#' all the three number 1 indicate one individual which has been sampled over three different time courses.
#' See example data sets or refer \href{../doc/vignette.html}{\code{vignette}}, section "data structure and column arrangement" for details.
#' 
#' @param numOfFactors Integer. Number of experimental factor columns (excluding optional \code{block}).
#' @param repeatedFactor
#' Character string specifying the factor for which fold changes are analysed
#' (commonly \code{"time"}).
#' @param numberOfrefGenes Integer. Number of reference genes. Each reference gene
#'   must be represented by two columns (E and Ct).
#' @param calibratorLevel
#' A level of \code{repeatedFactor} to be used as the calibrator (reference level) which is the reference level or sample that all others are compared to. Examples are untreated 
#' or time 0.
#' @param block
#' Optional blocking factor column name. If supplied, block effects are treated
#' as random effects.
#' @param analyseAllTarget Logical or character. If \code{TRUE} (default), all 
#' detected target genes are analysed. Alternatively, a character 
#' vector specifying the names (names of their Efficiency columns) of target genes to be analysed.
#' @param p.adj Method for p-value adjustment. See \code{\link[stats]{p.adjust}}.
#' @param plot Logical; if \code{FALSE}, plots are not produced.
#' 
#' @importFrom stats setNames
#'
#' @details
#' Column layout requirements for \code{x}:
#' \itemize{
#'   \item Target gene columns: E/Ct column pairs located between design and reference columns
#'   \item Reference gene columns: E/Ct column pairs located at the end
#' }
#'
#'
#' @return
#' An object containing expression table, lm models, residuals, raw data and ANOVA table for each gene.
#' \describe{ 
#' \item{ΔΔCt combined expression table}{\code{object$Relative_Expression_table}}
#' \item{ANOVA table}{\code{object$perGene$gene_name$ANOVA_table}}
#' \item{lm ANOVA}{\code{object$perGene$gene_name$lm}}
#' \item{Residuals}{\code{resid(object$perGene$gene_name$lm)}}
#' \item{log2FC_Plot}{\code{object$perGene$gene_name$log2FC_Plot}}
#' \item{RE_Plot}{\code{object$perGene$gene_name$RE_Plot}}
#' }
#' @export
#' 
#'
#' @examples
#' data1 <- read.csv(system.file("extdata", "data_repeated_measure_1.csv", package = "rtpcr"))
#' REPEATED_DDCt(
#'   data1,
#'   numOfFactors = 1,
#'   numberOfrefGenes = 1,
#'   repeatedFactor = "time",
#'   calibratorLevel = "1",
#'   block = NULL)
#'
#'
#'
#' data2 <- read.csv(system.file("extdata", "data_repeated_measure_2.csv", package = "rtpcr"))
#' REPEATED_DDCt(
#'   data2,
#'   numOfFactors = 2,
#'   numberOfrefGenes = 1,
#'   repeatedFactor = "time", 
#'   calibratorLevel = "1",
#'   block = NULL,
#'   p.adj = "none",
#'   plot = FALSE,
#'   analyseAllTarget = TRUE)
#' 


REPEATED_DDCt <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    repeatedFactor,
    calibratorLevel,
    block = NULL,
    p.adj = "none",
    plot = FALSE,
    analyseAllTarget = TRUE
) {
  
  
  col_to_rename <- if (is.null(block)) numOfFactors + 1 else numOfFactors + 2
  colnames(x)[col_to_rename] <- "id"
  
  
  # Basic checks
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!is.numeric(numOfFactors) || numOfFactors < 1)
    stop("numOfFactors must be a positive integer")
  if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1)
    stop("numberOfrefGenes must be a positive integer")
  if (!(isTRUE(analyseAllTarget) || is.character(analyseAllTarget)))
    stop("analyseAllTarget must be TRUE or a character vector")
  if (!repeatedFactor %in% colnames(x))
    stop("The specified repeatedFactor is not in x; Or numOfFactors is not correct.")
  if (!is.null(block) && !block %in% colnames(x))
    stop("The specified block column was not found in x")
  
  colnames(x)[colnames(x) == repeatedFactor] <- "Time_"
  
  n <- ncol(x)
  
  
  # Reference gene columns (ALWAYS last)
  nRefCols <- 2 * numberOfrefGenes
  if (nRefCols >= n)
    stop("Not enough columns for reference genes")
  
  
  targetCols <- (numOfFactors + 2 + !is.null(block)) : (n - nRefCols)
  
  
  #refCols <- (n - nRefCols + 1):n
  refCols <- (numOfFactors + 2 + length(targetCols) + !is.null(block)) : n
  
  # Target gene columns (just before ref genes)
  preRefCols <- if (min(refCols) > 1) {
    seq_len(min(refCols) - 1)
  } else {
    integer(0)
  }
  
  # number of target columns must be even (E/Ct pairs)
  nTargetCols <- length(targetCols)
  if (nTargetCols <= 0 || nTargetCols %% 2 != 0) {
    stop(
      "Gene columns must be provided as E/Ct pairs",
      call. = FALSE
    )
  }
  
  
  
  # Design columns (everything before target pairs)
  designCols <- setdiff(seq_len(n), c(targetCols, refCols))
  
  expectedDesign <- numOfFactors + 1 + !is.null(block)
  if (length(designCols) != expectedDesign) {
    stop(
      "Mismatch between numOfFactors and detected design columns.\n",
      "Expected: ", expectedDesign,
      " | Found: ", length(designCols),
      call. = FALSE
    )
  }
  
  
  move_col <- targetCols[1] - 1
  
  # columns before target genes
  before_target <- seq_len(move_col)
  
  # reorder columns: moved column first, then remaining before-target columns,
  # then target + reference columns unchanged
  new_order <- c(
    move_col,
    setdiff(before_target, move_col),
    (targetCols[1]):ncol(x)
  )
  
  x <- x[, new_order, drop = FALSE]
  
  
  # arrange rows based on columns before target genes
  ord <- do.call(order, x[, seq_len(length(before_target)), drop = FALSE])
  x <- x[ord, , drop = FALSE]
  
  rownames(x) <- NULL
  
  
  # Target gene pairing and names
  
  targetPairs <- split(targetCols, ceiling(seq_along(targetCols) / 2))
  targetNames <- vapply(
    targetPairs,
    function(tc) colnames(x)[tc[1]],
    character(1)
  )
  
  
  # Subset target genes if requested
  if (!isTRUE(analyseAllTarget)) {
    keep <- targetNames %in% analyseAllTarget
    if (!any(keep))
      stop("None of the specified target genes were found.")
    targetPairs <- targetPairs[keep]
    targetNames <- targetNames[keep]
  }
  
  
  ## Analyse each target gene
  perGene <- lapply(seq_along(targetPairs), function(i) {
    
    tc <- targetPairs[[i]]
    gene_name <- targetNames[i]
    
    gene_df <- x[, c(designCols, tc, refCols), drop = FALSE]
    
    target_in_gene_df <- (length(designCols) + 1):(length(designCols) + length(tc))
    
    if (all(is.na(gene_df[, target_in_gene_df]))) {
      warning("Skipping target gene ", gene_name, " (all NA)")
      return(NULL)
    }
    
    res <- .REPEATED_DDCt_uniTarget(
      x = gene_df,
      numberOfrefGenes = numberOfrefGenes,
      repeatedFactor = repeatedFactor,
      block = block,
      calibratorLevel = calibratorLevel,
      p.adj = p.adj,
      plot = plot
    )
    
    res$Relative_Expression_table$gene <- gene_name
    res
  })
  
  perGene <- Filter(Negate(is.null), perGene)
  if (length(perGene) == 0)
    stop("No target genes could be analysed")
  
  
  # Combine Relative Expression tables
  combinedFoldChange <- do.call(
    rbind,
    lapply(perGene, function(g) g$Relative_Expression_table)
  )
  
  rownames(combinedFoldChange) <- NULL
  combinedFoldChange <- combinedFoldChange[
    , c(ncol(combinedFoldChange), 1:(ncol(combinedFoldChange) - 1))
  ]
  
  cat("\nCombined Relative Expression Table (all genes)\n")
  print(combinedFoldChange)
  
  
  # Return object
  invisible(list(
    perGene = setNames(perGene, targetNames),
    combinedFoldChange = combinedFoldChange
  ))
}
