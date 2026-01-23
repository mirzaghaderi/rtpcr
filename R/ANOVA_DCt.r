#' Delta Ct ANOVA analysis
#'
#' Performs Delta Ct (dCt) analysis of the data from a 1-, 2-, or 3-factor experiment. 
#' Per-gene statistical grouping is also performed for all treatment (T) combinations.
#' @details
#' The function returns analysis of variance components and the expression table which include these columns: gene (name 
#' of target genes), factor columns, dCt (mean weighted delta Ct for each 
#' treatment combination), RE (relative expression = 2^-dCt), log2FC 
#' (log(2) of relative expression), 
#' LCL (95\% lower confidence level), UCL (95\% upper confidence level),
#' se (standard error of the mean calculated from the weighted delta Ct (wDCt) values of each treatment combination),
#' Lower.se.RE (The lower limit error bar for RE which is 2^(log2(RE) - se)), 
#' Upper.se.RE (The upper limit error bar for RE which is 2^(log2(RE) + se)),
#' Lower.se.log2FC (The lower limit error bar for log2 RE), 
#' Upper.se.log2FC (The upper limit error bar for log2 RE), and sig (per-gene significance grouping letters).
#' 
#' 
#' @param x The input data frame containing experimental design columns, target gene
#' E/Ct column pairs, and reference gene E/Ct column pairs. Reference gene
#' columns must be located at the end of the data frame. See "Input data 
#' structure" in vignettes for details about data structure.
#' @param numOfFactors Integer. Number of experimental factor columns 
#' (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes. Each reference gene
#'   must be represented by two columns (E and Ct).
#' @param analyseAllTarget Logical or character.
#' If \code{TRUE} (default), all detected target genes are analysed.
#' Alternatively, a character vector specifying the names (names of their 
#' Efficiency columns) of target genes to be analysed.
#' @param block  Character. Block column name or \code{NULL}.
#' When a qPCR experiment is done in multiple qPCR plates, 
#' variation resulting from the plates may interfere with the actual amount of 
#' gene expression. One solution is to conduct each plate as a randomized block 
#' so that at least one replicate of each treatment and control is present 
#' on a plate. Block effect is usually considered as random and its interaction 
#' with any main effect is not considered.
#' @param p.adj
#' Method for p-value adjustment. See \code{\link[stats]{p.adjust}}.
#' @param alpha
#' statistical level for comparisons
#' 
#' @importFrom stats setNames
#'
#' @return
#' An object containing expression table, lm models, ANOVA table, residuals, and raw data for each gene:
#' \describe{
#' \item{dCt expression table for all treatment combinations along with the per-gene statistical grouping}{\code{object$relativeExpression}}
#' \item{ANOVA table for treatments}{\code{object$perGene$gene_name$ANOVA_T}}
#' \item{ANOVA table factorial}{\code{object$perGene$gene_name$ANOVA_factorial}}
#' \item{lm ANOVA for tratments}{\code{object$perGene$gene_name$lm_T}}
#' \item{lm ANOVA factorial}{\code{object$perGene$gene_name$lm_factorial}}
#' \item{Residuals}{\code{resid(object$perGene$gene_name$lm_T)}}
#' }
#' @export
#' 
#' @examples
#' data <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
#' res <- ANOVA_DCt(
#'   data,
#'   numOfFactors = 3,
#'   numberOfrefGenes = 1,
#'   block = NULL)



ANOVA_DCt <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    block,
    alpha = 0.05,
    p.adj = "none",
    analyseAllTarget = TRUE
) {
  

  # Basic checks
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!is.numeric(numOfFactors) || numOfFactors < 1) stop("numOfFactors must be a positive integer")
  if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1) stop("numberOfrefGenes must be a positive integer")
  if (!(isTRUE(analyseAllTarget) || is.character(analyseAllTarget))) {
    stop("analyseAllTarget must be TRUE or a character vector of target gene names")
  }
  
  n <- ncol(x)
  

  # Design columns
  nDesign <- if (is.null(block)) numOfFactors + 1 else numOfFactors + 2
  if (nDesign >= n) stop("Not enough columns for target and reference genes")
  designCols <- seq_len(nDesign)
  

  # Reference gene columns
  nRefCols <- 2 * numberOfrefGenes
  refCols  <- (n - nRefCols + 1):n
  

  # Target gene columns
  targetCols <- setdiff(seq_len(n), c(designCols, refCols))
  if (length(targetCols) == 0 || length(targetCols) %% 2 != 0) {
    stop("Target genes must be supplied as E/Ct column pairs")
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
    
    
    
    res <- .ANOVA_DCt_uniTarget(
      x = gene_df,
      numOfFactors = numOfFactors,
      numberOfrefGenes = numberOfrefGenes,
      block = block,
      alpha = alpha,
      p.adj = p.adj
    )
    
    # Add gene name to results
    res$Results$gene <- gene_name
    res
  })
  

  # Combine results for all genes
  relativeExpression <- do.call(rbind, lapply(perGene, function(g) g$Results))
  rownames(relativeExpression) <- NULL
  
  re_col <- which(names(relativeExpression) == "RE")
  
  relativeExpression <- relativeExpression[, c(ncol(relativeExpression), 1:(ncol(relativeExpression) - 1))]
  

  # Print combined results automatically
  cat("\nRelative Expression\n")
  print(relativeExpression)
  

  # Return full structured object invisibly
  invisible(list(
    perGene = setNames(perGene, targetNames),  # All individual gene outputs with models
    relativeExpression = relativeExpression          # Combined table
  ))
}
