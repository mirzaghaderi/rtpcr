#' Cleaning data and weighted delta Ct (wDCt) calculation
#'
#' The \code{compute_wDCt} function cleans the data and computes wDCt. This function is 
#' automatically applied to the expression analysis functions like \code{ANOVA_DDCt}, 
#' \code{TTEST_DDCt}, etc. So it should not be applied in advance of expression analysis functions.
#'
#' @param x A data frame containing experimental design columns, replicates (integer), target gene
#'   E/Ct column pairs, and reference gene E/Ct column pairs. Reference gene
#'   columns must be located at the end of the data frame. 
#' @param numOfFactors Integer. Number of experimental factor columns
#'   (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes.
#' @param block Character or \code{NULL}. Name of the blocking factor column.
#' When a qPCR experiment is done in multiple qPCR plates, each plate is considered as a random block 
#' so that at least one replicate of each treatment and control is present 
#' on a plate.
#' 
#' @importFrom stats setNames
#'
#' @details
#' The \code{compute_wDCt} function computes weighted delta Ct (wDCt) for the input data. 
#' Missing data can be denoted by NA in the input data frame. 
#' Values such as '0' and 'undetermined' (for any E and Ct) are
#' automatically converted to NA. For target genes, NA for E or Ct measurements cause returning NA for 
#' the corresponding delta Ct for that replicate (row). 
#' If there are more than one reference gene, NA in the place of the E or the Ct value cause
#' skipping that gene and remaining references are geometrically averaged.
#' The \code{compute_wDCt} function is automatically applied to the expression analysis
#' functions.
#' @return
#' The original data frame along with the weighted delta Ct column.
#' @export
#' 
#' @examples
#'            
#' data <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
#' data
#' compute_wDCt(x = data,
#'              numOfFactors = 2,
#'              numberOfrefGenes = 3,
#'              block = "block")
#' 

compute_wDCt <- function(x, 
                          numOfFactors,
                          numberOfrefGenes, 
                          block) {
  
  if (is.null(block)) {
    x[seq_len(numOfFactors)] <- lapply(
      x[seq_len(numOfFactors)],
      factor
    )
  } else {
    x[seq_len(numOfFactors + 1)] <- lapply(
      x[seq_len(numOfFactors + 1)],
      factor
    )
  }
  
  x <- .cleanup(x, numOfFactors, block)
  
  stopifnot(numberOfrefGenes >= 1)
  nRef <- numberOfrefGenes
  nc   <- ncol(x)
  
  # Identify columns
  ref_E_cols  <- seq(nc - 2 * nRef + 1, nc, by = 2)
  ref_Ct_cols <- seq(nc - 2 * nRef + 2, nc, by = 2)
  
  target_E_col  <- ref_E_cols[1] - 2
  target_Ct_col <- ref_E_cols[1] - 1
  
  # Target term
  target_term <- log2(x[[target_E_col]]) * x[[target_Ct_col]]
  
  # Reference matrices
  E_mat  <- as.matrix(x[, ref_E_cols])
  Ct_mat <- as.matrix(x[, ref_Ct_cols])
  
  # Row-wise geometric mean of reference efficiencies (NA-aware)
  geoMeanE <- apply(E_mat, 1, function(z) {
    k <- sum(!is.na(z))
    if (k > 0) prod(z, na.rm = TRUE)^(1 / k) else NA_real_
  })
  
  # Reference term (Excel-consistent)
  ref_term <- numeric(nrow(x))
  for (r in seq_len(nrow(x))) {
    
    tmp <- log2(geoMeanE[r]) * Ct_mat[r, ]
    k   <- sum(!is.na(tmp))
    
    if (k > 0) {
      ref_term[r] <- prod(tmp, na.rm = TRUE)^(1 / k)
    } else {
      ref_term[r] <- NA_real_
    }
  }
  
  # Final wDCt
  x$wDCt <- target_term - ref_term
  x
}
