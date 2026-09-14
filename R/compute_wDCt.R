#' Cleaning data and weighted delta Ct (wDCt) calculation
#'
#' The \code{Compute_wDCt} function cleans the data and computes wDCt for
#' every target gene present in the data. This function is
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
#' @param set_missing_target_Ct_to_40 If \code{TRUE}, missing target gene Ct
#'   values become 40; if \code{FALSE} (default), they become NA.
#'
#' @importFrom stats setNames
#'
#' @details
#' The \code{Compute_wDCt} function computes weighted delta Ct (wDCt) for every
#' target gene in the input data. Missing data can be denoted by NA in the
#' input data frame. Values such as '0' and 'undetermined' (for any E and Ct) are
#' automatically converted to NA. For target genes, NA for E or Ct measurements cause returning NA for
#' the corresponding delta Ct for that replicate (row).
#' If there are more than one reference gene, NA in the place of the E or the Ct value cause
#' skipping that gene and remaining references are geometrically averaged.
#' The \code{Compute_wDCt} function is automatically applied to the expression analysis
#' functions.
#' @return
#' The original data frame along with one weighted delta Ct column per target
#' gene (named \code{wDCt_<gene>}, derived from the target gene's E/Ct column
#' names).
#' @export
#'
#' @examples
#'
#' data <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
#' data
#' Compute_wDCt(x = data,
#'              numOfFactors = 2,
#'              numberOfrefGenes = 3,
#'              block = "block")
#'



Compute_wDCt <- function(x,
                         numOfFactors,
                         numberOfrefGenes,
                         block,
                         set_missing_target_Ct_to_40 = FALSE) {
  
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
  
  x <- .cleanup(x = x, numOfFactors = numOfFactors, numberOfrefGenes = numberOfrefGenes, block = block,
                set_missing_target_Ct_to_40 = set_missing_target_Ct_to_40)
  
  stopifnot(numberOfrefGenes >= 1)
  nRef <- numberOfrefGenes
  nc   <- ncol(x)
  
  # Identify reference gene columns (last 2 * nRef columns)
  ref_E_cols  <- seq(nc - 2 * nRef + 1, nc, by = 2)
  ref_Ct_cols <- seq(nc - 2 * nRef + 2, nc, by = 2)
  
  # Identify design columns (factors + optional block + rep); everything
  # before that and before the reference columns is target gene E/Ct pairs
  nDesign         <- numOfFactors + (if (!is.null(block)) 1 else 0)
  target_start_E  <- nDesign + 2         # skip design cols + rep col
  target_end_Ct   <- ref_E_cols[1] - 1   # last column before references
  
  stopifnot(target_start_E <= target_end_Ct)
  
  target_E_cols  <- seq(target_start_E, target_end_Ct, by = 2)
  target_Ct_cols <- target_E_cols + 1
  nTargets       <- length(target_E_cols)
  
  # Reference matrices
  E_mat  <- as.matrix(x[, ref_E_cols])
  Ct_mat <- as.matrix(x[, ref_Ct_cols])
  
  # Row-wise geometric mean of reference efficiencies (NA-aware)
  geoMeanE <- apply(E_mat, 1, function(z) {
    k <- sum(!is.na(z))
    if (k > 0) {
      # Check for zeros or negative values
      if (any(z <= 0, na.rm = TRUE)) {
        return(NA_real_)
      }
      prod(z, na.rm = TRUE)^(1 / k)
    } else {
      NA_real_
    }
  })
  
  # Reference term (Excel-consistent) - shared across all target genes
  ref_term <- numeric(nrow(x))
  for (r in seq_len(nrow(x))) {
    
    # Check if geoMeanE is valid (positive and finite)
    if (is.na(geoMeanE[r]) || geoMeanE[r] <= 0 || !is.finite(geoMeanE[r])) {
      ref_term[r] <- NA_real_
      next
    }
    
    logE <- log2(geoMeanE[r])
    
    # Check if logE is finite
    if (!is.finite(logE)) {
      ref_term[r] <- NA_real_
      next
    }
    
    tmp <- logE * Ct_mat[r, ]
    
    # Remove NA values
    tmp <- tmp[!is.na(tmp)]
    k   <- length(tmp)
    
    if (k > 0) {
      # Check if any tmp values are infinite
      if (any(!is.finite(tmp))) {
        ref_term[r] <- NA_real_
      } else {
        # Calculate geometric mean
        result <- prod(tmp)^(1 / k)
        # Check if result is finite
        ref_term[r] <- if (is.finite(result)) result else NA_real_
      }
    } else {
      ref_term[r] <- NA_real_
    }
  }
  
  # Derive a base name for each target gene from its E/Ct column names
  target_names <- sub("(_E|_Ct|\\.E|\\.Ct)$", "", colnames(x)[target_Ct_cols])
  # Fall back to generic names if stripping the suffix produced duplicates
  # or empty strings
  if (any(duplicated(target_names)) || any(target_names == "")) {
    target_names <- paste0("target", seq_len(nTargets))
  }
  
  # Compute wDCt for every target gene
  for (i in seq_len(nTargets)) {
    target_term <- log2(x[[target_E_cols[i]]]) * x[[target_Ct_cols[i]]]
    # Convert Inf/-Inf to NA
    target_term[!is.finite(target_term)] <- NA_real_
    
    wDCt_i <- target_term - ref_term
    wDCt_i[!is.finite(wDCt_i)] <- NA_real_
    
    x[[paste0("wDCt_", target_names[i])]] <- wDCt_i
  }
  
  x
}
