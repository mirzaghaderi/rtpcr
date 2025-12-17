#' @title Relative expression analysis using the \eqn{\Delta C_T} method with ANOVA
#'
#' @description
#' The \code{ANOVA_DCt} function performs analysis of variance (ANOVA) on
#' relative expression values calculated using the \eqn{\Delta C_T} method.
#' Expression levels are normalized using one or more reference genes and
#' analyzed across all combinations of experimental factor levels.
#'
#' @details
#' Relative expression (RE) values are calculated using the \eqn{\Delta C_T}
#' method, where Ct values of target genes are normalized to reference gene(s).
#' The resulting weighted Delta Ct (wDCt) values are then analyzed using ANOVA.
#'
#' The function supports uni- and multi-factorial experimental designs.
#' Blocking factors (e.g. qPCR plates) can optionally be included to account
#' for technical variation. Each row of the input data represents an independent
#' biological individual, corresponding to a non-repeated-measures experiment.
#'
#' @author Ghader Mirzaghaderi
#'
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import lmerTest
#' @import multcomp
#' @import emmeans
#' @import multcompView
#'
#' @param x
#' A data frame structured as described in the package vignette, containing
#' experimental condition columns, amplification efficiency (E) and Ct values
#' for target and reference genes. Each Ct value should represent the mean of
#' technical replicates.
#'
#' \strong{NOTE:} Each row corresponds to a separate biological individual,
#' reflecting a non-repeated-measures experimental design.
#'
#' @param numberOfrefGenes
#' Integer specifying the number of reference genes used for normalization
#' (must be \eqn{\ge 1}).
#'
#' @param block
#' Character string or \code{NULL}. If provided, this specifies the column name
#' in \code{x} corresponding to a blocking factor (e.g. qPCR plate).
#' Blocking is used to reduce technical variation arising from experimental
#' conditions such as plate-to-plate differences.
#'
#' @param alpha
#' Significance level used for compact letter display (CLD);
#' default is \code{0.05}.
#'
#' @param adjust
#' P-value adjustment method passed to \code{emmeans} and \code{cld}.
#'
#' @return
#' A list containing the following components:
#' \describe{
#'   \item{Final_data}{Input data frame augmented with weighted Delta Ct (wDCt) values.}
#'   \item{lm}{Fitted linear model object, including ANOVA results.}
#'   \item{ANOVA}{ANOVA table based on a completely randomized design (CRD).}
#'   \item{Result}{Result table containing treatment and factor levels, relative
#'   expression (RE), log2 fold change (log2FC), confidence limits (LCL, UCL),
#'   compact letter display for pairwise comparisons, and standard errors with
#'   corresponding lower and upper limits.}
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
#' # If the data include technical replicates, calculate means first:
#' # df <- meanTech(data_3factor, groups = 1:3)
#'
#' # One-factor or multi-factor ANOVA without blocking
#' ANOVA_DCt(
#'   data_3factor,
#'   numberOfrefGenes = 1,
#'   block = NULL
#' )
#'
#' # ANOVA with blocking factor
#' ANOVA_DCt(
#'   data_2factorBlock,
#'   numberOfrefGenes = 1,
#'   block = "Block"
#' )


ANOVA_DCt <- function(x,
                      numberOfrefGenes,
                      block,
                      alpha = 0.05,
                      adjust = "none") {
  
  ## basic argument checks
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
    if (!is.numeric(numberOfrefGenes) || length(numberOfrefGenes) != 1) {
    stop("`numberOfrefGenes` must be a single numeric value")
  }
  
  
  x <- .compute_wDCt(x, numberOfrefGenes, block)
  # get names of factor columns
  factors <- names(x)[vapply(x, is.factor, logical(1))]
  
  
  ##build treatment factor T and fit lm
    if (is.null(block)) {
      x$T <- do.call(paste, c(x[1:length(factors)], sep = ":"))
      x$T <- as.factor(x$T)
      lm_fit <- stats::lm(wDCt ~ T, data = x)
      anovaCRD <- stats::anova(lm_fit)
    } else {
      x$T <- do.call(paste, c(x[1:length(factors)], sep = ":"))
      x$T <- as.factor(x$T)
      lm_fit <- stats::lm(wDCt ~ block + T, data = x)
      anovaCRD <- stats::anova(lm_fit)
    }

  
  ## emmeans / multiple comparisons
  emg <- emmeans::emmeans(lm_fit, pairwise ~ T, mode = "satterthwaite")
  # use cld() on the emmeans object (the first element)
  meanPairs <- multcomp::cld(emg[[1]], adjust = adjust, alpha = alpha, reversed = FALSE, Letters = letters)
  # meanPairs typically contains columns: T, emmean, lower.CL, upper.CL, .group
  ROWS <- as.character(meanPairs[[1]])     # treatment labels in the same order
  diffs <- meanPairs$emmean
  ucl <- meanPairs$upper.CL
  lcl <- meanPairs$lower.CL
  letters_grp <- meanPairs$.group
  
  # compute group-wise means and SE (base R, robust)
  bwDCt <- x$wDCt
  means_by_T <- tapply(bwDCt, x$T, function(z) mean(z, na.rm = TRUE))
  sds_by_T   <- tapply(bwDCt, x$T, function(z) stats::sd(z, na.rm = TRUE))
  n_by_T     <- tapply(bwDCt, x$T, function(z) sum(!is.na(z)))
  se_by_T    <- sds_by_T / sqrt(n_by_T)
  
  se_df <- data.frame(T = names(means_by_T),
                      mean = as.numeric(means_by_T),
                      se = as.numeric(se_by_T),
                      stringsAsFactors = FALSE)
  
  # match se to the order used by emmeans/cld (ROWS)
  se_matched <- se_df$se[match(ROWS, se_df$T)]
  
  # build Results table
  Results <- data.frame(row.names = ROWS,
                        RE = 2^(-diffs),
                        log2FC = log2(2^(-diffs)),
                        LCL = 2^(-lcl),
                        UCL = 2^(-ucl),
                        se = se_matched,
                        letters = trimws(letters_grp), 
                        stringsAsFactors = FALSE)
  
  # preserve rownames as a column for splitting
  Results$RowNames <- rownames(Results)


  ## ---- split RowNames back to factor columns (base R) ----
  parts <- strsplit(Results$RowNames, ":", fixed = TRUE)
  parts_mat <- do.call(rbind, lapply(parts, function(p) {
    # ensure length matches number of factor columns
    length(p) <- length(factors)
    p
  }))
  parts_df <- as.data.frame(parts_mat, stringsAsFactors = FALSE)
  names(parts_df) <- factors
  
  # combine parts_df (factor columns) with Results
  Results_combined <- cbind(parts_df, Results)
  rownames(Results_combined) <- NULL
  
  # compute Lower/Upper SE for RE and attach to Results
  Results_combined$Lower.se.RE <- 2^(log2(Results_combined$RE) - Results_combined$se)
  Results_combined$Upper.se.RE <- 2^(log2(Results_combined$RE) + Results_combined$se)
  
  # compute log2FC SE bounds (vectorized)
  # initialize
  Results_combined$Lower.se.log2FC <- 0
  Results_combined$Upper.se.log2FC <- 0
  
  # vectorized computation replacing the for loop
  idx_less1 <- Results_combined$RE < 1
  idx_ge1   <- !idx_less1
  
  # when RE < 1
  Results_combined$Lower.se.log2FC[idx_less1] <- (Results_combined$Upper.se.RE[idx_less1] *
                                                    log2(Results_combined$RE[idx_less1])) / Results_combined$RE[idx_less1]
  Results_combined$Upper.se.log2FC[idx_less1] <- (Results_combined$Lower.se.RE[idx_less1] *
                                                    log2(Results_combined$RE[idx_less1])) / Results_combined$RE[idx_less1]
  
  # when RE >= 1
  Results_combined$Lower.se.log2FC[idx_ge1] <- (Results_combined$Lower.se.RE[idx_ge1] *
                                                  log2(Results_combined$RE[idx_ge1])) / Results_combined$RE[idx_ge1]
  Results_combined$Upper.se.log2FC[idx_ge1] <- (Results_combined$Upper.se.RE[idx_ge1] *
                                                  log2(Results_combined$RE[idx_ge1])) / Results_combined$RE[idx_ge1]
  
  # reorder columns: put lower/upper SEs before letters
  # find current column positions
  # We'll place Lower.se.RE, Upper.se.RE, Lower.se.log2FC, Upper.se.log2FC before 'letters'
  cols <- colnames(Results_combined)
  letters_pos <- match("letters", cols)
  new_order <- c(cols[1:(letters_pos - 1)],
                 "Lower.se.RE", "Upper.se.RE", "Lower.se.log2FC", "Upper.se.log2FC",
                 cols[letters_pos:length(cols)])
  # deduplicate and keep only existing columns
  new_order <- unique(new_order[new_order %in% cols])
  Results_final <- Results_combined[, new_order, drop = FALSE]
  Results_final$RowNames <- NULL
  
  # prepare output
  # remove the temporary column T from original data for output (you created x$T earlier)
  xx <- x[, setdiff(names(x), "T"), drop = FALSE]
  
  Results_final <- Results_final %>%
    mutate_if(is.numeric, ~ round(., 4))
  
  outlist2 <- structure(list(Final_data = xx,
                             lmCRD = lm_fit,
                             ANOVA = anovaCRD,
                             Results = Results_final),
                        class = "XX")
  
  print.XX <- function(outlist2) {
    print(outlist2$ANOVA)
    cat("\n", sep = '', "Relative expression (DCt method)", "\n")
    print(outlist2$Results)
    invisible(outlist2)
  }
  
  print.XX(outlist2)
}
