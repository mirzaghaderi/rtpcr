#' Delta Ct ANOVA analysis with flexible model specification
#'
#' Performs Delta Ct (dCt) analysis of the data from a 1-, 2-, or 3-factor experiment 
#' with support for both fixed effects and mixed effects models. Per-gene statistical 
#' grouping is performed for all treatment combinations.
#'
#' @details
#' The function performs ANOVA analysis on weighted delta Ct (wDCt) values and returns 
#' variance components along with an expression table containing:
#' \itemize{
#' \item \code{gene}: Name of target genes
#' \item Factor columns: Experimental design factors
#' \item \code{dCt}: Mean weighted delta Ct for each treatment combination
#' \item \code{RE}: Relative expression = 2^-dCt
#' \item \code{log2FC}: log2 of relative expression
#' \item \code{LCL}: 95\% lower confidence level
#' \item \code{UCL}: 95\% upper confidence level
#' \item \code{se}: Standard error of the mean calculated from wDCt values
#' \item \code{Lower.se.RE}: Lower limit error bar for RE (2^(log2(RE) - se))
#' \item \code{Upper.se.RE}: Upper limit error bar for RE (2^(log2(RE) + se))
#' \item \code{Lower.se.log2FC}: Lower limit error bar for log2 RE
#' \item \code{Upper.se.log2FC}: Upper limit error bar for log2 RE
#' \item \code{sig}: Per-gene significance grouping letters
#' }
#'
#'
#' @param x The input data frame containing experimental design columns, target gene
#'   E/Ct column pairs, and reference gene E/Ct column pairs. Reference gene
#'   columns must be located at the end of the data frame. See "Input data 
#'   structure" in vignettes for details about data structure.
#' @param numOfFactors Integer. Number of experimental factor columns 
#'   (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes. Each reference gene
#'   must be represented by two columns (E and Ct).
#' @param block Character. Block column name or \code{NULL}. When a qPCR experiment 
#'   is done in multiple qPCR plates, variation resulting from the plates may 
#'   interfere with the actual amount of gene expression. One solution is to 
#'   conduct each plate as a randomized block so that at least one replicate of 
#'   each treatment and control is present on a plate. Block effect is usually 
#'   considered as random and its interaction with any main effect is not considered.
#'   Note: This parameter is ignored if \code{model} is provided.
#' @param alpha Statistical level for comparisons (default: 0.05).
#' @param p.adj Method for p-value adjustment. See \code{\link[stats]{p.adjust}}.
#' @param analyseAllTarget Logical or character. If \code{TRUE} (default), all 
#'   detected target genes are analysed. Alternatively, a character vector 
#'   specifying the names (names of their Efficiency columns) of target genes 
#'   to be analysed.
#' @param model Optional model formula. If provided, this overrides the automatic formula (CRD or RCBD 
#'   based on \code{block} and \code{numOfFactors}). The formula uses 
#'   \code{wDCt} as the response variable. 
#'   For mixed models, random effects can be defined using \code{lmer} syntax 
#'   (e.g., \code{"wDCt ~ Treatment + (1|Block)"}). When using \code{model}, 
#'   the \code{block} and \code{numOfFactors} arguments are ignored for model 
#'   specification, but still used for data structure identification.
#'   
#' for fixed effects only, the \code{"lm"} (ordinary least squares) is used. 
#' \code{"lmer"} is used for mixed effects models 
#' (requires the \code{lmerTest} package). If a custom formula is provided with 
#' random effects, the function will use \code{lmerTest::lmer()}; otherwise 
#' it will use \code{stats::lm()}. Note that \code{emmeans} supports both 
#' model types and will use appropriate degrees of freedom methods (Satterthwaite by default).
#'
#' @return
#' An object containing expression tables, lm/lmer models, ANOVA tables, 
#' residuals, and raw data for each gene:
#' \describe{
#' \item{\code{relativeExpression}}{dCt expression table for all treatment 
#'   combinations along with per-gene statistical grouping}
#' \item{\code{perGene}}{Nested list containing detailed results for each target gene:
#'   \itemize{
#'   \item \code{ANOVA_table}: Full factorial ANOVA table
#'   \item \code{lm}: lm/lmer model for factorial design
#'   \item \code{Final_data}: Processed data with wDCt values
#'   \item \code{resid(object$perGene$gene_name$lm)}: Residuals
#'   }}
#' }
#' 
#' @examples
#' # Default usage with fixed effects
#' result <- ANOVA_DCt(data_2factorBlock3ref, numOfFactors = 2, numberOfrefGenes = 3, 
#'                     block = "block")
#'
#' # Mixed model with random block effect
#' result_mixed <- ANOVA_DCt(data_2factorBlock3ref, numOfFactors = 2, numberOfrefGenes = 3,
#'                           block = "block")
#'
#' # Custom mixed model formula with nested random effects
#' result_custom <- ANOVA_DCt(data_repeated_measure_2, numOfFactors = 2, numberOfrefGenes = 1,
#'                             block = NULL,
#'                             model = wDCt ~ treatment * time + (1 | id))
#'
#'
#' @export


ANOVA_DCt <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    block = NULL,
    alpha = 0.05,
    p.adj = "none",
    analyseAllTarget = TRUE,
    model = NULL
) {

  
  # Store default model formula for messaging
  default_model_formula <- NULL
  
  # If model is provided, ignore other formula-related arguments
  if (!is.null(model)) {
    if (!inherits(model, "formula")) {
      model <- as.formula(model)
    }
    message("Using custom formula. Ignoring block and numOfFactors for model specification.")
  } else {
    # Create default model formula for message
    factors <- colnames(x)[1:numOfFactors]
    if (is.null(block)) {
      default_model_formula <- paste("wDCt ~", paste(factors, collapse = " * "))
    } else {
      default_model_formula <- paste("wDCt ~", block, "+", paste(factors, collapse = " * "))
    }
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
      p.adj = p.adj,
      model = model
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
  
  # Add default model message if no user-defined model was provided
  if (is.null(model) && !is.null(default_model_formula)) {
    cat("\nNote: Using default model for analysis:", default_model_formula, "\n")
  }
  
  # Return full structured object invisibly
  invisible(list(
    perGene = setNames(perGene, targetNames),  # All individual gene outputs with models
    relativeExpression = relativeExpression,    # Combined table
    default_model_formula = default_model_formula  # Store default formula in output
  ))
}

# Internal helper function
.ANOVA_DCt_uniTarget <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    block,
    alpha = 0.05,
    p.adj = "none",
    model = NULL,
    verbose = FALSE
) {
  
  # basic argument checks
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame")
  }
  if (!is.numeric(numberOfrefGenes) || length(numberOfrefGenes) != 1) {
    stop("`numberOfrefGenes` must be a single numeric value")
  }
  
  # Compute wDCt values (these are what we'll compare statistically)
  x <- compute_wDCt(x, numOfFactors, numberOfrefGenes, block)
  
  # Convert all factor columns to character
  x[] <- lapply(x, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  # Convert the first numOfFactors columns to factor with proper levels
  if (is.null(block)) {
    x[seq_len(numOfFactors)] <- lapply(
      x[seq_len(numOfFactors)],
      function(col) factor(col, levels = unique(col)))
  } else {
    x[seq_len(numOfFactors + 1)] <- lapply(
      x[seq_len(numOfFactors + 1)],
      function(col) factor(col, levels = unique(col)))
  }
  
  # Get names of factor columns
  factors <- colnames(x)[1:numOfFactors]
  
  # Create T factor for grouping (treatment combinations)
  x$T <- do.call(paste, c(x[1:length(factors)], sep = ":"))
  x$T <- as.factor(x$T)
  
  # Handle custom formula or generate default
  if (!is.null(model)) {
    model <- model
  } else {
    if (is.null(block)) {
      model <- as.formula(paste("wDCt ~", paste(factors, collapse = " * ")))
    } else {
      model <- as.formula(paste("wDCt ~", block, "+", paste(factors, collapse = " * ")))
    }
  }
  
  # detect model_type based on formula
  formula_text <- as.character(model)[3]  # Get RHS of formula
  has_random_effects <- grepl("\\|", formula_text)
  
  if (has_random_effects) {
    lm <- lmerTest::lmer(model, data = x)    # suppressMessages
  } else {
    lm <- stats::lm(model, data = x)
  }
  
  # Get ANOVA table
  ANOVA_table <- stats::anova(lm)
  lm_formula <- formula(lm)
  
  # Create emmeans specification
  ABC <- as.formula(paste("pairwise ~", paste(factors, collapse = " * ")))
  
  # Get emmeans from the model for statistical comparisons
  emg <- suppressMessages(
    emmeans::emmeans(lm, specs = ABC,
                     adjust = p.adj, mode = "satterthwaite")
  )
  
  emm_obj <- emg[[1]]
  emm_df <- as.data.frame(emm_obj)
  # Get treatment levels in the order they appear in emmeans
  ROWS <- do.call(paste, c(emm_df[factors], sep = ":"))
  
  # Use cld() to get statistical groupings
  meanPairs <- multcomp::cld(emm_obj, adjust = p.adj, alpha = alpha,
                             reversed = FALSE, Letters = letters)
  meanPairs_df <- as.data.frame(meanPairs)
  
  # Extract model estimates and confidence intervals
  model_estimates <- meanPairs_df$emmean
  model_lcl <- meanPairs_df$lower.CL
  model_ucl <- meanPairs_df$upper.CL
  letters_grp <- trimws(meanPairs_df$.group)
  
  # Calculate OBSERVED wDCt values
  bwDCt <- x$wDCt
  T <- do.call(paste, c(x[factors], sep = ":"))
  observed_means <- tapply(bwDCt, x$T, function(z) mean(z, na.rm = TRUE))
  observed_sds <- tapply(bwDCt, x$T, function(z) stats::sd(z, na.rm = TRUE))
  observed_n <- tapply(bwDCt, x$T, function(z) sum(!is.na(z)))
  observed_se <- observed_sds / sqrt(observed_n)
  
  # Match observed values to the order used by emmeans (ROWS)
  dCt_observed <- observed_means[match(ROWS, names(observed_means))]
  se_observed <- observed_se[match(ROWS, names(observed_se))]
  
  # Calculate Relative Expression (RE) from wDCt
  RE_observed <- 2^(-dCt_observed)
  log2FC_observed <- log2(RE_observed)
  
  # Transform model confidence intervals for RE
  RE_LCL <- 2^(-model_ucl)
  RE_UCL <- 2^(-model_lcl)
  
  # build Results table
  Results <- data.frame(row.names = ROWS,
                        dCt = dCt_observed,
                        RE = RE_observed,
                        log2FC = log2FC_observed,
                        LCL = RE_LCL,
                        UCL = RE_UCL,
                        se = se_observed,
                        stringsAsFactors = FALSE)
  
  Results <- Results[order(Results[[1]]), ]
  Results <- data.frame(Results, sig = letters_grp)
  
  # Add factor columns by splitting T
  Results$RowNames <- rownames(Results)
  
  # split RowNames back to factor columns
  parts <- strsplit(Results$RowNames, ":", fixed = TRUE)
  parts_mat <- do.call(rbind, lapply(parts, function(p) {
    length(p) <- length(factors)
    p
  }))
  parts_df <- as.data.frame(parts_mat, stringsAsFactors = FALSE)
  names(parts_df) <- factors
  
  # combine parts_df with Results
  Results_combined <- cbind(parts_df, Results)
  rownames(Results_combined) <- NULL
  
  # compute Lower/Upper SE for RE and log2FC (based on observed SE)
  Results_combined$Lower.se.RE <- 2^(log2(Results_combined$RE) - Results_combined$se)
  Results_combined$Upper.se.RE <- 2^(log2(Results_combined$RE) + Results_combined$se)
  
  Results_combined$Lower.se.log2FC <- 0
  Results_combined$Upper.se.log2FC <- 0
  
  idx_less1 <- Results_combined$RE < 1
  idx_ge1   <- !idx_less1
  
  # Calculate SE bounds for log2FC
  Results_combined$Lower.se.log2FC[idx_less1] <- (Results_combined$Upper.se.RE[idx_less1] *
                                                    log2(Results_combined$RE[idx_less1])) / Results_combined$RE[idx_less1]
  Results_combined$Upper.se.log2FC[idx_less1] <- (Results_combined$Lower.se.RE[idx_less1] *
                                                    log2(Results_combined$RE[idx_less1])) / Results_combined$RE[idx_less1]
  
  Results_combined$Lower.se.log2FC[idx_ge1] <- (Results_combined$Lower.se.RE[idx_ge1] *
                                                  log2(Results_combined$RE[idx_ge1])) / Results_combined$RE[idx_ge1]
  Results_combined$Upper.se.log2FC[idx_ge1] <- (Results_combined$Upper.se.RE[idx_ge1] *
                                                  log2(Results_combined$RE[idx_ge1])) / Results_combined$RE[idx_ge1]
  
  # Reorder columns: put factor columns first, then results
  cols <- colnames(Results_combined)
  
  # Identify which columns are factor columns
  factor_cols <- factors
  other_cols <- setdiff(cols, factor_cols)
  
  # Create new order
  # Remove RowNames if present
  other_cols <- other_cols[other_cols != "RowNames"]
  
  # Define preferred order for result columns
  result_cols_order <- c("dCt", "RE", "log2FC", "LCL", "UCL", "se",
                         "Lower.se.RE", "Upper.se.RE",
                         "Lower.se.log2FC", "Upper.se.log2FC", "sig")
  
  # Keep only columns that exist
  result_cols_order <- result_cols_order[result_cols_order %in% other_cols]
  
  # Add any remaining columns
  remaining_cols <- setdiff(other_cols, result_cols_order)
  final_cols_order <- c(factor_cols, result_cols_order, remaining_cols)
  
  Results_final <- Results_combined[, final_cols_order, drop = FALSE]
  
  # remove the temporary column T from original data
  xx <- x[, setdiff(names(x), "T"), drop = FALSE]
  
  # Round numeric columns
  Results_final <- Results_final %>%
    dplyr::mutate_if(is.numeric, ~ round(., 5))
  
  # Return output list
  list(Final_data = xx,
       lm = lm,
       lm_formula = lm_formula,
       ANOVA_table = ANOVA_table,
       Results = Results_final)
}