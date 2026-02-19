#' Delta Ct ANOVA analysis with optional model specification
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
#' @param set_missing_target_Ct_to_40 If \code{TRUE}, missing target gene Ct values become 40; if \code{FALSE} (default), they become NA. 
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
    model = NULL,
    modelBased_se = TRUE,
    set_missing_target_Ct_to_40 = FALSE
) {
  
  default_model_formula <- NULL
  
  ## Model specification
  if (!is.null(model)) {
    if (!inherits(model, "formula")) model <- as.formula(model)
    message("Using user defined formula. Ignoring block and numOfFactors for model specification.")
  } else {
    factors <- colnames(x)[1:numOfFactors]
    rhs <- paste(factors, collapse = " * ")
    default_model_formula <- if (is.null(block)) {
      paste("wDCt ~", rhs)
    } else {
      paste("wDCt ~", block, "+", rhs)
    }
  }
  
  n <- ncol(x)
  nDesign <- numOfFactors + if (is.null(block)) 1 else 2
  if (nDesign >= n) stop("Not enough columns for target and reference genes")
  
  designCols <- seq_len(nDesign)
  refCols <- (n - 2 * numberOfrefGenes + 1):n
  targetCols <- setdiff(seq_len(n), c(designCols, refCols))
  
  if (length(targetCols) == 0 || length(targetCols) %% 2 != 0) {
    stop("Target genes must be supplied as E/Ct column pairs")
  }
  
  targetPairs <- split(targetCols, ceiling(seq_along(targetCols) / 2))
  targetNames <- vapply(targetPairs, function(tc) colnames(x)[tc[1]], character(1))
  
  if (!isTRUE(analyseAllTarget)) {
    keep <- targetNames %in% analyseAllTarget
    if (!any(keep)) stop("None of the specified target genes were found in the data.")
    targetPairs <- targetPairs[keep]
    targetNames <- targetNames[keep]
  }
  
  perGene <- lapply(seq_along(targetPairs), function(i) {
    
    gene_name <- targetNames[i]
    gene_df <- x[, c(designCols, targetPairs[[i]], refCols), drop = FALSE]
    
    gene_df <- compute_wDCt(
      gene_df,
      numOfFactors,
      numberOfrefGenes,
      block,
      set_missing_target_Ct_to_40 = set_missing_target_Ct_to_40
    )
    
    gene_df[] <- lapply(gene_df, function(z)
      if (is.factor(z)) as.character(z) else z)
    
    nFac <- if (is.null(block)) numOfFactors else numOfFactors + 1
    gene_df[seq_len(nFac)] <- lapply(
      gene_df[seq_len(nFac)],
      function(col) factor(col, levels = unique(col))
    )
    
    factors <- colnames(gene_df)[1:numOfFactors]
    
    ## Treatment ID (single source of truth)
    gene_df$T <- do.call(paste, c(gene_df[factors], sep = ":"))
    
    ## Model
    if (is.null(model)) {
      rhs <- paste(factors, collapse = " * ")
      model_i <- if (is.null(block)) {
        as.formula(paste("wDCt ~", rhs))
      } else {
        as.formula(paste("wDCt ~", block, "+", rhs))
      }
    } else {
      model_i <- model
    }
    
    #has_random_effects <- grepl("\\|", as.character(model_i)[3])
    has_random_effects <- grepl("\\|", paste(deparse(model_i), collapse = " "))
    
    is_singular <- FALSE
    
    lm <- if (has_random_effects) {
      fit <- suppressMessages(lmerTest::lmer(model_i, data = gene_df, na.action = na.exclude))
      is_singular <- lme4::isSingular(fit)
      fit
    } else {
      stats::lm(model_i, data = gene_df, na.action = na.exclude)
    }
    
    ANOVA_table <- stats::anova(lm)
    lm_formula <- formula(lm)
    
    ## Estimated marginal means
    emm_formula <- as.formula(
      paste("pairwise ~", paste(factors, collapse = " * "))
    )
    
    emm_obj <- suppressMessages(
      emmeans::emmeans(
        lm,
        specs = emm_formula,
        adjust = p.adj,
        mode = "satterthwaite"
      )
    )[[1]]
    
    emm_df <- as.data.frame(emm_obj)
    
    
    
    
    if (!modelBased_se) {
      wDCt <- gene_df$wDCt
    } else {
      wDCt <- residuals(lm, type = "response")
    }
    
    
    ## Raw means and se
    obs_df <- stats::aggregate(
      wDCt ~ T,
      data = gene_df,
      FUN = function(z) c(
        mean = mean(z, na.rm = TRUE),
        se   = stats::sd(z, na.rm = TRUE) / sqrt(sum(!is.na(z)))
      )
    )
    
    obs_df <- do.call(data.frame, obs_df)
    colnames(obs_df) <- c("T", "dCt", "se")
    
    ## Treatment ID for emmeans
    emm_df$T <- do.call(paste, c(emm_df[factors], sep = ":"))
    
    merged <- merge(
      emm_df,
      obs_df,
      by = "T",
      all.x = TRUE,
      sort = FALSE
    )
    
    ## Fold change
    RE <- 2^(-merged$dCt)
    log2FC <- log2(RE)
    RE_LCL <- 2^(-merged$upper.CL)
    RE_UCL <- 2^(-merged$lower.CL)
    
    ## Compact letter display
    cld_df <- as.data.frame(
      multcomp::cld(
        emm_obj,
        adjust = p.adj,
        alpha = alpha,
        reversed = FALSE,
        Letters = letters
      )
    )
    
    cld_df$T <- do.call(paste, c(cld_df[factors], sep = ":"))
    cld_df$.group <- trimws(cld_df$.group)
    
    merged <- merge(merged,
                    cld_df[, c("T", ".group")],
                    by = "T",
                    all.x = TRUE,
                    sort = FALSE)
    
    merged$sig <- merged$.group
    merged$.group <- NULL
    
    
    
    Results <- merged[, c("T", factors, "dCt")]
    Results$RE <- RE
    Results$log2FC <- log2FC
    Results$LCL <- RE_LCL
    Results$UCL <- RE_UCL
    Results$se <- merged$se
    Results$sig <- merged$sig
    
    ## SE propagation
    Results$Lower.se.RE <- 2^(log2(Results$RE) - Results$se)
    Results$Upper.se.RE <- 2^(log2(Results$RE) + Results$se)
    
    idx_less1 <- Results$RE < 1
    idx_ge1   <- !idx_less1
    
    Results$Lower.se.log2FC <- Results$Upper.se.log2FC <- 0
    
    Results$Lower.se.log2FC[idx_less1] <-
      (Results$Upper.se.RE[idx_less1] * log2(Results$RE[idx_less1])) /
      Results$RE[idx_less1]
    
    Results$Upper.se.log2FC[idx_less1] <-
      (Results$Lower.se.RE[idx_less1] * log2(Results$RE[idx_less1])) /
      Results$RE[idx_less1]
    
    Results$Lower.se.log2FC[idx_ge1] <-
      (Results$Lower.se.RE[idx_ge1] * log2(Results$RE[idx_ge1])) /
      Results$RE[idx_ge1]
    
    Results$Upper.se.log2FC[idx_ge1] <-
      (Results$Upper.se.RE[idx_ge1] * log2(Results$RE[idx_ge1])) /
      Results$RE[idx_ge1]
    
    Results <- Results %>%
      dplyr::mutate_if(is.numeric, ~ round(., 5))
    
    Results$gene <- gene_name
    
    list(
      Final_data = gene_df[, -ncol(gene_df)],
      lm = lm,
      lm_formula = lm_formula,
      ANOVA_table = ANOVA_table,
      Results = Results,
      is_singular = is_singular
    )
  })
  
  relativeExpression <- do.call(rbind, lapply(perGene, `[[`, "Results"))
  rownames(relativeExpression) <- NULL
  
  
  relativeExpression <- relativeExpression[, c("gene",
                                               setdiff(colnames(relativeExpression), c("gene", "T", "sig")),
                                               "sig"), drop = FALSE]
  
  cat("\nRelative Expression\n\n")
  print(relativeExpression)
  
  
  singular_vec <- vapply(perGene, `[[`, logical(1), "is_singular")
  singular_genes <- targetNames[singular_vec]
  
  if (any(singular_vec)) {
    warning(
      "Singular fit detected for the following genes:\n  ",
      paste(singular_genes, collapse = ", ")
    )
  }
  
  if (is.null(model) && !is.null(default_model_formula)) {
    cat("\nNote: Using default model for statistical analysis:",
        default_model_formula, "\n")
  }
  
  invisible(list(
    perGene = setNames(perGene, targetNames),
    relativeExpression = relativeExpression,
    singular_genes = singular_genes,
    default_model_formula = default_model_formula
  ))
}
