#' Delta Delta Ct ANOVA analysis with optional model specification
#'
#' Apply Delta Delta Ct (ddCt) analysis to each target gene
#' and performs per-gene statistical analysis.
#'
#' @param x The input data frame containing experimental design columns, replicates (integer), target gene
#' E/Ct column pairs, and reference gene E/Ct column pairs. Reference gene
#' columns must be located at the right end of the data frame. See "Input data 
#' structure" in vignettes for details about data structure.
#' @param numOfFactors Integer. Number of experimental factor columns
#' (excluding \code{rep} and optional \code{block}).
#' @param numberOfrefGenes Integer. Number of reference genes.
#' @param block Character. Block column name or \code{NULL}. 
#' When a qPCR experiment is done in multiple qPCR plates, 
#' variation resulting from the plates may interfere with the actual amount of 
#' gene expression. One solution is to conduct each plate as a randomized block 
#' so that at least one replicate of each treatment and control is present 
#' on a plate. Block effect is usually considered as random and its interaction 
#' with any main effect is not considered.
#' @param specs Example: "A", "A|B" or A|B*C" if A, B and C are name of factor columns in the input data
#' The first name (here A) is the factor for which the relative expression analysis is applied.
#' @param calibratorLevel NULL or one of the levels of the first selected factor in specs argument. If NULL the first level of that factor is used as calibrator.
#' Optional character vector specifying the order of levels for the main factor.
#' If \code{NULL}, the first observed level is used as the calibrator.
#' If provided, the first element of the vector is used as the calibrator level.
#' @param analyseAllTarget Logical or character.
#' If \code{TRUE} (default), all target genes are analysed.
#' Alternatively, a character vector specifying the names (names of their Efficiency columns) of target genes
#' to be analysed.
#' @param model Optional model formula. If provided, this overrides the automatic formula (uni - or multi-factorial CRD or RCBD 
#' based on \code{block} and \code{numOfFactors}). The formula uses 
#' \code{wDCt} as the response variable. 
#' For mixed models, random effects can be defined using \code{lmer} syntax 
#' (e.g., \code{"wDCt ~ Treatment + (1 | id)"}). When using \code{model}, 
#' the \code{block} and \code{numOfFactors} arguments are ignored for model 
#' specification, but still used for data structure identification.
#'   
#' for fixed effects only, the \code{"lm"} (ordinary least squares) is used. 
#' \code{"lmer"} is used for mixed effects models 
#' (requires the \code{lmerTest} package). If a custom formula is provided with 
#' random effects, the function will use \code{lmerTest::lmer()}; otherwise 
#' it will use \code{stats::lm()}. Note that \code{emmeans} supports both 
#' model types and will use appropriate degrees of freedom methods (Satterthwaite by default).
#' @param p.adj
#' Method for p-value adjustment. See \code{\link[stats]{p.adjust}}.
#' @param set_missing_target_Ct_to_40 If \code{TRUE}, missing target gene Ct values become 40; if \code{FALSE} (default), they become NA. 
#' @param se.type Character string specifying how standard error is calculated. 
#' One of \code{"paired.group"}, \code{"two.group"}, or \code{"single.group"}. 
#' \code{"paired.group"} computes SE from paired differences (used when a random 
#' \code{id} effect is present), \code{"two.group"} uses the unpaired two-group 
#' t-test standard error against the reference level, and \code{"single.group"} 
#' computes SE within each level using a one-group t-test.
#' @param modelBased_se Logical. If \code{TRUE} (default), standard errors are  
#' calculated from model-based residuals. If \code{FALSE}, standard errors are calculated directly from the observed 
#' \code{wDCt} values within each treatment group according to the selected \code{se.type}.  
#' For single factor data, both methods are the same. It is recommended to use \code{modelBased_se = TRUE} (default).
#' @param ... Additional arguments. Included for backward compatibility with deprecated \code{mainFactor.column}.
#' 
#' @importFrom stats setNames
#'
#' @details
#' ddCt analysis of variance (ANOVA) is performed for 
#' the \code{mainFactor.column} based on a full model factorial 
#' experiment by default. However, if \code{ANCOVA_DDCt} function is used, 
#' analysis of covariance is performed for the levels of the \code{mainFactor.column} and the other factors are 
#' treated as covariates. if the interaction between the main factor and the covariate is significant, ANCOVA is not appropriate.
#' 
#' All the functions for relative expression analysis (including \code{TTEST_DDCt()}, 
#' \code{WILCOX_DDCt()}, \code{ANOVA_DDCt()}, and \code{ANOVA_DCt()}) return the 
#' relative expression table which include fold change and corresponding 
#' statistics. The output of \code{ANOVA_DDCt()}, 
#' and \code{ANOVA_DCt()} also include lm models, residuals, raw data and ANOVA table 
#' for each gene. 
#' 
#' The expression table returned by \code{TTEST_DDCt()}, 
#' \code{WILCOX_DDCt()}, and \code{ANOVA_DDCt()} functions 
#' include these columns: gene (name of target genes), 
#' contrast (calibrator level and contrasts for which the relative expression is computed), 
#' ddCt (mean of weighted delta delta Ct values), RE (relative expression or 
#' fold change = 2^-ddCt),  log2FC (log(2) of relative expression or fold change), 
#' pvalue, sig (per-gene significance), LCL (95\% lower confidence level), UCL (95\% upper confidence level),
#' se (standard error of mean calculated from the weighted delta Ct values of each of the main factor levels),
#' Lower.se.RE (The lower limit error bar for RE which is 2^(log2(RE) - se)), 
#' Upper.se.RE (The upper limit error bar for RE which is 2^(log2(RE) + se)),
#' Lower.se.log2FC (The lower limit error bar for log2 RE), and 
#' Upper.se.log2FC (The upper limit error bar for log2 RE)
#'
#' @import emmeans
#' @return
#' An object containing expression table, lm model, residuals, raw data and ANOVA table for each gene:
#' \describe{  
#' \item{ddCt expression table along with per-gene statistical comparison outputs}{\code{object$relativeExpression}}
#' \item{ANOVA table}{\code{object$perGene$gene_name$ANOVA_table}}
#' \item{lm ANOVA}{\code{object$perGene$gene_name$lm}}
#' \item{lm_formula}{\code{object$perGene$gene_name$lm_formula}}
#' \item{Residuals}{\code{resid(object$perGene$gene_name$lm)}}
#' }
#' @export
#' 
#' @references
#' LivakKJ, Schmittgen TD (2001).
#' Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR
#' and the Double Delta CT Method.
#' \emph{Methods}, 25(4), 402–408.
#' doi:10.1006/meth.2001.1262
#'
#' Ganger MT, Dietz GD, and Ewing SJ (2017).
#' A common base method for analysis of qPCR data and the application of
#' simple blocking in qPCR experiments.
#' \emph{BMC Bioinformatics}, 18, 1–11.
#' 
#' Taylor SC, Nadeau K, Abbasi M, Lachance C, Nguyen M, Fenrich, J. (2019). 
#' The ultimate qPCR experiment: producing publication quality, reproducible 
#' data the first time. \emph{Trends in Biotechnology}, 37, 761-774. 
#' 
#' Yuan JS, Reed A, Chen F, Stewart N (2006).
#' Statistical Analysis of Real-Time PCR Data.
#' \emph{BMC Bioinformatics}, 7, 85.
#'
#' 
#' @examples
#' data1 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
#' ANOVA_DDCt(data1,
#'            numOfFactors = 2,
#'            numberOfrefGenes = 3,
#'            block = "block",
#'            specs = "Concentration",
#'            p.adj = "none")
#'            
#' data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
#' ANOVA_DDCt(data2,
#'            numOfFactors = 1,
#'            numberOfrefGenes = 1,
#'            block = NULL,
#'            specs = "Condition",
#'            p.adj = "none",
#'            se.type = "single.group")
#'   
#' # Repeated measure analysis
#' a <- ANOVA_DDCt(data_repeated_measure_1,
#'            numOfFactors = 1,
#'            numberOfrefGenes = 1,
#'            block = NULL,
#'            specs = "time",
#'            p.adj = "none", model = wDCt ~ time + (1 | id))
#'
#' a$perGene$Target$ANOVA_table
#'
#'
#' # Repeated measure analysis: split-plot in time
#' a <- ANOVA_DDCt(data_repeated_measure_2,
#'            numOfFactors = 2, numberOfrefGenes = 1,
#'            specs = "time", block = NULL,
#'            model = wDCt ~ treatment * time + (1 | id))
#'    
ANOVA_DDCt <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    specs,           
    block,           
    calibratorLevel = NULL, 
    p.adj = "none",
    analyseAllTarget = TRUE,
    model = NULL,
    set_missing_target_Ct_to_40 = FALSE,
    se.type = c("single.group", "paired.group", "two.group"),
    modelBased_se = TRUE, ...) {
  
  se.type <- match.arg(se.type)

  
  dots <- list(...)
  deprecated_args <- c("mainFactor.column", "mainFactor.level.order")
  for (arg in deprecated_args) {
    if (arg %in% names(dots)) {
      lifecycle::deprecate_warn(
        when = "2.1.5", 
        what = paste0("ANOVA_DDCt(", arg, ")"), 
        with = "ANOVA_DDCt(specs)") } }
  
  
  n <- ncol(x)
  nDesign <- if (is.null(block)) numOfFactors + 1 else numOfFactors + 2
  designCols <- seq_len(nDesign)
  nRefCols <- 2 * numberOfrefGenes
  refCols <- (n - nRefCols + 1):n
  targetCols <- setdiff(seq_len(n), c(designCols, refCols))
  
  if (length(targetCols) == 0 || length(targetCols) %% 2 != 0) {
    stop("The supplied arguments don't match the input data.")
  }
  
  
  detect_rep_id_random <- function(model, numOfFactors, block, x) {
    if (is.null(model)) return(FALSE)
    ftxt <- paste(deparse(formula(model)), collapse = " ")
    rep_id_col <- if (is.null(block)) numOfFactors + 1 else numOfFactors + 2
    rep_id_name <- colnames(x)[rep_id_col]
    grepl(paste0("\\|\\s*", rep_id_name, "\\b"), ftxt)
  }
  
  
  
  targetPairs <- split(targetCols, ceiling(seq_along(targetCols)/2))
  targetNames <- vapply(targetPairs, function(tc) colnames(x)[tc[1]], character(1))
  
  if (!isTRUE(analyseAllTarget)) {
    keep <- targetNames %in% analyseAllTarget
    if (!any(keep)) stop("None of the specified target genes were found in the data.")
    targetPairs <- targetPairs[keep]
    targetNames <- targetNames[keep]
  }
  
  perGene <- list()
  relativeExpression_list <- list()
  
  for (i in seq_along(targetPairs)) {
    
    tc <- targetPairs[[i]]
    gene_name <- targetNames[i]
    
    gene_df <- x[, c(designCols, tc, refCols), drop = FALSE]
    user_defined_model <- !is.null(model)
    
    
    
    if (is.null(block)) {
      gene_df[seq_len(numOfFactors)] <- lapply(
        gene_df[seq_len(numOfFactors)],
        function(col) factor(col, levels = unique(col)))
    } else {
      gene_df[seq_len(numOfFactors + 1)] <- lapply(
        gene_df[seq_len(numOfFactors + 1)],
        function(col) factor(col, levels = unique(col)))
    }
    
    mainFactor.level.order <- unique(gene_df[,1])
    mainFactor <- trimws(strsplit(specs, "\\|")[[1]][1])
    mainFactor_col_number <- which(names(x) == mainFactor)
    gene_df <- gene_df[, c(mainFactor_col_number, (1:ncol(gene_df))[-mainFactor_col_number])]
    if(is.null(calibratorLevel)){gene_df} else {
      gene_df <- gene_df %>% arrange(gene_df[[1]] != calibratorLevel, gene_df[[1]])
    }
    
    
    gene_df <- compute_wDCt(gene_df, numOfFactors, numberOfrefGenes, block,
                            set_missing_target_Ct_to_40 = set_missing_target_Ct_to_40)
    
    gene_df[] <- lapply(gene_df, function(x) {
      if (is.factor(x)) as.character(x) else x
    })
    
    if (is.null(block)) {
      gene_df[seq_len(numOfFactors)] <- lapply(
        gene_df[seq_len(numOfFactors)],
        function(col) factor(col, levels = unique(col)))
    } else {
      gene_df[seq_len(numOfFactors + 1)] <- lapply(
        gene_df[seq_len(numOfFactors + 1)],
        function(col) factor(col, levels = unique(col)))
    }
    
    factors <- colnames(gene_df)[1:numOfFactors]
    default_model_formula <- NULL
    is_mixed_model <- FALSE
    is_singular <- FALSE
    
    
    if (user_defined_model) {
      if (is.character(model)) {
        formula_str <- paste("wDCt ~", model)
        formula_obj <- as.formula(formula_str)
      } else if (inherits(model, "formula")) {
        formula_obj <- model
      } else {
        stop("model must be either a character string or a formula object")
      }
      
      has_random_effects <- grepl("\\|", paste(deparse(formula_obj), collapse = " "))
      
      if (has_random_effects) {
        is_mixed_model <- TRUE
        lm_fit <- suppressMessages(lmerTest::lmer(formula_obj, data = gene_df, na.action = na.exclude))
        # singularity check
        is_singular <- lme4::isSingular(lm_fit, tol = 1e-4)
      } else {
        lm_fit <- lm(formula_obj, data = gene_df, na.action = na.exclude)
      }
      
      lm_formula <- paste(deparse(formula(lm_fit)), collapse = " ")
      ANOVA_table <- stats::anova(lm_fit)
      
    } else {
      
      if (is.null(block)) {
        formula_ANOVA <- as.formula(
          paste("wDCt ~", paste(factors, collapse = " * "))
        )
        default_model_formula <- deparse(formula_ANOVA)
        lm_fit <- lm(formula_ANOVA, data = gene_df, na.action = na.exclude)
        
      } else {
        formula_ANOVA <- as.formula(
          paste("wDCt ~", block, "+", paste(factors, collapse = " * "))
        )
        default_model_formula <- deparse(formula_ANOVA)
        lm_fit <- lm(formula_ANOVA, data = gene_df, na.action = na.exclude)
      }
      
      lm_formula <- paste(deparse(formula(lm_fit)), collapse = " ")
      ANOVA_table <- stats::anova(lm_fit)
    }
    
    ############################################################################
    specs2 <- stats::as.formula(paste("~ ", specs))
    pp1 <- suppressMessages(emmeans::emmeans(lm_fit, specs2, adjust = p.adj))    # mode = "satterthwaite" data = gene_df,
    cld_results <- multcomp::cld(pp1, alpha = 0.05, Letters = letters)
    if(is.null(calibratorLevel)) {calibratorLevel <- levels(factor(gene_df[[mainFactor]]))[1]}
    stats::confint(graphics::pairs(pp1), ref = calibratorLevel)
    contrast <- emmeans::contrast(pp1, "trt.vs.ctrl", ref = calibratorLevel, type = "response", adjust = "none")
    conf <- stats::confint(contrast)
    conf$RE <- 2^(-conf$estimate)
    conf$LCL <- 2^(-conf$upper.CL)
    conf$UCL <- 2^(-conf$lower.CL)
    conf <- as.data.frame(conf)
    contrast <- as.data.frame(contrast)
    conf <- data.frame(conf, p.value = contrast$p.value)
    
    
    
    all_cols <- names(conf)
    stat_cols <- c("contrast", "estimate", "SE", "df", "lower.CL", "upper.CL", "RE", "LCL", "UCL", "p.value")
    grouping_cols <- setdiff(all_cols, stat_cols)
    conf <- conf %>%
      mutate(contrast = as.character(contrast)) %>%
      group_by(across(all_of(grouping_cols))) %>% 
      group_modify(~ {
        ref_name <- sub(".*[ -]", "", .x$contrast[1])
        ref_name <- trimws(ref_name)
        ref_row <- .x[1, ]
        ref_row$contrast <- ref_name
        ref_row$p.value  <- 1
        ref_row$RE       <- 1
        cols_to_zero <- c("estimate", "SE", "df", "lower.CL", "upper.CL", "LCL", "UCL")
        existing_numeric <- intersect(cols_to_zero, names(ref_row))
        ref_row[existing_numeric] <- 0
        bind_rows(ref_row, .x)
      }) %>%
      ungroup()
    ############################################################################
    
    
    # 1. Parse specs to separate the primary factor from conditioning factors
    # Split "A | B * C" into primary ("A") and context ("B", "C")
    spec_parts <- strsplit(specs, "\\|")[[1]]
    primary_var <- trimws(all.vars(as.formula(paste("~", spec_parts[1]))))
    context_vars <- if(length(spec_parts) > 1) {
      trimws(all.vars(as.formula(paste("~", spec_parts[2]))))
    } else {
      NULL
    }
    
    if (!modelBased_se) {
      gene_df$residual <- gene_df$wDCt
    } else {
      gene_df$residual <- residuals(lm_fit, type = "response")
    }
    
    
    idRand <- detect_rep_id_random(model = model,
                                   numOfFactors = numOfFactors,
                                   block = block, x = x)
    
    
    # Add ID column if it's a random effect/paired design
    id_col_name <- NULL
    if (idRand || se.type == "paired.group") {
      id_col_idx <- min(tc) - 1
      id_col_name <- "pair_id"
      gene_df[[id_col_name]] <- x[[id_col_idx]]
    }
    # 3. Define the calculation logic per subgroup
    calculate_se <- function(sub_df, type, primary) {
      levs <- levels(factor(sub_df[[primary]]))
      if (length(levs) == 0) return(NULL)
      ref_lev <- levs[1]
      
      if (type == "single.group") {
        return(sub_df %>%
                 group_by(!!sym(primary)) %>%
                 summarise(se = stats::sd(residual, na.rm = TRUE) / sqrt(sum(!is.na(residual))), .groups = "drop"))
      } else {
        # For two.group or paired.group
        ref_vals <- sub_df[sub_df[[primary]] == ref_lev, ]
        
        results <- lapply(levs, function(l) {
          if (l == ref_lev) return(data.frame(primary = l, se = 0))
          comp_vals <- sub_df[sub_df[[primary]] == l, ]
          
          val <- tryCatch({
            if (type == "paired.group") {
              # Strict pairing by ID
              paired <- merge(ref_vals, comp_vals, by = id_col_name)
              if(nrow(paired) < 2) return(NA_real_)
              stats::sd(paired$residual.x - paired$residual.y, na.rm = TRUE) / sqrt(nrow(paired))
            } else {
              # Independent samples t-test stderr
              stats::t.test(comp_vals$residual, ref_vals$residual, paired = FALSE)$stderr
            }
          }, error = function(e) NA_real_)
          
          data.frame(primary = l, se = val)
        })
        res_df <- do.call(rbind, results)
        colnames(res_df)[1] <- primary
        return(res_df)
      }
    }
    # 4. Execute grouped by context_vars
    if (is.null(context_vars)) {
      # Simple case: No pipe in specs
      se_final <- calculate_se(gene_df, se.type, primary_var)
    } else {
      # Nested case: Calculate within each combination of context_vars
      se_final <- gene_df %>%
        group_by(across(all_of(context_vars))) %>%
        group_modify(~ calculate_se(.x, se.type, primary_var)) %>%
        ungroup()
    }
    
    
    
    
    
    
    
    
    sig <- .convert_to_character(conf$p.value)
    post_hoc_test <- data.frame(conf,
                                log2FC = log2(1/(2^-(conf$estimate))),
                                sig = sig,
                                se = se_final$se
    )
    post_hoc_test$RE[post_hoc_test$pvalue == "NaN"] <- 0
    post_hoc_test$log2FC[post_hoc_test$pvalue == "NaN"] <- 0
    post_hoc_test$sig[post_hoc_test$pvalue == "NaN"] <- "ND"
    
    
    
    post_hoc_test$contrast <- sapply(
      strsplit(as.character(post_hoc_test$contrast), " - "),
      function(x) paste(x, collapse = " vs ")
    )
    
    post_hoc_test <- data.frame(
      post_hoc_test,
      Lower.se.RE = 2^(log2(post_hoc_test$RE) - post_hoc_test$se),
      Upper.se.RE = 2^(log2(post_hoc_test$RE) + post_hoc_test$se),
      Lower.se.log2FC = log2(post_hoc_test$RE) - post_hoc_test$se,
      Upper.se.log2FC = log2(post_hoc_test$RE) + post_hoc_test$se
    )
    
    
    post_hoc_test$gene <- gene_name
    
    res <- list(
      Final_data = gene_df,
      lm = lm_fit,
      ANOVA_table = ANOVA_table,
      Fold_Change = post_hoc_test,
      lm_formula = lm_formula,
      user_defined_model = user_defined_model,
      default_model_formula = default_model_formula,
      is_mixed_model = is_mixed_model,
      singular = is_singular
    )
    
    perGene[[gene_name]] <- res
    relativeExpression_list[[i]] <- post_hoc_test
  }
  
  relativeExpression <- do.call(rbind, relativeExpression_list)
  rownames(relativeExpression) <- NULL
  
  relativeExpression <- relativeExpression %>%
    dplyr::select(-any_of(c("SE", "df", "lower.CL", "upper.CL", "row_number..", "ddCt"))) %>%
    dplyr::select(-p.value, -sig, everything(), p.value, sig)
  names(relativeExpression)[names(relativeExpression) == "estimate"] <- "ddCt"
  relativeExpression <- relativeExpression[, c("gene", setdiff(names(relativeExpression), "gene"))]
  
  for (col in names(relativeExpression)) {
    if (is.numeric(relativeExpression[[col]])) {
      relativeExpression[[col]] <- round(relativeExpression[[col]], 5)
    }
  }
  
  
  
  
  cat("\nRelative Expression\n")
  print(relativeExpression)
  cat("\n")
  
  first_gene_res <- perGene[[1]]
  calibrator_level <- strsplit(first_gene_res$Fold_Change$contrast[1], " vs ")[[1]][1]
  cat(paste("The", calibrator_level, "level was used as calibrator.\n"))
  
  singular_genes <- names(perGene)[
    vapply(perGene, function(z) isTRUE(z$singular), logical(1))
  ]
  
  if (length(singular_genes) > 0) {
    warning(
      "Singular fit detected for the following genes:\n  ",
      paste(singular_genes, collapse = ", ")
    )
  }
  
  if (is.null(model) && length(perGene) > 0) {
    default_formula <- perGene[[1]]$default_model_formula
    if (!is.null(default_formula)) {
      cat("Note: Using default model for statistical analysis:", default_formula, "\n")
    }
  }
  
  
  invisible(list(
    perGene = perGene,
    relativeExpression = relativeExpression
  ))
}
