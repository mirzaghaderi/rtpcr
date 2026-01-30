#' Delta Delta Ct ANOVA analysis
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
#' @param mainFactor.column
#' Integer. Column index of the factor for which the relative expression analysis is applied.
#' @param mainFactor.level.order
#' Optional character vector specifying the order of levels for the main factor.
#' If \code{NULL}, the first observed level is used as the calibrator.
#' If provided, the first element of the vector is used as the calibrator level.
#' @param analyseAllTarget Logical or character.
#' If \code{TRUE} (default), all target genes are analysed.
#' Alternatively, a character vector specifying the names (names of their Efficiency columns) of target genes
#' to be analysed.
#' @param model Optional model formula. If provided, this overrides the automatic formula (CRD or RCBD 
#' based on \code{block} and \code{numOfFactors}). The formula uses 
#' \code{wDCt} as the response variable. 
#' For mixed models, random effects can be defined using \code{lmer} syntax 
#' (e.g., \code{"wDCt ~ Treatment + (1|Block)"}). When using \code{model}, 
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
#' All the functions for relative expression analysis (including `TTEST_DDCt()`, 
#' `WILCOX_DDCt()`, `ANOVA_DDCt()`, and `ANOVA_DCt()`) return the 
#' relative expression table which include fold change and corresponding 
#' statistics. The output of `ANOVA_DDCt()`, 
#' and `ANOVA_DCt()` also include lm models, residuals, raw data and ANOVA table 
#' for each gene. 
#' 
#' The expression table returned by `TTEST_DDCt()`, 
#' `WILCOX_DDCt()`, and `ANOVA_DDCt()` functions 
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
#' ANOVA_DDCt(x = data1,
#'            numOfFactors = 2,
#'            numberOfrefGenes = 3,
#'            block = "block",
#'            mainFactor.column = 2,
#'            p.adj = "none")
#'            
#' data2 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))          
#' ANOVA_DDCt(x = data2,
#'            numOfFactors = 1,
#'            numberOfrefGenes = 1,
#'            block = NULL,
#'            mainFactor.column = 1,
#'            p.adj = "none")
#'   
#' # Repeated measure analysis         
#' a <- ANOVA_DDCt(data_repeated_measure_1,
#'            numOfFactors = 1,
#'            numberOfrefGenes = 1,
#'            block = NULL,
#'            mainFactor.column = 1,
#'            p.adj = "none", model = wDCt ~ time + (1 | id))
#' 
#' a$perGene$Target$ANOVA_table
#' 
#' 
#' # Repeated measure analysis: split-plot in time
#' a <- ANOVA_DDCt(data_repeated_measure_2,
#'            numOfFactors = 2, numberOfrefGenes = 1,
#'            mainFactor.column = 2, block = NULL,
#'            model = wDCt ~ treatment * time + (1 | id))
#'            

ANOVA_DDCt <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    mainFactor.column,
    block,
    mainFactor.level.order = NULL,
    p.adj = "none",
    analyseAllTarget = TRUE,
    model = NULL
) {
  
  n <- ncol(x)
  nDesign <- if (is.null(block)) numOfFactors + 1 else numOfFactors + 2
  designCols <- seq_len(nDesign)
  nRefCols <- 2 * numberOfrefGenes
  refCols <- (n - nRefCols + 1):n
  targetCols <- setdiff(seq_len(n), c(designCols, refCols))
  
  if (length(targetCols) == 0 || length(targetCols) %% 2 != 0) {
    stop("Target genes must be supplied as E/Ct column pairs")
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
    
    gene_df <- gene_df[, c(mainFactor.column, (1:ncol(gene_df))[-mainFactor.column])]
    
    if (is.null(block)) {
      gene_df[seq_len(numOfFactors)] <- lapply(
        gene_df[seq_len(numOfFactors)],
        function(col) factor(col, levels = unique(col)))
    } else {
      gene_df[seq_len(numOfFactors + 1)] <- lapply(
        gene_df[seq_len(numOfFactors + 1)],
        function(col) factor(col, levels = unique(col)))
    }
    
    if (is.null(mainFactor.level.order)) {
      mainFactor.level.order <- unique(gene_df[,1])
      calibrator <- gene_df[,1][1]
    } else if (any(is.na(match(unique(gene_df[,1]), mainFactor.level.order)))) {
      stop("The `mainFactor.level.order` doesn't match main factor levels.")
    } else {
      gene_df <- gene_df[order(match(gene_df[,1], mainFactor.level.order)), ]
      calibrator <- gene_df[,1][1]
    }
    
    if (!exists("compute_wDCt")) {
      stop("compute_wDCt function is required but not found")
    }
    gene_df <- compute_wDCt(gene_df, numOfFactors, numberOfrefGenes, block)
    
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
      
      has_random_effects <- grepl("\\|", as.character(formula_obj)[3])
      
      if (has_random_effects) {
        if (!requireNamespace("lmerTest", quietly = TRUE)) {
          stop("lmerTest package is required for mixed models")
        }
        if (!requireNamespace("lme4", quietly = TRUE)) {
          stop("lme4 package is required for singularity checks")
        }
        
        is_mixed_model <- TRUE
        lm_fit <- suppressMessages(lmerTest::lmer(formula_obj, data = gene_df))
        
        # singularity check 
        is_singular <- lme4::isSingular(lm_fit, tol = 1e-4)
        
      } else {
        lm_fit <- lm(formula_obj, data = gene_df)
      }
      
      lm_formula <- formula(lm_fit)
      ANOVA_table <- stats::anova(lm_fit)
      
    } else {
      
      if (is.null(block)) {
        formula_ANOVA <- as.formula(
          paste("wDCt ~", paste(factors, collapse = " * "))
        )
        default_model_formula <- deparse(formula_ANOVA)
        lm_fit <- lm(formula_ANOVA, data = gene_df)
        
      } else {
        formula_ANOVA <- as.formula(
          paste("wDCt ~", block, "+", paste(factors, collapse = " * "))
        )
        default_model_formula <- deparse(formula_ANOVA)
        lm_fit <- lm(formula_ANOVA, data = gene_df)
      }
      
      lm_formula <- formula(lm_fit)
      ANOVA_table <- stats::anova(lm_fit)
    }
    
    if (!requireNamespace("emmeans", quietly = TRUE)) {
      stop("emmeans package is required for post-hoc tests")
    }
    
    pp1 <- suppressMessages(
      emmeans::emmeans(lm_fit, colnames(gene_df)[1],
                       data = gene_df, adjust = p.adj,
                       mode = "satterthwaite")
    )
    
    pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
    pp3 <- pp2[1:length(mainFactor.level.order)-1,]
    ci  <- as.data.frame(stats::confint(graphics::pairs(pp1)),
                         adjust = p.adj)[1:length(unique(gene_df[,1]))-1,]
    pp  <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
    
    bwDCt <- gene_df$wDCt
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("dplyr package is required")
    }
    
    se <- dplyr::summarise(
      dplyr::group_by(data.frame(factor = gene_df[,1], bwDCt = bwDCt),
                      gene_df[,1]),
      se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt))
    )
    
    sig <- .convert_to_character(pp$p.value)
    
    post_hoc_test <- data.frame(
      contrast = pp$contrast,
      ddCt = - pp$estimate,
      RE = 1/(2^-(pp$estimate)),
      log2FC = log2(1/(2^-(pp$estimate))),
      pvalue = pp$p.value,
      sig = sig,
      LCL = 1/(2^-pp$lower.CL),
      UCL = 1/(2^-pp$upper.CL),
      se = se$se[-1]
    )
    
    reference <- data.frame(
      contrast = mainFactor.level.order[1],
      ddCt = 0, RE = 1, log2FC = 0,
      pvalue = 1, sig = " ",
      LCL = 0, UCL = 0,
      se = se$se[1]
    )
    
    tableC <- rbind(reference, post_hoc_test)
    
    tableC$contrast <- sapply(
      strsplit(as.character(tableC$contrast), " - "),
      function(x) paste(rev(x), collapse = " vs ")
    )
    
    tableC <- data.frame(
      tableC,
      Lower.se.RE = 2^(log2(tableC$RE) - tableC$se),
      Upper.se.RE = 2^(log2(tableC$RE) + tableC$se),
      Lower.se.log2FC = 0,
      Upper.se.log2FC = 0
    )
    
    for (j in seq_len(nrow(tableC))) {
      if (tableC$RE[j] < 1) {
        tableC$Lower.se.log2FC[j] <- (tableC$Upper.se.RE[j]*log2(tableC$RE[j]))/tableC$RE[j]
        tableC$Upper.se.log2FC[j] <- (tableC$Lower.se.RE[j]*log2(tableC$RE[j]))/tableC$RE[j]
      } else {
        tableC$Lower.se.log2FC[j] <- (tableC$Lower.se.RE[j]*log2(tableC$RE[j]))/tableC$RE[j]
        tableC$Upper.se.log2FC[j] <- (tableC$Upper.se.RE[j]*log2(tableC$RE[j]))/tableC$RE[j]
      }
    }
    
    tableC[] <- lapply(tableC, function(x) if(is.numeric(x)) round(x, 5) else x)
    tableC$gene <- gene_name
     
    res <- list(
      Final_data = gene_df,
      lm = lm_fit,
      ANOVA_table = ANOVA_table,
      Fold_Change = tableC,
      lm_formula = lm_formula,
      user_defined_model = user_defined_model,
      default_model_formula = default_model_formula,
      is_mixed_model = is_mixed_model,
      singular = is_singular
    )
    
    perGene[[gene_name]] <- res
    relativeExpression_list[[i]] <- tableC
  }
  
  relativeExpression <- do.call(rbind, relativeExpression_list)
  rownames(relativeExpression) <- NULL
  
  relativeExpression <- relativeExpression[, c(
    "gene","contrast","ddCt","RE","log2FC",
    "LCL","UCL","se",
    "Lower.se.RE","Upper.se.RE",
    "Lower.se.log2FC","Upper.se.log2FC",
    "pvalue","sig"
  )]
  
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
      cat("Note: Using default model for analysis:", default_formula, "\n")
    }
  }
  

  invisible(list(
    perGene = perGene,
    relativeExpression = relativeExpression
  ))
}
