.convert_to_character <- function(numbers) {
  characters <- rep(" ", length(numbers))
  ok <- !is.na(numbers)
  characters[ok & numbers < 0.001] <- "***"
  characters[ok & numbers >= 0.001 & numbers < 0.01] <- "**"
  characters[ok & numbers >= 0.01  & numbers < 0.05] <- "*"
  characters[ok & numbers >= 0.05  & numbers < 0.1 ] <- "."
  return(characters)
}




.cleanup <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    block = NULL,
    set_missing_target_Ct_to_40 = FALSE  # Default FALSE = NA, TRUE = 40
) {
  
  x <- as.data.frame(x)
  
  # Set missing value based on parameter
  missing_value <- if (set_missing_target_Ct_to_40) 40 else NA
  
  # Number of non-gene columns
  n_meta <- numOfFactors + if (!is.null(block)) 2 else 1
  
  # Total gene columns
  gene_cols <- (n_meta + 1):ncol(x)
  
  if (length(gene_cols) %% 2 != 0) stop("Gene columns must be in E/Ct pairs")
  
  # Split into pairs
  gene_pairs <- split(gene_cols, ceiling(seq_along(gene_cols) / 2))
  
  # Reference vs target genes
  ref_pairs <- tail(gene_pairs, numberOfrefGenes)
  target_pairs <- utils::head(gene_pairs, length(gene_pairs) - numberOfrefGenes)
  
  # Helper: check invalid (NA or any string/non-numeric)
  is_invalid <- function(v) {
    is.na(v) | is.na(suppressWarnings(as.numeric(v)))
  }
  
  # Reference genes: NA for invalid values 
  for (pair in ref_pairs) {
    E_col <- pair[1]
    Ct_col <- pair[2]
    
    bad_E <- is_invalid(x[[E_col]])
    x[[E_col]][bad_E] <- NA
    x[[E_col]] <- as.numeric(x[[E_col]])
    
    bad_Ct <- is_invalid(x[[Ct_col]])
    x[[Ct_col]][bad_Ct] <- NA  # Always NA for reference genes
    x[[Ct_col]] <- as.numeric(x[[Ct_col]])
  }
  
  # Target genes - use parameter choice ----
  for (pair in target_pairs) {
    E_col  <- pair[1]
    Ct_col <- pair[2]
    
    bad_E <- is_invalid(x[[E_col]])
    x[[E_col]][bad_E] <- NA
    x[[E_col]] <- as.numeric(x[[E_col]])
    
    bad_Ct <- is_invalid(x[[Ct_col]])
    x[[Ct_col]][bad_Ct] <- missing_value  # NA or 40 based on parameter
    x[[Ct_col]] <- as.numeric(x[[Ct_col]])
  }
  
  return(x)
}






.wide_to_long <- function(df) {
  
  if (ncol(df) < 6) {
    stop("Data frame must contain at least 6 columns.")
  }
  
  # metadata (first two columns)
  meta <- df[, 1:2, drop = FALSE]
  
  # remaining columns (paired by position)
  data_cols <- df[, -(1:2), drop = FALSE]
  
  if (ncol(data_cols) %% 2 != 0) {
    stop("After the first two columns, remaining columns must be in pairs.")
  }
  
  n_pairs <- ncol(data_cols) / 2
  
  out <- do.call(
    rbind,
    lapply(seq_len(n_pairs), function(i) {
      
      e_col  <- data_cols[, 2*i - 1]
      ct_col <- data_cols[, 2*i]
      
      gene_name <- colnames(data_cols)[2*i - 1]
      
      data.frame(
        Condition = meta[[1]],
        Gene = gene_name,
        E = e_col,
        Ct = ct_col,
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  out[[1]] <- factor(out[[1]])
  out
}



.long_to_wide <- function(df) {
  
  if (ncol(df) < 4) {
    stop("Data frame must contain at least 4 columns: Condition, Gene, E, Ct")
  }
  
  # standardize first 4 columns internally
  tmp <- data.frame(
    Condition = df[[1]],
    Gene      = df[[2]],
    E         = df[[3]],
    Ct        = df[[4]],
    stringsAsFactors = FALSE
  )
  
  # replicate number within Condition C Gene
  tmp$Rep <- ave(
    seq_len(nrow(tmp)),
    tmp$Condition,
    tmp$Gene,
    FUN = seq_along
  )
  
  # reshape to wide
  wide <- reshape(
    tmp,
    idvar   = c("Condition", "Rep"),
    timevar = "Gene",
    direction = "wide"
  )
  wide[[1]] <- factor(wide[[1]])
  .cleanup(wide)
}



.geom_pub_cols <- function(col_width = 0.8,
                           err_width = 0.15,
                           fill_colors = NULL,
                           dodge_width = 0.8,
                           alpha = 1,
                           ...) {
  
  pos <- position_dodge(width = dodge_width)
  
  layers <- list(
    geom_col(width = col_width, position = pos, alpha = alpha, ...),
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = err_width, position = pos)
  )
  
  if (!is.null(fill_colors)) {
    layers <- c(layers, scale_fill_manual(values = fill_colors))
  }
  
  layers
}


.theme_pub <- function(base_size = 12,
                       base_family = "sans",
                       legend_position = "right",
                       ...) {
  
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      legend.position = legend_position,
      ...
    )
}






.ANOVA_DDCtx <- function(
    x,
    numOfFactors,
    numberOfrefGenes,
    mainFactor.column,
    block,
    mainFactor.level.order = NULL,
    p.adj = "none",
    analyseAllTarget = TRUE,
    model = NULL,
    set_missing_target_Ct_to_40 = FALSE,
    se.type = c("paired.group", "two.group", "single.group"),
    modelBased_se = TRUE
) {
  
  
  se.type <- match.arg(se.type)
  
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
      
      # has_random_effects <- grepl("\\|", as.character(formula_obj)[3])
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
    
    
    pp1 <- suppressMessages(
      emmeans::emmeans(lm_fit, colnames(gene_df)[1],
                       adjust = p.adj))    # mode = "satterthwaite"      data = gene_df, 
    
    pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
    pp3 <- pp2[1:length(mainFactor.level.order)-1,]
    ci  <- as.data.frame(stats::confint(graphics::pairs(pp1)),   # , na.action = stats::na.pass
                         adjust = p.adj)[1:length(unique(gene_df[,1]))-1,]
    pp  <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
    
    
    
    
    
    
    
    
    if (!modelBased_se) {
      bwDCt <- gene_df$wDCt
    } else {
      bwDCt <- residuals(lm_fit, type = "response")
    }
    
    
    
    
    idRand <- detect_rep_id_random(model = model,
                                   numOfFactors = numOfFactors,
                                   block = block, x = x)
    
    se.type <- match.arg(se.type, c("single.group", "two.group", "paired.group"))
    
    grp_levels <- levels(factor(gene_df[[1]]))
    n_groups  <- length(grp_levels)
    ref_level <- grp_levels[1]
    
    
    #if (idRand) 
      if (se.type == "paired.group") {
        warning("paired se is only appropriate for paired groups or repeated measure experiments!")  
      # id column is immediately before E/Ct gene columns
        df_se <- data.frame(
          factor = gene_df[[1]],
          wDCt   = bwDCt
        )
        
        ref_vals <- df_se$wDCt[df_se$factor == ref_level]
        
        se <- data.frame(
          factor = grp_levels,
          se = sapply(grp_levels, function(g) {
            if (g == ref_level) return(0)
            grp_vals <- df_se$wDCt[df_se$factor == g]
            
            tryCatch(
              stats::t.test(grp_vals, ref_vals, paired = TRUE)$stderr,
              error = function(e) NA_real_
            )
          })
        )
      
    } else if (se.type == "single.group") {
        df_se <- data.frame(
          factor = gene_df[[1]],
          wDCt   = bwDCt
        )
        
        se <- dplyr::summarise(
          dplyr::group_by(df_se, factor),
          se = tryCatch(
            stats::t.test(wDCt)$stderr,
            error = function(e) NA_real_
          ),
          .groups = "drop"
        )
        
      } else if (se.type == "two.group") {
        df_se <- data.frame(
          factor = gene_df[[1]],
          wDCt   = bwDCt
        )
        
        ref_vals <- df_se$wDCt[df_se$factor == ref_level]
        
        se <- data.frame(
          factor = grp_levels,
          se = sapply(grp_levels, function(g) {
            if (g == ref_level) return(0)
            grp_vals <- df_se$wDCt[df_se$factor == g]
            
            tryCatch(
              stats::t.test(grp_vals, ref_vals, paired = FALSE)$stderr,
              error = function(e) NA_real_
            )
          })
        )
      }
    
    
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
    #post_hoc_test$sig[post_hoc_test$RE < 0.001] <- "ND"
    post_hoc_test$RE[post_hoc_test$pvalue == "NaN"] <- 0
    post_hoc_test$log2FC[post_hoc_test$pvalue == "NaN"] <- 0
    post_hoc_test$sig[post_hoc_test$pvalue == "NaN"] <- "ND"
    
    
    reference <- data.frame(
      contrast = mainFactor.level.order[1],
      ddCt = 0, RE = 1, log2FC = 0,
      pvalue = 1, sig = " ",
      LCL = 0, UCL = 0,
      se = se$se[1])
    
    tableC <- rbind(reference, post_hoc_test)
    
    tableC$contrast <- sapply(
      strsplit(as.character(tableC$contrast), " - "),
      function(x) paste(rev(x), collapse = " vs ")
    )
    
    tableC <- data.frame(
      tableC,
      Lower.se.RE = 2^(log2(tableC$RE) - tableC$se),
      Upper.se.RE = 2^(log2(tableC$RE) + tableC$se),
      Lower.se.log2FC = log2(tableC$RE) - tableC$se,
      Upper.se.log2FC = log2(tableC$RE) + tableC$se
    )
    
    
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
  
  
  for (col in names(relativeExpression)) {
    if (is.numeric(relativeExpression[[col]])) {
      relativeExpression[[col]] <- round(relativeExpression[[col]], 5)
    }
  }
  
  
  
  
  # cat("\nRelative Expression\n")
  # print(relativeExpression)
  # cat("\n")
  
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
