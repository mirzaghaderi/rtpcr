#' @title Delta Delta Ct pairwise comparisons using fitted ANOVA_DDCt object
#'
#' @description
#' Performs relative expression (fold change) pairwise comparisons for all
#' genes using a fitted object returned by \code{ANOVA_DDCt()}.
#'
#' Standard errors are computed from the raw data grouped according to the
#' \code{specs} argument (same structure used by \pkg{emmeans}).
#'
#' @export
#'
#' @param object Object returned by \code{ANOVA_DDCt()}
#' @param specs emmeans specification (e.g. "Factor", "A | B")
#' @param p.adj p-value adjustment method
#' @param modelBased_se Logical; if TRUE use model residuals, else raw wDCt will be used for computing se.
#' @param se.type Character string specifying how standard error is calculated. 
#' One of \code{"paired.group"}, \code{"two.group"}, or \code{"single.group"}. 
#' \code{"paired.group"} computes SE from paired differences (used when a random 
#' \code{id} effect is present), \code{"two.group"} uses the unpaired two-group 
#' t-test standard error against the reference level, and \code{"single.group"} 
#' computes SE within each level using a one-group t-test.
#' @param numOfFactors Number of factors in the input data set excluding block.
#' @param block Character. Block column name or \code{NULL}.
#' @return Data frame of fold changes and SE for all genes
#'
#' @examples
#' a <- ANOVA_DDCt(data_repeated_measure_2,
#' numOfFactors = 2, numberOfrefGenes = 1, model = wDCt ~ treatment * time + (1 | id),
#' mainFactor.column = 2, block = NULL, se.type = "two.group")
#' 
#' Means_DDCt2(a,
#'             specs = "time | treatment",
#'             p.adj = "none",
#'             block = NULL,
#'             numOfFactors = 2, se.type = "single.group",
#'             modelBased_se = TRUE)
#' 
#' 
#' res <- ANOVA_DDCt(
#'   data_2factorBlock3ref,
#'   numOfFactors = 2,
#'   numberOfrefGenes = 1,
#'   mainFactor.column = 1,
#'   block = "block")
#' 
#' Means_DDCt2(res, specs = "Concentration | Type",
#'             block = "block",
#'             numOfFactors = 2, se.type = "two.group",
#'             modelBased_se = TRUE)
#' 
Means_DDCt2 <- function(x,
                        specs,
                        numOfFactors,
                        block,
                        se.type = c("paired.group", "two.group", "single.group"),
                        modelBased_se = TRUE,
                        p.adj = "none") {
  

  if (!is.list(object) || is.null(object$perGene)) {
    stop("object must be created by ANOVA_DDCt()")
  }
  
  if (!grepl("\\|", specs)) {
    stop("No subgroup was determined in specs. The 'specs' argument must contain a '|' sign to specify grouping variables (e.g., 'A | B' or 'A | B * C').")
  }
  
  se.type <- match.arg(se.type)
  specs_formula <- stats::as.formula(paste("~", specs))
  

  # Extract grouping variables (after | if present)
  split_specs <- strsplit(specs, "\\|")[[1]]
  grouping_vars <- character(0)
  if (length(split_specs) > 1) {
    rhs_group <- trimws(split_specs[2])
    grouping_vars <- all.vars(stats::as.formula(paste("~", rhs_group)))
  }
  
  detect_rep_id_random <- function(model, numOfFactors, block, x) {
    if (is.null(model)) return(FALSE)
    ftxt <- paste(deparse(formula(model)), collapse = " ")
    if (is.null(numOfFactors)) stop("numOfFactors must be provided to detect random ID")
    rep_id_col <- if (is.null(block)) numOfFactors + 1 else numOfFactors + 2
    rep_id_name <- colnames(x)[rep_id_col]
    grepl(paste0("\\|\\s*", rep_id_name, "\\b"), ftxt)
  }
  
  results_all <- lapply(names(object$perGene), function(gene) {
    
    gene_obj <- object$perGene[[gene]]
    model    <- gene_obj$lm
    gene_df  <- gene_obj$Final_data
    
    if (!any(c("lmerModLmerTest", "lmerMod", "lm") %in% class(model))) {
      stop("Model class not supported for gene: ", gene)
    }
    

    # EMMEANS
    pair_formula <- stats::as.formula(paste("pairwise ~", specs))
    suppressMessages(emm_pair <- emmeans::emmeans(model, pair_formula, adjust = p.adj))
    suppressMessages(emm_ci <- stats::confint(emm_pair))
    contr_df <- as.data.frame(emm_ci$contrasts)
    contr_df$p.value <- as.data.frame(emm_pair$contrasts)$p.value
    

    # Back-transform ddCt -> fold change
    contr_df$RE  <- 2^(contr_df$estimate)
    contr_df$LCL <- 2^(contr_df$lower.CL)
    contr_df$UCL <- 2^(contr_df$upper.CL)
    
    contr_df$contrast <- as.character(contr_df$contrast)
    contr_df$contrast <- sapply(
      strsplit(contr_df$contrast, " - "),
      function(x) paste(rev(x), collapse = " vs ")
    )
    

    # Compute bwDCt
    if (!modelBased_se) {
      gene_df$bwDCt <- gene_df$wDCt
    } else {
      gene_df$bwDCt <- residuals(model, type = "response")
    }
    
    vars <- all.vars(specs_formula)
    idRand <- detect_rep_id_random(model, numOfFactors, block, gene_df)
    if (idRand && se.type != "paired.group") {
      warning("Random ID detected: using paired.group SE")
    }

    
    # Compute SE per group
    # Single-group SE
    if (se.type == "single.group") {
      se_df <- gene_df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(vars))) %>%
        dplyr::summarise(
          n = dplyr::n(),
          se = ifelse(n >= 2, stats::sd(bwDCt, na.rm = TRUE)/sqrt(n), NA_real_),
          .groups = "drop"
        )
      
    } else if (se.type == "two.group") {
      # Two-group SE: each level vs first level
      se_df <- gene_df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(vars))) %>%
        dplyr::summarise(
          bwDCt_list = list(bwDCt),
          .groups = "drop"
        )
      
      se_vec <- numeric(nrow(se_df))
      ref_vals <- se_df$bwDCt_list[[1]]
      se_vec[1] <- NA_real_  # SE for reference vs itself is NA
      
      for (k in 2:nrow(se_df)) {
        grp_vals <- se_df$bwDCt_list[[k]]
        # Two-sample unpaired t-test SE
        se_vec[k] <- tryCatch(
          stats::t.test(grp_vals, ref_vals, paired = FALSE)$stderr,
          error = function(e) NA_real_
        )
      }
      
      se_df$se <- se_vec
      
    } else if (se.type == "paired.group" || idRand) {
      # Paired SE requires id column
      if (!"id" %in% colnames(gene_df)) stop("Paired SE requested but no id column found.")
      
      se_df <- gene_df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(vars))) %>%
        dplyr::summarise(
          bwDCt_list = list(bwDCt),
          id_list   = list(id),
          .groups = "drop"
        )
      
      se_vec <- numeric(nrow(se_df))
      se_vec[1] <- NA_real_
      ref_vals <- se_df$bwDCt_list[[1]]
      ref_ids  <- se_df$id_list[[1]]
      
      for (k in 2:nrow(se_df)) {
        grp_vals <- se_df$bwDCt_list[[k]]
        grp_ids  <- se_df$id_list[[k]]
        merged <- merge(
          data.frame(id=ref_ids, bwDCt=ref_vals),
          data.frame(id=grp_ids, bwDCt=grp_vals),
          by="id"
        )
        if (nrow(merged) < 2) {
          se_vec[k] <- NA_real_
        } else {
          d <- merged$bwDCt.x - merged$bwDCt.y
          se_vec[k] <- sqrt(stats::var(d, na.rm=TRUE)/nrow(merged))
        }
      }
      
      se_df$se <- se_vec
    }


    # Map SE to contrasts using the computed SE per group
    emm_grid <- as.data.frame(emm_pair$emmeans)
    # Create keys to match groups
    grid_key <- do.call(paste, emm_grid[vars])
    se_key   <- do.call(paste, se_df[vars])
    se_ordered <- se_df$se[match(grid_key, se_key)]
    
    coef_mat <- as.matrix(emm_pair$contrasts@linfct)
    se_contrast <- apply(coef_mat, 1, function(w) {
      idx <- which(w != 0)
      sqrt(sum((w[idx]^2) * (se_ordered[idx]^2), na.rm=TRUE))
    })
    
    contr_df$se_group <- se_contrast
    contr_df$log2FC <- log2(contr_df$RE)
    contr_df$Lower.se.RE <- 2^(log2(contr_df$RE) - contr_df$se_group)
    contr_df$Upper.se.RE <- 2^(log2(contr_df$RE) + contr_df$se_group)
    contr_df$L.se.log2FC <- log2(contr_df$RE) - contr_df$se_group
    contr_df$U.se.log2FC <- log2(contr_df$RE) + contr_df$se_group
    
    contr_df$sig <- .convert_to_character(contr_df$p.value)
    contr_df$gene <- gene
    rownames(contr_df) <- NULL
    

    # Final output
    final_cols <- c(
      "gene",
      grouping_vars[grouping_vars %in% colnames(contr_df)],
      "contrast",
      "RE", "log2FC", "LCL", "UCL",
      "se_group", "Lower.se.RE", "Upper.se.RE", "L.se.log2FC", "U.se.log2FC",
      "p.value", "sig"
    )
    
    contr_df[, final_cols, drop = FALSE]
  })
  
  dplyr::bind_rows(results_all)
}
