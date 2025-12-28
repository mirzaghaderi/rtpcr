# convert_to_character function
.convert_to_character <- function(numbers) {
  characters <- character(length(numbers))  # Initialize a character vector to store the results
  
  for (i in seq_along(numbers)) {
    if (numbers[i] < 0.001) {
      characters[i] <- "***"
    } else if (numbers[i] < 0.01) {
      characters[i] <- "**"
    } else if (numbers[i] < 0.05) {
      characters[i] <- "*"
    } else if (numbers[i] < 0.1) {
      characters[i] <- "."
    } else {
      characters[i] <- " "
    }
  }
  
  return(characters)
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
  out
}



.rearrange_repeatedMeasureData <- function(
    df,
    column_name,
    level
) {
  
  if (!column_name %in% names(df)) {
    stop("column_name not found in data frame")
  }
  
  col_idx <- match(column_name, names(df))
  if (col_idx == 1) {
    stop("No previous columns to define groups")
  }
  
  # grouping columns = all columns before column_name
  group_cols <- seq_len(col_idx - 1)
  
  # groups in order of appearance
  grp <- interaction(df[group_cols], drop = TRUE)
  grp_levels <- unique(grp)
  
  out <- df[0, , drop = FALSE]
  
  # safe comparison (numeric / factor / character)
  level_chr <- as.character(level)
  
  for (g in grp_levels) {
    block <- df[grp == g, , drop = FALSE]
    
    is_level <- as.character(block[[column_name]]) == level_chr
    
    block <- rbind(block[is_level, , drop = FALSE],
                   block[!is_level, , drop = FALSE])
    
    out <- rbind(out, block)
  }
  
  rownames(out) <- NULL
  out
}





.compute_wDCt <- function(x,
                          numberOfrefGenes,
                          block = NULL) {
  
  stopifnot(numberOfrefGenes >= 1)
  
  nc   <- ncol(x)
  nRef <- numberOfrefGenes
  
  # core columns: rep + target + refs
  nTargetCols <- 2
  nRefCols    <- 2 * nRef
  nCoreCols   <- 1 + nTargetCols + nRefCols
  
  if (!is.null(block)) {
    nCoreCols <- nCoreCols + 1
  }
  
  # rename from the end
  idx <- nc
  
  # reference genes
  for (i in nRef:1) {
    colnames(x)[idx - 1] <- paste0("Eref",  if (i == 1) "" else i)
    colnames(x)[idx]     <- paste0("Ctref", if (i == 1) "" else i)
    idx <- idx - 2
  }
  
  # target
  colnames(x)[idx - 1] <- "Etarget"
  colnames(x)[idx]     <- "Cttarget"
  idx <- idx - 2
  
  # replicate
  colnames(x)[idx] <- "rep"
  idx <- idx - 1
  
  # block (optional)
  if (!is.null(block)) {
    colnames(x)[idx] <- "block"
  }
  
  #FACTOR CONVERSION
  if (!is.null(block)) {
    factor_cols <- seq_len(which(colnames(x) == "block") - 1)
  } else {
    factor_cols <- seq_len(which(colnames(x) == "rep") - 1)
  }
  
  x[factor_cols] <- lapply(
    x[factor_cols],
    function(z) {
      if (is.numeric(z)) {
        z <- as.factor(z)
      } else {
        z <- factor(z)
      }
      z
    }
  )
  
  #wDCt calculation
  target_term <- log2(x$Etarget) * x$Cttarget
  
  ref_terms <- vapply(
    seq_len(nRef),
    function(i) {
      E  <- x[[paste0("Eref",  if (i == 1) "" else i)]]
      Ct <- x[[paste0("Ctref", if (i == 1) "" else i)]]
      log2(E) * Ct
    },
    numeric(nrow(x))
  )
  
  ref_mean <- rowMeans(ref_terms, na.rm = FALSE)
  
  y <- data.frame(
    x,
    wDCt = target_term - ref_mean
  )
  
  return(y)
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








.ANOVA_DDCt_uniTarget <- function(
    x, 
    numberOfrefGenes, 
    mainFactor.column, 
    analysisType = "anova",
    mainFactor.level.order = NULL, 
    block = NULL, 
    p.adj = "none",  
    plot = TRUE,
    plotType = "RE", 
    verbose = FALSE) {
  
  
  x <- x[, c(mainFactor.column, (1:ncol(x))[-mainFactor.column])]
  
  
  # basic checks
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame")
  }
  if (!is.numeric(numberOfrefGenes) || length(numberOfrefGenes) != 1) {
    stop("`numberOfrefGenes` must be a single numeric value")
  }
  if (!is.numeric(mainFactor.column) || length(mainFactor.column) != 1) {
    stop("`mainFactor.column` must be a single numeric value")
  }
  
  
  if (missing(numberOfrefGenes)) {
    stop("argument 'numberOfrefGenes' is missing, with no default")
  }
  if (missing(mainFactor.column)) {
    stop("argument 'mainFactor.column' is missing, with no default")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
  
  
  
  
  
  if (is.null(mainFactor.level.order)) {
    mainFactor.level.order <- unique(x[,1])
    calibrartor <- x[,1][1]
    on.exit(cat(structure(paste("*** The", calibrartor, "level was used as calibrator.\n"))))
  } else if (any(is.na(match(unique(x[,1]), mainFactor.level.order))) == TRUE){
    stop("The `mainFactor.level.order` doesn't match main factor levels.")
  } else {
    x <- x[order(match(x[,1], mainFactor.level.order)), ]
    #x <- x[order(match(as.character(x[,1]), as.character(mainFactor.level.order))), ]
    calibrartor <- x[,1][1]
    on.exit(cat(structure(paste("*** The", calibrartor, "level was used as calibrator.\n"))))
  }
  
  x <- .compute_wDCt(x, numberOfrefGenes, block)
  x[,1] <- factor(
    x[,1],
    levels = mainFactor.level.order
  )
  
  # get names of factor columns
  factors <- names(x)[vapply(x, is.factor, logical(1))]


  
  
  # Check if there is block
  if (is.null(block)) {
    
    # ANOVA based on factorial design
    formula_ANOVA <- as.formula(
      paste("wDCt ~", paste(factors, collapse = " * "))
    )
    lmf <- lm(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmf)
    
    # ANCOVA (other factors as covariates)
    formula_ANCOVA <- as.formula(
      paste("wDCt ~", paste(rev(factors), collapse = " + "))
    )
    lmc <- lm(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmc)
    
  } else {
    
    # ANOVA with blocking factor (block treated as fixed)
    formula_ANOVA <- as.formula(
      paste("wDCt ~ block +", paste(factors, collapse = " * "))
    )
    lmfb <- lm(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmfb)
    
    # ANCOVA with blocking factor
    formula_ANCOVA <- as.formula(
      paste("wDCt ~ block +", paste(rev(factors), collapse = " + "))
    )
    lmcb <- lm(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmcb)
    
  }
  
  
  
  
  
  
  # Type of analysis: ancova or anova
  if (is.null(block)) {
    if(analysisType == "ancova") {
      lm <- lmc
    } 
    else{
      lm <- lmf
    }
  } else {
    if(analysisType == "ancova") {
      lm <- lmcb
    } 
    else{
      lm <- lmfb
    } 
  }
  
  
  pp1 <- emmeans(lm, colnames(x)[1], data = x, adjust = p.adj, mode = "satterthwaite")
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  pp3 <- pp2[1:length(mainFactor.level.order)-1,]
  ci <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(unique(x[,1]))-1,]
  pp <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
  
  
  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(factor = x[,1], bwDCt = bwDCt), x[,1]),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))  
  
  
  sig <- .convert_to_character(pp$p.value)
  contrast <- pp$contrast
  post_hoc_test <- data.frame(contrast, 
                              RE = 1/(2^-(pp$estimate)),
                              log2FC = log2(1/(2^-(pp$estimate))),
                              pvalue = pp$p.value,
                              sig = sig,
                              LCL = 1/(2^-pp$lower.CL),
                              UCL = 1/(2^-pp$upper.CL),
                              se = se$se[-1])
  
  reference <- data.frame(contrast = mainFactor.level.order[1],
                          RE = 1,
                          log2FC = 0,
                          pvalue = 1, 
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])
  
  tableC <- rbind(reference, post_hoc_test)
  
  FINALDATA <- x
  
  tableC$contrast <- as.character(tableC$contrast)
  tableC$contrast <- sapply(strsplit(tableC$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))
  
  
  
  tableC$contrast <- factor(tableC$contrast, levels = unique(tableC$contrast))
  contrast <- tableC$contrast
  LCL <- tableC$LCL
  UCL <- tableC$UCL
  REp <- as.numeric(tableC$RE)
  FCp <- as.numeric(tableC$log2FC)
  significance <- tableC$sig
  se <- tableC$se
  
  tableC <- data.frame(tableC, 
                       Lower.se.RE = 2^(log2(tableC$RE) - tableC$se), 
                       Upper.se.RE = 2^(log2(tableC$RE) + tableC$se))  
  ##################################################
  a <- data.frame(tableC, d = 0)
  
  for (i in 1:length(tableC$RE)) {
    if (tableC$RE[i] < 1) {
      a$Lower.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] - 0.2
    } else {
      a$Lower.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] + 0.2
    }
  }
  pfc1 <- ggplot(a, aes(contrast,RE)) + 
    geom_col() +
    geom_errorbar(aes(ymin = tableC$Lower.se.RE, ymax=tableC$Upper.se.RE), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = tableC$Upper.se.RE + 0.2)) +
    ylab("Relative Expression (DDCt)")
  pfc2 <- ggplot(a, aes(contrast,log2FC)) +
    geom_col() +
    geom_errorbar(aes(ymin = Upper.se, ymax=Lower.se), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = d)) +
    ylab("log2FC")
  
  tableC <- data.frame(tableC, Lower.se.log2FC = a$Lower.se, Upper.se.log2FC = a$Upper.se)
  ##################################################  
  
  if (is.null(block)) {
    lm_ANOVA <- lmf
    lm_ANCOVA <- lmc
  } else {
    lm_ANOVA <- lmfb
    lm_ANCOVA <- lmcb
  }
  
  #round tableC to 4 decimal places
  tableC <- tableC %>%
    mutate_if(is.numeric, ~ round(., 4))
  
  outlist2 <- structure(list(Final_data = x,
                             lm_ANOVA = lm_ANOVA,
                             lm_ANCOVA = lm_ANCOVA,
                             ANOVA_table = ANOVA,
                             ANCOVA_table = ANCOVA,
                             Fold_Change  = tableC,
                             RE_Plot = pfc1,
                             log2FC_Plot = pfc2), class = "XX")
  
  print.XX <- function(outlist2){
    cat("ANOVA table", "\n")
    print(outlist2$ANOVA_table)
    cat("\n", sep = '', "ANCOVA table", "\n")
    print(outlist2$ANCOVA_table)
    cat("\n", sep = '', "Expression table", "\n")
    print(outlist2$Fold_Change)
    
    
    if (plot == TRUE){
      if(plotType == "RE"){
        cat("RE_Plot\n")
        print(outlist2$RE_Plot)
      }else{
        cat("log2FC_Plot\n")
        print(outlist2$log2FC_Plot)
      }
    }
    invisible(outlist2)
  }
  print.XX(outlist2)
}







.ANOVA_DCt_uniTarget <- function(x,
                      numberOfrefGenes,
                      block,
                      alpha,
                      adjust,
                      verbose = FALSE) {
  
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
                        sig = trimws(letters_grp), 
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
  cols <- colnames(Results_combined)
  letters_pos <- match("sig", cols)
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








.REPEATED_DDCt_uniTarget <- function(x, 
                          numberOfrefGenes,
                          repeatedFactor, 
                          calibratorLevel,
                          block,
                          p.adj,
                          plot,
                          plotType){
  
  # basic checks
  # if (!is.data.frame(x)) stop("`x` must be a data.frame")
  # if (missing(repeatedFactor)) stop("argument 'repeatedFactor' is missing")
  # if (missing(calibratorLevel)) stop("argument 'calibratorLevel' is missing")
  # if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1)
  #   stop("`numberOfrefGenes` must be >= 1")
  # if (missing(block)) stop("argument 'block' is missing")
  
  # rearrange_repeatedMeasureData
  x <- .rearrange_repeatedMeasureData(x, column_name = "Time_", level = calibratorLevel) 
  
  
  ## validate number of target genes
  ## validate that only ONE target gene exists
  expr_cols_expected <- if (is.null(block)) {
    3 + 2 * numberOfrefGenes   # time + target + refs
  } else {
    4 + 2 * numberOfrefGenes   # block + time + target + refs
  }
  
  non_expr_cols <- ncol(x) - expr_cols_expected
  
  if (non_expr_cols < 1) {
    stop(
      "Input data structure error:\n",
      "At least one non-expression column (id) must exist before expression columns.",
      call. = FALSE
    )
  }
  
  ## if expression columns are MORE than expected → extra target genes
  actual_expr_cols <- ncol(x) - non_expr_cols
  
  if (actual_expr_cols != expr_cols_expected) {
    stop(
      sprintf(
        paste0(
          "Exactly ONE target gene is allowed.\n\n",
          "Expected expression columns:\n",
          "  %d  (= time + 1 target + %d reference gene(s)%s)\n\n",
          "But detected:\n",
          "  %d expression-related columns\n\n",
          "This usually means:\n",
          "  more than one target gene is present, or\n",
          "  numberOfrefGenes is incorrect, or\n",
          "  expression columns are not at the end of the data frame."
        ),
        expr_cols_expected,
        numberOfrefGenes,
        if (is.null(block)) "" else " + block",
        actual_expr_cols
      ),
      call. = FALSE
    )
  }
  

  id <- colnames(x)[1]
  
  # column parsing
  if (is.null(block)) {
    
    n_expr <- 3 + 2 * numberOfrefGenes
    factors <- if ((ncol(x) - n_expr) <= 1) NULL else colnames(x)[2:(ncol(x) - n_expr)]
    
    colnames(x)[(ncol(x) - n_expr + 1)] <- "Time_"
    colnames(x)[(ncol(x) - n_expr + 2)] <- "Etarget"
    colnames(x)[(ncol(x) - n_expr + 3)] <- "Cttarget"
    
    ref_start <- ncol(x) - (2 * numberOfrefGenes) + 1
    ref_cols <- ref_start:ncol(x)
    
  } else {
    
    n_expr <- 4 + 2 * numberOfrefGenes
    factors <- if ((ncol(x) - n_expr) <= 1) NULL else colnames(x)[2:(ncol(x) - n_expr)]
    
    colnames(x)[(ncol(x) - n_expr + 1)] <- "Time_"
    colnames(x)[(ncol(x) - n_expr + 2)] <- "block"
    colnames(x)[(ncol(x) - n_expr + 3)] <- "Etarget"
    colnames(x)[(ncol(x) - n_expr + 4)] <- "Cttarget"
    
    ref_start <- ncol(x) - (2 * numberOfrefGenes) + 1
    ref_cols <- ref_start:ncol(x)
  }
  
  # compute wDCt (GENERALIZED)
  target_part <- log2(x$Etarget) * x$Cttarget
  
  ref_matrix <- matrix(
    mapply(
      function(E, Ct) log2(E) * Ct,
      x[, ref_cols[seq(1, length(ref_cols), 2)]],
      x[, ref_cols[seq(2, length(ref_cols), 2)]]
    ),
    ncol = numberOfrefGenes
  )
  
  ref_part <- rowMeans(ref_matrix)
  x <- data.frame(x, wDCt = target_part - ref_part)
  
  # convert factors
  for (i in 2:which(names(x) == "Time_")) {
    x[[i]] <- factor(x[[i]], levels = unique(x[[i]]))
  }

  
  
  # model formula 
  if (is.null(block)) {
    if (is.null(factors)) {
      formula <- wDCt ~ Time_ + (1 | id)
    } else {
      formula <- as.formula(
        paste("wDCt ~ Time_ *", paste(factors, collapse = " * "), "+ (1 | id)")
      )
    }
  } else {
    if (is.null(factors)) {
      formula <- wDCt ~ Time_ + (1 | id) + (1 | block)
    } else {
      formula <- as.formula(
        paste("wDCt ~ Time_ *", paste(factors, collapse = " * "),
              "+ (1 | id) + (1 | block)")
      )
    }
  }
  
  lm <- lmerTest::lmer(formula, data = x)
  ANOVA <- stats::anova(lm)
  
  #post hoc
  v <- match(colnames(x), "Time_") 
  n <- which(!is.na(v))
  repeatedFactor <- colnames(x)[n]
  lvls <- unique(x[,n])
  calibrartor <- lvls[1]
  
  on.exit(cat(paste("The level", calibrartor, " of the selected factor was used as calibrator.\n")))
  pp1 <- emmeans(lm, repeatedFactor, data = x, adjust = p.adj, mode = "satterthwaite")
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  if (length(lvls) >= 3){
    pp3 <- pp2[1:length(lvls) - 1,] 
  } else {
    pp3 <- pp2
  }
  ci <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(lvls)-1,]
  pp <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
  
  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(repeatedFactor = x[n], bwDCt = bwDCt), x[n]),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))  
  
  
  sig <- .convert_to_character(pp$p.value)
  contrast <- pp$contrast
  post_hoc_test <- data.frame(contrast, 
                              RE = 1/(2^-(pp$estimate)),
                              log2FC = log2(1/(2^-(pp$estimate))),
                              pvalue = pp$p.value,
                              sig = sig,
                              LCL = 1/(2^-pp$lower.CL),
                              UCL = 1/(2^-pp$upper.CL),
                              se = se$se[-1])
  
  
  words <- strsplit(as.character(contrast[1]), " ")[[1]]
  referencelevel <- words[1]
  
  
  reference <- data.frame(contrast = as.character(referencelevel),
                          RE = 1,
                          log2FC = 0,
                          pvalue = 1, 
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])
  
  tableC  <- rbind(reference, post_hoc_test)

  FINALDATA <- x
  
  tableC$contrast <- sapply(strsplit(tableC$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))
  
  
  
  tableC$contrast <- factor(tableC$contrast, levels = unique(tableC$contrast))
  contrast <- tableC$contrast
  LCL <- tableC$LCL
  UCL <- tableC$UCL
  FCp <- as.numeric(tableC$FC)
  significance <- tableC$sig
  se <- tableC$se
  
  
  tableC <- data.frame(tableC, 
                       Lower.se.RE = 2^(log2(tableC$RE) - tableC$se), 
                       Upper.se.RE = 2^(log2(tableC$RE) + tableC$se))  

  
  
  a <- data.frame(tableC, d = 0)
  
  for (i in 1:length(tableC$RE)) {
    if (tableC$RE[i] < 1) {
      a$Lower.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] - 0.2
    } else {
      a$Lower.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] + 0.2
    }
  }
  pfc1 <- ggplot(a, aes(contrast,RE)) + 
    geom_col() +
    geom_errorbar(aes(ymin = tableC$Lower.se.RE, ymax=tableC$Upper.se.RE), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = tableC$Upper.se.RE + 0.2)) +
    ylab("Relative Expression (DDCt)")
  pfc2 <- ggplot(a, aes(contrast,log2FC)) +
    geom_col() +
    geom_errorbar(aes(ymin = Upper.se, ymax=Lower.se), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = d)) +
    ylab("log2FC")
  
  tableC <- data.frame(tableC, Lower.se.log2FC = a$Lower.se, Upper.se.log2FC = a$Upper.se)

  
  
  tableC <- tableC %>%
    mutate_if(is.numeric, ~ round(., 4))
  
  outlist2 <- structure(list(Final_data = x,
                             lm = lm,
                             ANOVA_table = ANOVA,
                             Relative_Expression_table  = tableC,
                             RE_Plot = pfc1,
                             log2FC_Plot = pfc2), class = "XX")
  
  print.XX <- function(outlist2){
    print(outlist2$ANOVA_table)
    cat("\n", sep = '',"Expression table", "\n")
    print(outlist2$Relative_Expression_table)
    
    if (plot == TRUE){
      if(plotType == "RE"){
        cat("\n", sep = '', "Expression plot", "\n")
        print(outlist2$RE_Plot)
      }else{
        cat("\n", sep = '', "Expression plot", "\n")
        print(outlist2$log2FC_Plot)
      }
    }
    
    invisible(outlist2)
  }
  print.XX(outlist2)
}
