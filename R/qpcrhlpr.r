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




.cleanup <- function(x, numOfFactors, block) {
  
  # define factor columns 
  if (is.null(block)) {
    x[seq_len(numOfFactors)] <- lapply(x[seq_len(numOfFactors)], factor)
  } else {
    x[seq_len(numOfFactors + 1)] <- lapply(x[seq_len(numOfFactors + 1)], factor)
  }

  # pairwise cleanup (from last to first) 
  n <- ncol(x)
  
  for (i in seq(from = n, to = 2, by = -2)) {
    
    colA <- x[[i - 1]]
    colB <- x[[i]]
    
    # skip if either column is a factor
    if (is.factor(colA) || is.factor(colB)) next
    
    badA <- suppressWarnings(is.na(as.numeric(colA))) | colA == 0
    badB <- suppressWarnings(is.na(as.numeric(colB))) | colB == 0
    
    colA[badB] <- NA
    colB[badA] <- NA
    
    x[[i - 1]] <- colA
    x[[i]]     <- colB
  }

  
  x[] <- lapply(x, function(col) {
    if (is.factor(col)) return(col)
    if (is.character(col)) {
      col[col %in% c("Undetermined", "undetermined")] <- NA
      suppressWarnings(col <- as.numeric(col))
    }
    if (is.numeric(col)) {
      col[col == 0] <- NA
    }
    col
  })
  x
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







.ANOVA_DCt_uniTarget <- function(x,
                                 numOfFactors,
                                 numberOfrefGenes,
                                 block,
                                 alpha,
                                 p.adj,
                                 verbose = FALSE) {
  
  # basic argument checks
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
  if (!is.numeric(numberOfrefGenes) || length(numberOfrefGenes) != 1) {
    stop("`numberOfrefGenes` must be a single numeric value")
  }
  
  
  x <- compute_wDCt(x, numOfFactors, numberOfrefGenes, block) 
  
  # Convert all factor columns to character
  x[] <- lapply(x, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  # convert the first numOfFactors columns to factor
  if (is.null(block)) {
    x[seq_len(numOfFactors)] <- lapply(
      x[seq_len(numOfFactors)],
      function(col) factor(col, levels = unique(col)))
  } else {
    x[seq_len(numOfFactors + 1)] <- lapply(
      x[seq_len(numOfFactors + 1)],
      function(col) factor(col, levels = unique(col)))
  }
  
  # get names of factor columns
  factors <- colnames(x)[1:numOfFactors]
  

  # build treatment factor T and fit lm
  if (is.null(block)) {
    x$T <- do.call(paste, c(x[1:length(factors)], sep = ":"))
    x$T <- as.factor(x$T)
    lm_fit <- stats::lm(wDCt ~ T, data = x)
    anovaCRD <- stats::anova(lm_fit)
  } else {
    x$T <- do.call(paste, c(x[1:length(factors)], sep = ":"))
    x$T <- as.factor(x$T)
    
    lm_fit <- as.formula(
      paste("wDCt ~", block, "+ T"))

    lm_fit <- stats::lm(lm_fit, data = x)
    anovaCRD <- stats::anova(lm_fit)
  }
  
  

  # LM factorial 
  if (is.null(block)) {
    
    # ANOVA based on factorial design      
    formula_ANOVA <- as.formula(
      paste("wDCt ~", paste(factors, collapse = " * "))
    )
    lm_factorial <- lm(formula_ANOVA, data = x)
    ANOVA_factorial <- stats::anova(lm_factorial)
    lm_formula <- formula(lm_factorial)
  } else {
    # ANOVA with blocking factor (block treated as fixed)
    formula_ANOVA <- as.formula(
      paste("wDCt ~", block, "+", paste(factors, collapse = " * "))
    )
    lm_factorial <- lm(formula_ANOVA, data = x)
    ANOVA_factorial <- stats::anova(lm_factorial)
    lm_formula <- formula(lm_factorial)
  }

  
  # emmeans / multiple comparisons
  emg <- suppressMessages(emmeans::emmeans(lm_fit, pairwise ~ T, mode = "satterthwaite"))
  # use cld() on the emmeans object (the first element)
  meanPairs <- multcomp::cld(emg[[1]], adjust = p.adj, alpha = alpha, reversed = FALSE, Letters = letters)
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
                        dCt = diffs,
                        RE = 2^(-diffs),
                        log2FC = log2(2^(-diffs)),
                        LCL = 2^(-lcl),
                        UCL = 2^(-ucl),
                        se = se_matched,
                        sig = trimws(letters_grp), 
                        stringsAsFactors = FALSE)
  
  # preserve rownames as a column for splitting
  Results$RowNames <- rownames(Results)
  
  
  # split RowNames back to factor columns (base R) ----
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
    mutate_if(is.numeric, ~ round(., 5))
  
  outlist2 <- list(Final_data = xx,
                   lm_T = lm_fit,
                   lm_factorial = lm_factorial,
                   ANOVA_T = anovaCRD,
                   ANOVA_factorial = ANOVA_factorial,
                   Results = Results_final)
  
  invisible(outlist2)
}

















.ANOVA <- function(
    x, 
    numOfFactors,
    numberOfrefGenes, 
    mainFactor.column, 
    analysisType,
    mainFactor.level.order = NULL, 
    block, 
    p.adj = "none",
    verbose = FALSE) {

  x <- x[, c(mainFactor.column, (1:ncol(x))[-mainFactor.column])]

  
  #convert the first numOfFactors columns to factor
  if (is.null(block)) {
    x[seq_len(numOfFactors)] <- lapply(
      x[seq_len(numOfFactors)],
      function(col) factor(col, levels = unique(col)))
  } else {
    x[seq_len(numOfFactors + 1)] <- lapply(
      x[seq_len(numOfFactors + 1)],
      function(col) factor(col, levels = unique(col)))
  }
  
  
  if (is.null(mainFactor.level.order)) {
    mainFactor.level.order <- unique(x[,1])
    calibrartor <- x[,1][1]
  } else if (any(is.na(match(unique(x[,1]), mainFactor.level.order))) == TRUE){
    stop("The `mainFactor.level.order` doesn't match main factor levels.")
  } else {
    x <- x[order(match(x[,1], mainFactor.level.order)), ]
    calibrartor <- x[,1][1]
  }
  
  x <- compute_wDCt(x, numOfFactors, numberOfrefGenes, block)
  
  x[] <- lapply(x, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  # convert the first numOfFactors columns to factor
  if (is.null(block)) {
    x[seq_len(numOfFactors)] <- lapply(
      x[seq_len(numOfFactors)],
      function(col) factor(col, levels = unique(col)))
  } else {
    x[seq_len(numOfFactors + 1)] <- lapply(
      x[seq_len(numOfFactors + 1)],
      function(col) factor(col, levels = unique(col)))
  }
  
  factors <- colnames(x)[1:numOfFactors]
  
  
  # Analysis of variance (no blocking factor)
  if (is.null(block)) {
    # ANOVA
    if(analysisType == "anova") {
      formula_ANOVA <- as.formula(
        paste("wDCt ~", paste(factors, collapse = " * "))
      )
      lm <- lm(formula_ANOVA, data = x)
      lm_formula <- formula(lm)
      ANOVA <- stats::anova(lm)
    }
    # ANCOVA (other factors as covariates)
    if (analysisType == "ancova") {
      formula_ANCOVA <- as.formula(
        paste("wDCt ~", paste(rev(factors), collapse = " + "))
      )
      lm <- lm(formula_ANCOVA, data = x)
      lm_formula <- formula(lm)
      ANOVA <- stats::anova(lm)
    }
    # Repeated measure
    if (analysisType == "repeated") {
        id <- colnames(x)[numOfFactors + 1]
        formula <- as.formula(
          paste("wDCt ~ ", paste(factors, collapse = " * "), " + (1 |", id, ")")
        )
        lm <- suppressMessages(lmerTest::lmer(formula, data = x))
        ANOVA <- stats::anova(lm)
        lm_formula <- formula(lm)
    }
  }
  
  
  # Analysis of variance (blocking factor)
  if (!is.null(block)) {
    # ANOVA
    if(analysisType == "anova") {
      formula_ANOVA <- as.formula(
        paste("wDCt ~", block, "+", paste(factors, collapse = " * "))
      )
      lm <- lm(formula_ANOVA, data = x)
      lm_formula <- formula(lm)
      ANOVA <- stats::anova(lm)
    }
    # ANCOVA with blocking factor
    if (analysisType == "ancova") {
      formula_ANCOVA <- as.formula(
        paste("wDCt ~", block, "+", paste(rev(factors), collapse = " + "))
      )
      lm <- lm(formula_ANCOVA, data = x)
      lm_formula <- formula(lm)
      ANOVA <- stats::anova(lm)
    }
    # Repeated measure
    if (analysisType == "repeated") {
        id <- colnames(x)[numOfFactors + 2]
        formula <- as.formula(
          paste("wDCt ~ ", paste(factors, collapse = " * "), " + (1 |", id, ") + (1 |", block,")")
        )
      lm <- suppressMessages(lmerTest::lmer(formula, data = x))
      ANOVA <- stats::anova(lm)
      lm_formula <- formula(lm)
    }
  }
  
  
  pp1 <- suppressMessages(emmeans(lm, colnames(x)[1], data = x, adjust = p.adj, mode = "satterthwaite"))
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  pp3 <- pp2[1:length(mainFactor.level.order)-1,]
  ci  <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(unique(x[,1]))-1,]
  pp  <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)


  bwDCt <- x$wDCt
  se <- summarise(
    group_by(data.frame(factor = x[,1], bwDCt = bwDCt), x[,1]),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))


  sig <- .convert_to_character(pp$p.value)
  contrast <- pp$contrast
  post_hoc_test <- data.frame(contrast,
                              ddCt = - pp$estimate,
                              RE = 1/(2^-(pp$estimate)),
                              log2FC = log2(1/(2^-(pp$estimate))),
                              pvalue = pp$p.value,
                              sig = sig,
                              LCL = 1/(2^-pp$lower.CL),
                              UCL = 1/(2^-pp$upper.CL),
                              se = se$se[-1])

  reference <- data.frame(contrast = mainFactor.level.order[1],
                          ddCt = 0,
                          RE = 1,
                          log2FC = 0,
                          pvalue = 1,
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])

  tableC <- rbind(reference, post_hoc_test)

  tableC$contrast <- as.character(tableC$contrast)
  tableC$contrast <- sapply(strsplit(tableC$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))



  contrast <- tableC$contrast
  LCL <- tableC$LCL
  UCL <- tableC$UCL
  REp <- as.numeric(tableC$RE)
  FCp <- as.numeric(tableC$log2FC)
  significance <- tableC$sig
  se <- tableC$se


  tableC <- data.frame(tableC,
                       Lower.se.RE = 2^(log2(tableC$RE) - tableC$se),
                       Upper.se.RE = 2^(log2(tableC$RE) + tableC$se),
                       Lower.se.log2FC = 0,
                       Upper.se.log2FC = 0)

  for (i in 1:length(tableC$RE)) {
    if (tableC$RE[i] < 1) {
      tableC$Lower.se.log2FC[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      tableC$Upper.se.log2FC[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
    } else {
      tableC$Lower.se.log2FC[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      tableC$Upper.se.log2FC[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
    }
  }



  #round tableC to 5 decimal places
  tableC <- tableC %>%
    mutate_if(is.numeric, ~ round(., 5))
  
  
  # Type of analysis: ancova or anova
    outlist2 <- list(Final_data = x,
                     lm = lm,
                     ANOVA_table = ANOVA,
                     Fold_Change  = tableC,
                     lm_formula = lm_formula)

  invisible(outlist2)
}
