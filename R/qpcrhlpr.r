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
        stop("Numeric column found among factor columns. Please check input.")
      }
      factor(z)
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





#' Publication-ready column + error-bar layers
#'
#' Adds geom_col and geom_errorbar layers with optional fill colors.
#'
#' @param col_width Numeric. Width of bars (default 0.8)
#' @param err_width Numeric. Width of error bars (default 0.15)
#' @param fill_colors Optional vector of fill colors
#' @param dodge_width Numeric. Horizontal spacing between groups (default 0.8)
#' @param alpha Numeric. Transparency of bars (default 1)
#' @param ... Additional valid ggplot2 geom_col arguments
#' @return A list of ggplot2 layers
#' @export
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

#' Publication-ready ggplot2 theme
#' @param base_size Numeric. Base font size (default 12)
#' @param base_family Character. Font family (default "sans")
#' @param legend_position Character or numeric vector. Position of legend (default "right")
#' @param ... Additional theme arguments
#' @return ggplot2 theme object
#' @export
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

#' Wrapper for .geom_pub_cols
#' @param ... Additional theme arguments
#' @export
geom_pub_cols <- function(...) {
  .geom_pub_cols(...)
}

#' Wrapper for .theme_pub
#' @param ... Additional theme arguments
#' @export
theme_pub <- function(...) {
  .theme_pub(...)
}
