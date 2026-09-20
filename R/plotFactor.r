#' @title Bar plot layout for 1- to 4-factor experiments
#'
#' @description
#' \code{plotFactor()} produces bar plot layout for 1- to 4-factor experiments.
#'
#' \code{split_col} is optional. If it is \code{NULL} (the default), no
#' splitting is done.
#'
#' @export
#'
#' @import ggplot2
#'
#' @param data Data frame containing expression results.
#' @param split_col Character. Column name for the 4th (splitting) factor
#'   (optional, default \code{NULL}) -- one \code{plotFactor()} panel is produced
#'   per level of this column, and the panels are combined with
#'   \code{multiplot()}. If \code{NULL}, the function behaves like
#'   \code{plotFactor()} (1-3 factors) and \code{split_levels},
#'   \code{split_titles} and \code{cols} are not used.
#' @param x_col Character. Column name for x-axis (1st factor).
#' @param y_col Character. Column name for bar height.
#' @param Lower.se_col Character. Column name for lower SE.
#' @param Upper.se_col Character. Column name for upper SE.
#' @param group_col Character. Column name for grouping bars (2nd factor, optional).
#' @param facet_col Character. Column name for faceting within each panel (3rd factor, optional).
#' @param facet_ncol Integer. Number of columns in the \code{facet_wrap()} layout of
#'   \code{facet_col}'s levels, within each \code{split_col} panel (optional; see
#'   \code{plotFactor()}).
#' @param facet_nrow Integer. Number of rows in that same \code{facet_wrap()} layout
#'   (optional; see \code{facet_ncol}).
#' @param letters_col Character. Column name for significance letters (optional).
#' @param letters_d Numeric. Vertical offset for letters (default \code{0.2}).
#' @param col_width Numeric. Width of bars (default \code{0.8}).
#' @param err_width Numeric. Width of error bars (default \code{0.15}).
#' @param dodge_width Numeric. Width of dodge for grouped bars (default \code{0.8}).
#' @param fill_colors Optional vector of fill colors to change the default colors.
#' @param color Optional color for the bar outline.
#' @param alpha Numeric. Transparency of bars (default \code{1}).
#' @param base_size Numeric. Base font size for theme (default \code{12}).
#' @param legend_position Character or numeric vector. Legend position (default \code{"right"}).
#' @param removeCalibratorCols Logical. Passed through to \code{plotFactor()}.
#' @param removeCalibratorText Logical. Passed through to \code{plotFactor()}.
#' @param split_levels Optional character vector giving which levels of
#'   \code{split_col} to plot and the order (left-to-right, top-to-bottom) they
#'   appear in the \code{multiplot()} layout. Defaults to all levels, in the
#'   order they first appear in \code{data}. Ignored (with a warning) when
#'   \code{split_col} is \code{NULL}.
#' @param split_titles Logical. If \code{TRUE} (default), each panel gets a
#'   \code{ggtitle()} of \code{"<split_col>: <level>"} so panels stay
#'   identifiable once combined.
#' @param cols Integer. Number of columns in the \code{multiplot()} layout
#'   (default \code{2}).
#' @param ... Additional ggplot2 layer arguments, forwarded to \code{plotFactor()}.
#'
#' @return
#' If \code{split_col} is \code{NULL}, a single \code{ggplot2} object, exactly as
#' returned by \code{plotFactor()}. Otherwise, invisibly returns a named list of the individual \code{ggplot2} objects (one
#' per level of \code{split_col}, in the order plotted), after drawing them
#' combined via \code{multiplot()}. Access an individual panel later with
#' e.g. \code{result[["<level>"]]}.
#'
#' @examples
#' # 3-factor experiment: no `split_col`.
#' data3 <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
#'
#' res3 <- ANOVA_DCt(
#'   data3,
#'   numOfFactors = 3,
#'   numberOfrefGenes = 1,
#'   block = NULL)
#'   
#'
#' p <- plotFactor(
#'   res3$relativeExpression,
#'   x_col = "SA",
#'   y_col = "log2FC",
#'   group_col = "Type",
#'   facet_col = "Conc",
#'   facet_ncol = 3,
#'   Lower.se_col = "Lower.se.log2FC",
#'   Upper.se_col = "Upper.se.log2FC",
#'   letters_col = "sig",
#'   letters_d = 0.3,
#'   col_width = 0.7,
#'   dodge_width = 0.7,
#'   base_size = 14)
#' p
#' 
#'
#' # A data_3factor_Multi_Target data example with data splitting. 
#' data2 <- read.csv(system.file("extdata", "data_3factorMultiTarget.csv", package = "rtpcr"))
#' #Perform analysis first
#' res <- ANOVA_DDCt(
#'   data2,
#'   numOfFactors = 3,
#'   numberOfrefGenes = 1,
#'   specs = "stress|genotype*stage",
#'   block = NULL)
#' df <- res$relativeExpression
#' 
#' plotFactor(df,
#'   split_col = "gene",
#'   x_col = "contrast",
#'   y_col = "RE",
#'   Lower.se_col = "Lower.se.RE",
#'   Upper.se_col = "Upper.se.RE",
#'   group_col = "genotype",
#'   facet_col = "stage",
#'   facet_ncol = 4,
#'   letters_col = "sig",
#'   letters_d = 0.2,
#'   col_width = 0.8,
#'   err_width = 0.15,
#'   dodge_width = 0.8,
#'   fill_colors = c("green", "blue"),
#'   color = NA,
#'   alpha = 1,
#'   base_size = 12,
#'   legend_position = "none",
#'   removeCalibratorCols = FALSE,
#'   removeCalibratorText = TRUE,
#'   split_levels = NULL,
#'   split_titles = TRUE,
#'   cols = 4)
plotFactor <- function(data,
                         split_col = NULL,
                         x_col,
                         y_col,
                         Lower.se_col,
                         Upper.se_col,
                         group_col = NULL,
                         facet_col = NULL,
                         facet_ncol = NULL,
                         facet_nrow = NULL,
                         letters_col = NULL,
                         letters_d = 0.2,
                         col_width = 0.8,
                         err_width = 0.15,
                         dodge_width = 0.8,
                         fill_colors = NULL,
                         color = NA,
                         alpha = 1,
                         base_size = 12,
                         legend_position = "right",
                         removeCalibratorCols = FALSE,
                         removeCalibratorText = FALSE,
                         split_levels = NULL,
                         split_titles = TRUE,
                         cols = 2,
                         ...) {

  # Calls plotFactor() on a (sub)set of the data with all the shared arguments
  make_panel <- function(d, ...) {
    .plotF(
      data = d,
      x_col = x_col,
      y_col = y_col,
      Lower.se_col = Lower.se_col,
      Upper.se_col = Upper.se_col,
      group_col = group_col,
      facet_col = facet_col,
      facet_ncol = facet_ncol,
      facet_nrow = facet_nrow,
      letters_col = letters_col,
      letters_d = letters_d,
      col_width = col_width,
      err_width = err_width,
      dodge_width = dodge_width,
      fill_colors = fill_colors,
      color = color,
      alpha = alpha,
      base_size = base_size,
      legend_position = legend_position,
      removeCalibratorCols = removeCalibratorCols,
      removeCalibratorText = removeCalibratorText,
      ...
    )
  }

  # No 4th factor: behave exactly like plotFactor() (1-3 factors) and return
  # the single ggplot object.
  if (is.null(split_col)) {
    if (!is.null(split_levels)) {
      warning("`split_levels` is ignored because `split_col` is NULL.")
    }
    return(make_panel(data, ...))
  }

  if (!is.character(split_col) || length(split_col) != 1L) {
    stop("`split_col` must be NULL or a single column name.")
  }

  if (!split_col %in% colnames(data)) {
    stop("`split_col` does not exist in `data`.")
  }

  split_values <- as.character(data[[split_col]])

  if (is.null(split_levels)) {
    # Preserve first-appearance order, matching plotFactor()'s own
    # factor(x, levels = unique(x)) convention for character columns.
    split_levels <- unique(split_values)
  } else {
    missing_levels <- setdiff(split_levels, unique(split_values))
    if (length(missing_levels) > 0) {
      stop("The following `split_levels` were not found in `data[[split_col]]`: ",
           paste(missing_levels, collapse = ", "))
    }
  }

  plots <- setNames(vector("list", length(split_levels)), split_levels)

  for (lev in split_levels) {
    sub_data <- data[split_values == lev, , drop = FALSE]

    p <- make_panel(sub_data, ...)

    if (split_titles) {
      p <- p + ggtitle(lev)
    }

    plots[[lev]] <- p
  }

  do.call(multiplot, c(unname(plots), list(cols = cols)))

  invisible(plots)
}
