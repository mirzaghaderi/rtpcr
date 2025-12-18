#' @title Bar plot of gene expression for single-factor experiments
#'
#' @description
#' Creates a bar plot of relative gene expression (fold change) values
#' from a single-factor experiment, including standard error or
#' confidence interval error bars and statistical significance.
#'
#' @details
#' The \code{plotOneFactor} function generates a bar plot of fold change
#' values for target genes using the output expression tables produced by
#' functions such as \code{ANOVA_DDCt()} or \code{ANOVA_DCt()}.
#' Error bars can represent either standard error (SE) or 95\% confidence
#' intervals, and optional grouping letters from post hoc statistical
#' comparisons can be displayed above the bars.
#'
#' @author
#' Ghader Mirzaghaderi
#'
#' @export
#'
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#'
#' @param data Data frame
#' @param x_col Numeric. Column index for x-axis
#' @param y_col Numeric. Column index for bar height
#' @param Lower.se_col Numeric. Column index for lower SE
#' @param Upper.se_col Numeric. Column index for upper SE
#' @param letters_col Optional column index for grouping letters
#' @param letters_d Numeric. Vertical offset for letters (default 0.2)
#' @param dodge_width Numeric. Spacing for bars (default 0.8)
#' @param col_width Numeric. Width of bars (default 0.8)
#' @param err_width Numeric. Width of error bars (default 0.15)
#' @param fill_colors Optional vector of fill colors
#' @param alpha Numeric. Transparency of bars (default 1)
#' @param base_size Numeric. Base font size for theme (default 12)
#' @param legend_position Character. Legend position (default "right")
#' @param ... Additional valid ggplot2 layer arguments
#' @return ggplot2 plot object
#' @export
plotOneFactor <- function(data,
                          x_col,
                          y_col,
                          Lower.se_col,
                          Upper.se_col,
                          letters_col = NULL,
                          letters_d = 0.2,
                          dodge_width = 0.8,
                          col_width = 0.8,
                          err_width = 0.15,
                          fill_colors = NULL,
                          alpha = 1,
                          base_size = 12,
                          legend_position = "right",
                          ...) {
  
  x_name <- names(data)[x_col]
  y_name <- names(data)[y_col]
  lower  <- names(data)[Lower.se_col]
  upper  <- names(data)[Upper.se_col]
  
  # Precompute ymin/ymax
  data$ymin <- data[[lower]]
  data$ymax <- data[[upper]]
  
  # Convert letters to character if provided
  if (!is.null(letters_col)) {
    letters_name <- names(data)[letters_col]
    data[[letters_name]] <- as.character(data[[letters_name]])
  }
  
  p <- ggplot(data, aes(x = .data[[x_name]], y = .data[[y_name]])) +
    geom_pub_cols(
      col_width = col_width,
      err_width = err_width,
      fill_colors = fill_colors,
      dodge_width = dodge_width,
      alpha = alpha,
      ...
    )
  
  # Add letters
  if (!is.null(letters_col)) {
    pos <- position_dodge(width = dodge_width)
    p <- p + geom_text(
      aes(
        label = .data[[letters_name]],
        y = ifelse(.data[[y_name]] < 0, ymin - letters_d, ymax + letters_d)
      ),
      position = pos,
      ...
    )
  }
  p + theme_pub(base_size = base_size, legend_position = legend_position) 
}