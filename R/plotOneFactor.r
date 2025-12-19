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
#' @param col_width Numeric. Width of bars (default 0.8)
#' @param err_width Numeric. Width of error bars (default 0.15)
#' @param fill_colors Optional vector of fill colors
#' @param alpha Numeric. Transparency of bars (default 1)
#' @param base_size Numeric. Base font size for theme (default 12)
#' @param legend_position Character. Legend position (default "right")
#' @param ... Additional valid ggplot2 layer arguments
#' @return ggplot2 plot object
#' 
#' @examples
#' 
#' res <- ANOVA_DCt(
#' data_1factor,
#' NumOfFactors = 1,
#' numberOfrefGenes = 1,
#' block = NULL
#' )
#' 
#' df <- res$combinedResults
#' 
#' plotOneFactor(
#'   df,
#'   x_col = 1,
#'   y_col = 2,
#'   Lower.se_col = 7,
#'   Upper.se_col = 8,
#'   letters_col = 11,
#'   letters_d = 0.1,
#'   col_width = 0.7,
#'   err_width = 0.15,
#'   fill_colors = "skyblue",
#'   alpha = 1,
#'   base_size = 16
#' )
#' 
#'
#' 
plotOneFactor <- function(data,
                          x_col,
                          y_col,
                          Lower.se_col,
                          Upper.se_col,
                          letters_col = NULL,
                          letters_d = 0.2,
                          col_width = 0.8,
                          err_width = 0.15,
                          fill_colors = "grey40",
                          alpha = 1,
                          base_size = 12,
                          legend_position = "none",
                          ...) {
  
  # Column names
  x_name <- names(data)[x_col]
  y_name <- names(data)[y_col]
  lower  <- names(data)[Lower.se_col]
  upper  <- names(data)[Upper.se_col]
  
  # preserve x order as in data
  # data[[x_name]] <- factor(data[[x_name]], levels = unique(data[[x_name]]))
  
  # Precompute error limits
  data$ymin <- data[[lower]]
  data$ymax <- data[[upper]]
  
  # Letters as character
  if (!is.null(letters_col)) {
    letters_name <- names(data)[letters_col]
    data[[letters_name]] <- as.character(data[[letters_name]])
  }
  
  # Base plot
  p <- ggplot(
    data,
    aes(x = .data[[x_name]], y = .data[[y_name]])
  ) +
    geom_col(
      width = col_width,
      fill  = fill_colors[1],  
      alpha = alpha,
      ...
    ) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = err_width
    )
  
  # Add letters
  if (!is.null(letters_col)) {
    p <- p +
      geom_text(
        aes(
          label = .data[[letters_name]],
          y = ifelse(
            .data[[y_name]] < 0,
            ymin - letters_d,
            ymax + letters_d
          )
        )
      )
  }
  
  # Final styling
  p +
    .theme_pub() 
}
