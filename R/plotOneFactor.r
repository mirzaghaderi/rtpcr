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
#' data <- read.csv(system.file("extdata", "data_1factor.csv", package = "rtpcr"))
#' res <- ANOVA_DCt(
#'   data,
#'   NumOfFactors = 1,
#'   numberOfrefGenes = 1,
#'   block = NULL
#' )
#' 
#' df <- res$combinedResults
#' 
#' plotOneFactor(
#'   df,
#'   x_col = "SA",
#'   y_col = "RE",
#'   Lower.se_col = "Lower.se.RE",
#'   Upper.se_col = "Upper.se.RE",
#'   letters_col = "sig",
#'   letters_d = 0.1,
#'   col_width = 0.7,
#'   err_width = 0.15,
#'   fill_colors = "skyblue",
#'   alpha = 1,
#'   base_size = 16)
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
  
  # checks
  required_cols <- c(x_col, y_col, Lower.se_col, Upper.se_col)
  
  if (!all(required_cols %in% colnames(data))) {
    stop("One or more specified columns do not exist in `data`.")
  }
  
  if (!is.null(letters_col) && !letters_col %in% colnames(data)) {
    stop("`letters_col` does not exist in `data`.")
  }
  
  # error bar columns
  data$ymin <- data[[Lower.se_col]]
  data$ymax <- data[[Upper.se_col]]
  
  # letters
  if (!is.null(letters_col)) {
    data[[letters_col]] <- as.character(data[[letters_col]])
  }
  
  # base plot
  p <- ggplot(
    data,
    aes(
      x = .data[[x_col]],
      y = .data[[y_col]]
    )
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
  
  # add letters
  if (!is.null(letters_col)) {
    p <- p +
      geom_text(
        aes(
          label = .data[[letters_col]],
          y = ifelse(
            .data[[y_col]] < 0,
            ymin - letters_d,
            ymax + letters_d
          )
        )
      )
  }
  
  # final styling
  p +
    .theme_pub(
      base_size       = base_size,
      legend_position = legend_position
    )
}
