#' @title Bar plot of gene expression for two-factor experiments
#'
#' @description
#' Creates a bar plot of relative gene expression (fold change) values
#' from a two-factor experiment, including error bars and statistical
#' significance annotations.
#'
#' @details
#' The \code{plotTwoFactor} function generates a bar plot of average
#' fold change (relative expression) values for target genes using
#' expression tables produced by functions such as
#' \code{ANOVA_DDCt()} or \code{ANOVA_DCt()}.
#' One factor is mapped to the x-axis and the second factor is used
#' to group the bars (fill aesthetic). Error bars represent standard
#' error (SE) or 95\% confidence intervals (CI). Optional grouping
#' letters from post hoc statistical comparisons can be displayed.
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
#' @param data
#' A data frame containing expression results, typically obtained from
#' \code{ANOVA_DDCt()} or \code{ANOVA_DCt()}.
#'
#' @param x_col
#' Integer specifying the column number used for the x-axis factor.
#'
#' @param y_col
#' Integer specifying the column number used for the bar height
#' (relative expression or fold change).
#'
#' @param group_col
#' Integer specifying the column number used for grouping bars
#' (fill aesthetic).
#' 
#' @param col_width 
#' Numeric. Width of bars (default 0.8)
#' 
#' @param err_width 
#' Numeric. Width of error bars (default 0.15)
#' 
#' @param fill_colors 
#' Optional vector of fill colors
#' 
#' @param alpha 
#' Numeric. Transparency of bars (default 1)
#' 
#' @param base_size Numeric. Base font size (default 12)
#'
#' @param Lower.se_col
#' Integer specifying the column number used for the lower limit
#' of the error bar.
#'
#' @param Upper.se_col
#' Integer specifying the column number used for the upper limit
#' of the error bar.
#'
#' @param letters_col
#' Integer specifying the column number containing grouping letters
#' from statistical comparisons.
#' 
#' @param letters_d
#' Numeric specifying the distance between sig letters and error bar.
#' 
#' @param dodge_width 
#' Numeric. Width of the dodge position adjustment for grouped bars. 
#' 
#' @param legend_position 
#' Character or numeric vector. Position of legend (default "right")
#' 
#' @param ... 
#' Additional ggplot2 layer arguments (e.g., fill, alpha, color)
#'
#' @return
#' A ggplot object
#' 
#'
#' @examples
#'
#' a <- ANOVA_DCt(data_2factorBlock, block = "Block", numberOfrefGenes = 1)
#' data <- a$Results
#' 
#' p1 <- plotTwoFactor(
#'   data = data,
#'   x_col = 2,
#'   y_col = 3,
#'   group_col = 1,
#'   Lower.se_col = 8,
#'   Upper.se_col = 9,
#'   letters_col = 12,
#'   letters_d = 0.2,
#'   fill_colors = c("aquamarine4", "gold2"),
#'   alpha = 1,
#'   col_width = 0.7,
#'   dodge_width = 0.7,
#'   base_size = 16, 
#'   legend_position = c(0.2, 0.8)
#' )
#' p1
#' 
#' 
#' p2 <- plotTwoFactor(
#'   data = data,
#'   x_col = 2,
#'   y_col = 4,
#'   group_col = 1,
#'   Lower.se_col = 10,
#'   Upper.se_col = 11,
#'   letters_col = 12,
#'   letters_d = 0.2,
#'   fill_colors = c("aquamarine4", "gold2"),
#'   alpha = 1,
#'   col_width = 0.7,
#'   dodge_width = 0.7,
#'   base_size = 16, 
#'   legend_position = c(0.2, 0.8)
#' )
#' p2



plotTwoFactor <- function(data,
                          x_col,
                          y_col,
                          group_col,
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
  group_name <- names(data)[group_col]
  lower <- names(data)[Lower.se_col]
  upper <- names(data)[Upper.se_col]
  
  data$ymin <- data[[lower]]
  data$ymax <- data[[upper]]
  
  if (!is.null(letters_col)) {
    letters_name <- names(data)[letters_col]
    data[[letters_name]] <- as.character(data[[letters_name]])
  }
  
  p <- ggplot(data, aes(x = .data[[x_name]], y = .data[[y_name]], fill = .data[[group_name]])) +
    .geom_pub_cols(
      col_width = col_width,
      err_width = err_width,
      fill_colors = fill_colors,
      dodge_width = dodge_width,
      alpha = alpha,
      ...
    )
  
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
  p + .theme_pub(base_size = base_size, legend_position = legend_position)
}
