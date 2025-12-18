#' @title Bar plot of gene expression for three-factor experiments
#'
#' @description
#' Creates a bar plot of relative gene expression (fold change) values
#' from a three-factor experiment, including error bars and statistical
#' significance annotations.
#'
#' @details
#' The \code{plotThreeFactor} function generates a bar plot of average
#' fold change (relative expression) values for target genes using the
#' output expression tables produced by functions such as
#' \code{ANOVA_DDCt()} or \code{ANOVA_DCt()}.
#' One factor is mapped to the x-axis, one to bar grouping (fill),
#' and one to faceting. Error bars represent standard error (SE)
#' or 95\% confidence intervals, and optional grouping letters from
#' post hoc statistical comparisons can be displayed.
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
#' @param group_col
#' Integer specifying the column number used for grouping bars
#' (fill aesthetic).
#'
#' @param facet_col
#' Integer specifying the column number used for faceting the plot.
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
#' @examples
#'
#' # Perform analysis first
#' res <- ANOVA_DCt(
#'   data_3factor,
#'   numberOfrefGenes = 1,
#'   block = NULL
#' )
#'
#' expr_data <- res$Result
#'
#' # Generate three-factor bar plot
#' plotThreeFactor(
#'   expr_data,
#'   x_col = 3,        # x-axis factor
#'   y_col = 5,        # bar height
#'   group_col = 1,    # grouping (fill)
#'   facet_col = 2,    # faceting factor
#'   Lower.se_col = 11,
#'   Upper.se_col = 12,
#'   letters_col = 13,
#'   letters_d = 0.4,
#'   dodge_width = 0.9,       # controls spacing
#'   fill_colors = c("blue", "brown"),
#'   alpha = 1,
#'   legend_position = c(0.1, 0.2)
#' )



plotThreeFactor <- function(data,
                            x_col,
                            y_col,
                            group_col,
                            facet_col,
                            Lower.se_col,
                            Upper.se_col,
                            letters_col = NULL,
                            letters_d = 0.2,
                            dodge_width = 0.9,
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
  facet_name <- names(data)[facet_col]
  lower <- names(data)[Lower.se_col]
  upper <- names(data)[Upper.se_col]
  
  data$ymin <- data[[lower]]
  data$ymax <- data[[upper]]
  
  if (!is.null(letters_col)) {
    letters_name <- names(data)[letters_col]
    data[[letters_name]] <- as.character(data[[letters_name]])
  }
  
  p <- ggplot(data, aes(x = .data[[x_name]],
                        y = .data[[y_name]],
                        fill = .data[[group_name]])) +
    geom_pub_cols(
      col_width = col_width,
      err_width = err_width,
      fill_colors = fill_colors,
      dodge_width = dodge_width,
      alpha = alpha,
      ...
    ) +
    facet_wrap(vars(.data[[facet_name]]))
  
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
