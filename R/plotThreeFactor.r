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
#' @param show.groupingLetters
#' Logical; if \code{TRUE}, grouping letters are displayed above the bars.
#'
#' @return
#' A \code{ggplot} object showing relative expression values arranged
#' by three experimental factors, with error bars and optional
#' significance annotations.
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
#'   show.groupingLetters = TRUE
#' )



plotThreeFactor <- function(data, 
                            x_col,        # x-axis factor
                            y_col,        # bar height
                            group_col,    # fill groups
                            facet_col,    # facet grid
                            Lower.se_col, # lower SE column
                            Upper.se_col, # upper SE column
                            letters_col = NULL,
                            show.groupingLetters = TRUE) {
  
  # Extract column names
  x_name     <- names(data)[x_col]
  y_name     <- names(data)[y_col]
  group_name <- names(data)[group_col]
  facet_name <- names(data)[facet_col]
  lower_name <- names(data)[Lower.se_col]
  upper_name <- names(data)[Upper.se_col]
  
  # compute ymin and ymax BEFORE ggplot
  data$ymin <- data[[lower_name]]
  data$ymax <- data[[upper_name]]


  
  # Base plot
  p <- ggplot(data,
              aes(x = .data[[x_name]],
                  y = .data[[y_name]],
                  fill = .data[[group_name]])) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.15,
                  position = position_dodge(width = 0.9)) +
    facet_wrap(vars(.data[[facet_name]])) +
    theme_bw() +
    labs(x = x_name,
         y = y_name,
         fill = group_name)
  
  # ---- Optional grouping letters ----
  if (show.groupingLetters) {
    if (is.null(letters_col)) {
      stop("letters_col must be provided when show.groupingLetters = TRUE")
    }
    
    letters_name <- names(data)[letters_col]
    
    p <- p +
      geom_text(
        aes(
          label = .data[[letters_name]],
          y = ifelse(
            .data[[y_name]] < 0,
            ymin - 0.5,   # negative bars
            ymax + 0.3   # positive bar
          )
        ),
        position = position_dodge(width = 0.9),
        vjust = 0
      )
  }
  
  return(p)
}
