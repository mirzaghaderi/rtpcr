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
#' @param data
#' A data frame containing expression results, typically obtained from
#' \code{ANOVA_DDCt()} or \code{ANOVA_DCt()}.
#'
#' @param show.groupingLetters
#' Logical; if \code{TRUE}, grouping letters from statistical comparisons
#' are added to the bars.
#'
#' @param x_col
#' Integer specifying the column number used for the x-axis (factor levels).
#'
#' @param y_col
#' Integer specifying the column number used for the y-axis (relative expression or fold change).
#'
#' @param Lower.se_col
#' Integer specifying the column number used for the lower limit of the error bar.
#'
#' @param Upper.se_col
#' Integer specifying the column number used for the upper limit of the error bar.
#'
#' @param letters_col
#' Integer specifying the column number containing grouping letters from
#' statistical comparisons.
#'
#' @return
#' A \code{ggplot} object showing a bar plot of relative expression values
#' with error bars and optional significance annotations.
#'
#' @examples
#'
#' # Extract expression results before plotting
#' res <- ANOVA_DCt(
#'   data_1factor,
#'   numberOfrefGenes = 1,
#'   block = NULL
#' )$Result
#'
#' # Generate bar plot
#' plotOneFactor(
#'   res,
#'   x_col = 1,
#'   y_col = 2,
#'   Lower.se_col = 7,
#'   Upper.se_col = 8,
#'   letters_col = 11,
#'   show.groupingLetters = TRUE
#' )



plotOneFactor <- function(data, 
                          x_col, 
                          y_col,
                          Lower.se_col,
                          Upper.se_col,
                          letters_col = NULL,
                          show.groupingLetters = TRUE) {

  
# column names from index
  x_name  <- names(data)[x_col]
  y_name  <- names(data)[y_col]
  lower   <- names(data)[Lower.se_col]
  upper   <- names(data)[Upper.se_col]
  
# compute ymin and ymax BEFORE ggplot
  data$ymin <- ifelse(
    data[[y_name]] < 0,
    data[[lower]],   # negative bars
    data[[lower]]    # positive bars
  )
  
  data$ymax <- ifelse(
    data[[y_name]] < 0,
    data[[upper]],   # negative bars
    data[[upper]]    # positive bars
  )
  
  # Base plot
  p <- ggplot(data, aes(x = .data[[x_name]], y = .data[[y_name]])) +
    geom_col() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1)
  
  # Add grouping letters
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
            .data[[lower]] - 0.2,   # negative bars
            .data[[upper]] + 0.2   # positive bar
          )
        )
      )
  }
  
  return(p)
}
