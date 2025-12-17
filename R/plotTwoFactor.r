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
#' by two experimental factors, with error bars and optional
#' significance annotations.
#'
#' @examples
#'
#' # Two-factor ANOVA example
#' res <- ANOVA_DCt(
#'   data_2factorBlock,
#'   block = "Block",
#'   numberOfrefGenes = 1
#' )
#'
#' expr_data <- res$Result
#'
#' plotTwoFactor(
#'   data = expr_data,
#'   x_col = 2,
#'   y_col = 3,
#'   group_col = 1,
#'   Lower.se_col = 8,
#'   Upper.se_col = 9,
#'   letters_col = 12,
#'   show.groupingLetters = TRUE
#' )
#'
#' # Combining fold-change results from two different genes
#' a <- REPEATED_DDCt(
#'   data_repeated_measure_1,
#'   numberOfrefGenes = 1,
#'   factor = "time",
#'   calibratorLevel = "1",
#'   block = NULL,
#'   plot = FALSE
#' )
#'
#' b <- REPEATED_DDCt(
#'   data_repeated_measure_2,
#'   numberOfrefGenes = 1,
#'   factor = "time",
#'   calibratorLevel = "1",
#'   block = NULL,
#'   plot = FALSE
#' )
#'
#' df_combined <- rbind(
#'   a$Relative_Expression_table,
#'   b$Relative_Expression_table
#' )
#'
#' df_combined$gene <- factor(rep(c("Gene1", "Gene2"), each = 3))
#'
#' plotTwoFactor(
#'   data = df_combined,
#'   x_col = 1,
#'   y_col = 2,
#'   group_col = 13,
#'   Lower.se_col = 9,
#'   Upper.se_col = 10,
#'   letters_col = 5,
#'   show.groupingLetters = TRUE
#' )



plotTwoFactor <- function(data, 
                          x_col,        # x-axis factor
                          y_col,        # bar height
                          group_col,    # fill groups
                          Lower.se_col, # lower SE column
                          Upper.se_col, # upper SE column
                          letters_col = NULL,
                          show.groupingLetters = TRUE){
  
  # Extract Column Names by Index
  x_name      <- names(data)[x_col]
  y_name      <- names(data)[y_col]
  group_name  <- names(data)[group_col]
  lower  <- names(data)[Lower.se_col]
  upper  <- names(data)[Upper.se_col]
  letter_name <- if (!is.null(letters_col)) names(data)[letters_col] else NULL
  
  # compute ymin and ymax BEFORE ggplot
  data$ymin <- data[[lower]]
  data$ymax <- data[[upper]]

  
  #Base ggplot
  p <- ggplot(data,
              aes(x = .data[[x_name]],
                  y = .data[[y_name]],
                  fill = .data[[group_name]])) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.8)
  

  
  
  #Error Bars
  p <- p + geom_errorbar(
    aes(ymin = ymin,
        ymax = ymax),
    width = 0.15,
    position = position_dodge(width = 0.8))
  
  #Add Grouping Letters if Provided
  if (show.groupingLetters & !is.null(letters_col)) {

    letters_name <- names(data)[letters_col]
    p <- p +
      geom_text(
        aes(
          label = .data[[letters_name]],
          y = ifelse(
            .data[[y_name]] < 0,
            ymin - 0.3,   # negative bars
            ymax + 0.3   # positive bar
          )
        ), position = position_dodge(width = 0.8)
      )
    
    
  }
  
  return(p)
}

