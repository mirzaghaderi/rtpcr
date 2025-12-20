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
#' data <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
#' #Perform analysis first
#' res <- ANOVA_DCt(
#'   data,
#'   NumOfFactors = 3,
#'   numberOfrefGenes = 1,
#'   block = NULL)
#'   
#' df <- res$combinedResults
#'  df
#'  # Generate three-factor bar plot
#'  p <- plotThreeFactor(
#'   df,
#'   x_col = "SA",       
#'   y_col = "log2FC",       
#'   group_col = "Type",   
#'   facet_col = "Conc",    
#'   Lower.se_col = "Lower.se.log2FC",
#'   Upper.se_col = "Upper.se.log2FC",
#'   letters_col = "sig",
#'   letters_d = 0.3,
#'   col_width = 0.7, 
#'   dodge_width = 0.7,
#'   fill_colors = c("blue", "brown"),
#'   base_size = 14, 
#'   alpha = 1,
#'   legend_position = c(0.1, 0.2))
#' p
#' library(ggplot2)
#' p + theme(
#'   panel.border = element_rect(color = "black", linewidth = 0.5))




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
  
  # checks
  required_cols <- c(
    x_col, y_col, group_col, facet_col,
    Lower.se_col, Upper.se_col
  )
  
  if (!all(required_cols %in% colnames(data))) {
    stop("One or more specified columns do not exist in `data`.")
  }
  
  if (!is.null(letters_col) && !letters_col %in% colnames(data)) {
    stop("`letters_col` does not exist in `data`.")
  }
  
  # error bar columns
  data$ymin <- data[[Lower.se_col]]
  data$ymax <- data[[Upper.se_col]]
  
  if (!is.null(letters_col)) {
    data[[letters_col]] <- as.character(data[[letters_col]])
  }
  
  # plot
  p <- ggplot(
    data,
    aes(
      x    = .data[[x_col]],
      y    = .data[[y_col]],
      fill = .data[[group_col]]
    )
  ) +
    .geom_pub_cols(
      col_width   = col_width,
      err_width   = err_width,
      fill_colors = fill_colors,
      dodge_width = dodge_width,
      alpha       = alpha,
      ...
    ) +
    facet_wrap(vars(.data[[facet_col]]))
  
  # letters
  if (!is.null(letters_col)) {
    pos <- position_dodge(width = dodge_width)
    
    p <- p + geom_text(
      aes(
        label = .data[[letters_col]],
        y = ifelse(
          .data[[y_col]] < 0,
          ymin - letters_d,
          ymax + letters_d
        )
      ),
      position = pos,
      ...
    )
  }
  
  p + .theme_pub(
    base_size       = base_size,
    legend_position = legend_position
  )
}
