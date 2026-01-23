#' @title Bar plot of gene expression for 1-, 2-, or 3-factor experiments
#'
#' @description
#' Creates a bar plot of relative gene expression (fold change) values
#' from 1-, 2-, or 3-factor experiments, including error bars and statistical
#' significance annotations.
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
#' @param data Data frame containing expression results
#' @param x_col Character. Column name for x-axis
#' @param y_col Character. Column name for bar height
#' @param Lower.se_col Character. Column name for lower SE
#' @param Upper.se_col Character. Column name for upper SE
#' @param group_col Character. Column name for grouping bars (optional)
#' @param facet_col Character. Column name for faceting (optional)
#' @param letters_col Character. Column name for significance letters (optional)
#' @param letters_d Numeric. Vertical offset for letters (default \code{0.2})
#' @param col_width Numeric. Width of bars (default \code{0.8})
#' @param err_width Numeric. Width of error bars (default \code{0.15})
#' @param dodge_width Numeric. Width of dodge for grouped bars (default \code{0.8})
#' @param fill_colors Optional vector of fill colors to change the default colors
#' @param color Optional color for the bar outline
#' @param alpha Numeric. Transparency of bars (default \code{1})
#' @param base_size Numeric. Base font size for theme (default \code{12})
#' @param legend_position Character or numeric vector. Legend position (default \code{right})
#' @param ... Additional ggplot2 layer arguments
#' 
#' @return ggplot2 plot object
#' 
#' @examples
#' data <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
#' 
#' res <- ANOVA_DDCt(x = data,
#'   numOfFactors = 2,
#'   numberOfrefGenes = 2,
#'   block = "block",
#'   mainFactor.column = 2,
#'   p.adj = "none")
#' 
#' df <- res$relativeExpression
#' 
#' p1 <- plotFactor(
#'   data = df,
#'   x_col = "contrast",
#'   y_col = "RE",
#'   group_col = "gene",
#'   facet_col = "gene",
#'   Lower.se_col = "Lower.se.RE",
#'   Upper.se_col = "Upper.se.RE",
#'   letters_col = "sig",
#'   letters_d = 0.2,
#'   alpha = 1,
#'   col_width = 0.7,
#'   dodge_width = 0.7,
#'   base_size = 14, 
#'   legend_position = "none")
#' 
#' p1
#' 
#' 
#' data2 <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
#' #Perform analysis first
#' res <- ANOVA_DCt(
#'   data2,
#'   numOfFactors = 3,
#'   numberOfrefGenes = 1,
#'   block = NULL)
#'   
#' df <- res$relativeExpression
#' # Generate three-factor bar plot
#'  p <- plotFactor(
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
#'   #fill_colors = c("blue", "brown"),
#'   color = "black",
#'   base_size = 14, 
#'   alpha = 1,
#'   legend_position = c(0.1, 0.2))
#' p
#'   
plotFactor <- function(data,
                       x_col,
                       y_col,
                       Lower.se_col,
                       Upper.se_col,
                       group_col = NULL,
                       facet_col = NULL,
                       letters_col = NULL,
                       letters_d = 0.2,
                       col_width = 0.8,
                       err_width = 0.15,
                       dodge_width = 0.8,
                       fill_colors = NULL,
                       color = "black",
                       alpha = 1,
                       base_size = 12,
                       legend_position = "right",
                       ...) {
  
  # required columns
  required_cols <- c(x_col, y_col, Lower.se_col, Upper.se_col)
  if (!is.null(group_col)) required_cols <- c(required_cols, group_col)
  if (!is.null(facet_col)) required_cols <- c(required_cols, facet_col)
  
  if (!all(required_cols %in% colnames(data))) {
    stop("One or more specified columns do not exist in `data`.")
  }
  
  if (!is.null(letters_col) && !letters_col %in% colnames(data)) {
    stop("`letters_col` does not exist in `data`.")
  }
  
  # add error columns
  data$ymin <- data[[Lower.se_col]]
  data$ymax <- data[[Upper.se_col]]
  
  if (!is.null(letters_col)) {
    data[[letters_col]] <- as.character(data[[letters_col]])
  }
  
  
  # To factor:
  data[] <- lapply(data, function(data) {
    if (is.character(data)) factor(data, levels = unique(data)) else data
  })

  
  # 1-factor plot
  if (is.null(group_col) && is.null(facet_col)) {
    p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
      geom_col(width = col_width, fill = fill_colors[1] %||% "grey40", alpha = alpha, color = color, ...)
    p <- p + geom_errorbar(aes(ymin = ymin, ymax = ymax), width = err_width)
    
    if (!is.null(letters_col)) {
      p <- p + geom_text(aes(
        label = .data[[letters_col]],
        y = ifelse(.data[[y_col]] < 0, ymin - letters_d, ymax + letters_d)
      ))
    }
    
  } else {
    # 2- or 3-factor plot
    p <- ggplot(data, aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      fill = .data[[group_col]]
    )) +
      .geom_pub_cols(
        col_width = col_width,
        err_width = err_width,
        fill_colors = fill_colors,
        dodge_width = dodge_width,
        alpha = alpha,
        color = color,
        ...
      )
    
    if (!is.null(facet_col)) {
      p <- p + facet_wrap(vars(.data[[facet_col]]))
    }
    
    if (!is.null(letters_col)) {
      pos <- position_dodge(width = dodge_width)
      p <- p + geom_text(aes(
        label = .data[[letters_col]],
        y = ifelse(.data[[y_col]] < 0, ymin - letters_d, ymax + letters_d)
      ), position = pos, ...)
    }
  }
  
  if(y_col == "log2FC" && any(data$log2FC < 0)){
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  }
  p + 
    .theme_pub(base_size = base_size, legend_position = legend_position) +
    xlab(NULL) +
    theme(axis.text.x = element_text(size = base_size, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = base_size,color = "black", angle = 0),
          panel.border = element_rect(color = "black"),
          legend.text = element_text(colour = "black", size = base_size),
          legend.background = element_rect(fill = "transparent"),
          strip.background = element_blank(),            # removes the faceting gray background
          strip.text = element_text(size = base_size))   # keeps the text visible
}
