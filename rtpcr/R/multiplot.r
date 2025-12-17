#' @title Combine multiple ggplot objects into a single layout
#'
#' @description
#' The \code{multiplot} function arranges multiple \code{ggplot2} objects
#' into a single plotting layout with a specified number of columns.
#'
#' @details
#' Multiple \code{ggplot2} objects can be provided either as separate
#' arguments via \code{...}.
#' The function uses the \code{grid} package to control the layout.
#'
#' @author
#' Pedro J. (adapted from \url{https://gist.github.com/pedroj/ffe89c67282f82c1813d})
#'
#' @export
#'
#' @import grid
#'
#' @param ...
#' One or more \code{ggplot2} objects.
#'
#'
#' @param cols
#' Integer specifying the number of columns in the layout.
#'
#' @return
#' A grid object displaying multiple plots arranged in the specified layout.
#'
#' @examples
#'
#' # Example using output from TTEST_DDCt
#' a <- TTEST_DDCt(
#'   data_1factor_one_ref,
#'   numberOfrefGenes = 1,
#'   plotType = "log2FC"
#' )
#' p1 <- a$plot
#'
#' # Example using output from ANOVA_DCt
#' out2 <- ANOVA_DCt(
#'   data_1factor,
#'   numberOfrefGenes = 1,
#'   block = NULL
#' )$Result
#'
#' p2 <- plotOneFactor(
#'   out2,
#'   x_col = 1,
#'   y_col = 2,
#'   Lower.se_col = 7,
#'   Upper.se_col = 8,
#'   letters_col = 11,
#'   show.groupingLetters = TRUE
#' )
#'
#' # Combine plots into a single layout
#' multiplot(p1, p2, cols = 2)
#'
#' multiplot(p1, p2, cols = 1)



multiplot <- function(..., cols=1) {
  
  # Make a list from the ... arguments
  plots <- c(list(...))
  
  numPlots = length(plots)
  
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow = TRUE)
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
