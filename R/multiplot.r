#' @title Multiple plot function
#' @description \code{multiplot} function combines multiple ggplot objects into a single plate.
#' @details Combining multiple ggplot objects into a single plate.
#' @author gist.github.com/pedroj/ffe89c67282f82c1813d
#' @export multiplot
#' @import grid
#' @param ... ggplot objects can be passed in ... or to plotlist (as a list of ggplot objects)
#' @param cols Number of columns in the panel
#' @return A  multiple-plots plate
#' @examples
#' 
#' p1 <- qpcrTTESTplot(data_ttest, 
#'                     numberOfrefGenes = 1,
#'                     ylab = "Average Fold Change (FC)",
#'                     width = 0.3)
#' 
#' 
#' out2 <- qpcrANOVARE(data_1factor, numberOfrefGenes = 1)$Result
#' p2 <- oneFACTORplot(out2,
#'                     width = 0.2,
#'                     fill = "skyblue",
#'                     y.axis.adjust = 0.5,
#'                     y.axis.by = 1,
#'                     errorbar = "ci",
#'                     show.letters = TRUE,
#'                     letter.position.adjust = 0.1,
#'                     ylab = "Relative Expression (RE)",
#'                     xlab = "Factor Levels",
#'                     fontsize = 12)
#'                     
#' multiplot(p1, p2, cols=2)
#' 
#' multiplot(p1, p2, cols=1)
#'
#'
#'
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
