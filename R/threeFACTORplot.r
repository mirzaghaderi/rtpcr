#' @title Bar plot of the relative gene expression (RE) from a three-factor experiment
#' @description Bar plot of the relative expression (RE) of a gene along with the confidence interval and significance
#' @details The \code{threeFACTORplot} function generates the bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @author Ghader Mirzaghaderi
#' @export threeFACTORplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x a data frame created by \link{qpcrANOVA} function via running \code{res <- qpcrANOVA(x)$Result} on a three factorial data such as \code{datathreegr} example data.
#' @param arrangement arrangement of the grouping columns in the output graph, for example c(2,3,1) or c(1,3,2). This affects grapg output.
#' @param bar.width a positive number determining bar width.
#' @param fill  a color vector specifying the fill color for the columns of the bar plot. One of the palettes in \code{\link[RColorBrewer]{display.brewer.all}} (e.g. "Reds" or "Blues", ...) can be applied.
#' @param xlab  the title of the x axis
#' @param ylab  the title of the y axis
#' @param errorbar Type of error bar, can be \code{std} or \code{ci}.
#' @param y.axis.adjust  a negative or positive number for reducing or increasing the length of the y axis.
#' @param y.axis.by determines y axis step length
#' @param letter.position.adjust adjust the distance between the grouping letters to the error bars
#' @param legend.title legend title
#' @param fontsize all fonts size of the plot
#' @param fontsizePvalue font size of the pvalue labels
#' @param legend.position a two digit vector specifying the legend position.
#' @param show.letters a logical variable. If TRUE, mean grouping letters are added to the bars. 
#' @param axis.text.x.angle angle of x axis text
#' @param axis.text.x.hjust horizontal justification of x axis text
#' @return Bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @examples
#' 
#' #' # See a sample data frame
#' data_3factor
#'
#' # Before plotting, the result needs to be extracted as below:
#' res <- qpcrANOVA(data_3factor, numberOfrefGenes = 1)
#' res
#'
#' # Arrange the first three colunms of the result table.
#' # This determines the columns order and shapes the plot output.
#' threeFACTORplot(res,
#'     arrangement = c(3, 1, 2),
#'     xlab = "condition")
#'
#'
#'
#' threeFACTORplot(res,
#'    arrangement = c(1, 2, 3),
#'    bar.width = 0.5,
#'    fill = "Greys",
#'    xlab = "Genotype",
#'    ylab = "Relative Expression")
#'
#'
#'
#' # Reordering factor levels to a desired order.
#' res$Conc <- factor(res$Conc, levels = c("L","M","H"))
#' res$Type <- factor(res$Type, levels = c("S","R"))
#'
#' # Producing the plot
#' threeFACTORplot(res,
#'    arrangement = c(2, 3, 1),
#'    bar.width = 0.5,
#'    fill = "Reds",
#'    xlab = "Drought",
#'    ylab = "Relative Expression",
#'    errorbar = "std",
#'    legend.title = "Genotype",
#'    legend.position = c(0.2, 0.8))
#'
#'
#' # When using ci as error, increase the 
#' # y.axis.adjust value to see the plot correctly!
#' threeFACTORplot(res,
#'    arrangement = c(2, 3, 1),
#'    bar.width = 0.8,
#'    fill = "Greens",
#'    xlab = "Drought",
#'    ylab = "Relative Expression",
#'    errorbar = "ci",
#'    y.axis.adjust = 8,
#'    y.axis.by = 2,
#'    letter.position.adjust = 0.6,
#'    legend.title = "Genotype",
#'    fontsize = 12,
#'    legend.position = c(0.2, 0.8),
#'    show.letters = TRUE)
#'
#'
#'


threeFACTORplot <- function(x,
                         arrangement = c(1, 2, 3),
                         bar.width = 0.5,
                         fill = "Reds",
                         xlab = "none",
                         ylab = "Relative Expression",
                         errorbar = "std",
                         y.axis.adjust = 0.5,
                         y.axis.by = 2,
                         letter.position.adjust = 0.3,
                         legend.title = "Legend Title",
                         legend.position = c(0.4, 0.8),
                         fontsize = 12,
                         fontsizePvalue = 7,
                         show.letters = TRUE,
                         axis.text.x.angle = 0,
                         axis.text.x.hjust = 0.5){
  
  x <- x$Result
  x <- x[, c(arrangement, 4:ncol(x))]
  RE <- x$RE
  std <- x$std
  LCL <- x$LCL
  UCL <- x$UCL
  
  
  pp1 <- ggplot(x, aes(x[,1], y = RE, fill = x[,2])) +
    geom_bar(stat = "identity", position = "dodge", width =  bar.width, col = "black") +
    geom_errorbar(aes(ymax = RE + std, ymin = RE),
                  position = position_dodge(bar.width), width = 0.15, color = "black") +
    facet_grid( ~ x[,3])+
    ylab(ylab) +
    theme_bw() +
    theme(legend.position = legend.position) +
    scale_fill_brewer(palette = fill) +
    scale_y_continuous(breaks = seq(0, max(x$RE) + max(x$std) +
                                      y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(x$RE) + max(x$std) + y.axis.adjust),
                       expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize)) +
    theme(legend.text = element_text(colour = "black", size = fontsize),
          legend.background = element_rect(fill = "transparent")) +
    guides(fill = guide_legend(title = legend.title, theme = theme(legend.title = element_text(size = 12, colour = "black")))) +
    theme(strip.background = element_rect(fill = "#F0FFFF")) +
    theme(strip.text = element_text(size = fontsize))
  
  if (show.letters) {
    pp1 <- pp1 + 
      geom_text(data = x, aes(label=letters, y = RE + std + letter.position.adjust), color = "black",
                show.legend = FALSE, position = position_dodge(bar.width), size = fontsizePvalue)
  }
  
  if(xlab == "none"){
    pp1 <- pp1 + 
      labs(x = NULL)
  }else{
    pp1 <- pp1 +
      xlab(xlab)
  }
  
  
  
  pp2 <- ggplot(x, aes(x[,1], y = RE, fill = x[,2])) +
    geom_bar(stat = "identity", position = "dodge", width =  bar.width, col = "black") +
    geom_errorbar(aes(ymax = RE + LCL, ymin = RE - UCL),
                  position = position_dodge(bar.width), width = 0.15, color = "black") +
    facet_grid( ~ x[,3])+
    ylab(ylab) +
    theme_bw() +
    theme(legend.position = legend.position) +
    scale_fill_brewer(palette = fill) +
    scale_y_continuous(breaks = seq(0, max(x$RE) + max(x$LCL) +
                                      y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(x$RE) + max(x$LCL) + y.axis.adjust),
                       expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize)) +
    theme(legend.text = element_text(colour = "black", size = fontsize),
          legend.background = element_rect(fill = "transparent")) +
    guides(fill = guide_legend(title = legend.title, theme = theme(legend.title = element_text(size = fontsize, colour = "black")))) +
    theme(strip.background = element_rect(fill = "#F0FFFF")) +
    theme(strip.text = element_text(size = fontsize))
  
  if (show.letters) {
    pp2 <- pp2 +
      geom_text(data = x, aes(label=letters, y = RE + LCL + letter.position.adjust), color = "black",
                show.legend = FALSE, position = position_dodge(bar.width), size = fontsizePvalue)
  }
  
  if(xlab == "none"){
    pp2 <- pp2 + 
      labs(x = NULL)
  }else{
    pp2 <- pp2 +
      xlab(xlab)
  }
  
  
  if(errorbar == "std") {
    out <- list(plot = pp1)
    
  } else if(errorbar == "ci") {
    out <- list(plot = pp2)
  }
  

  
  return(out)
}
