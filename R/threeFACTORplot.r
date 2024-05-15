#' @title Bar plot of the relative gene expression (\eqn{\Delta C_T} method) from the \code{qpcrANOVARE} output of a a three-factorial experiment data
#' 
#' @description Bar plot of the relative expression (\eqn{\Delta C_T} method) of a gene along with the confidence interval and significance
#' 
#' @details The \code{threeFACTORplot} function generates the bar plot of the average fold change for target genes along with the significance, standard error (se) and the 95\% confidence interval (ci).
#' 
#' @author Ghader Mirzaghaderi
#' @export threeFACTORplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @param res the FC data frame created by \code{qpcrANOVARE(x)$Result} function on a three factorial data such as \code{data_3factor} example data frame.
#' @param arrangement order based on the columns in the output table (e.g. c(2,3,1) or c(1,3,2)) affecting factor arrangement of the output graph.
#' @param bar.width a positive number determining bar width.
#' @param fill  a color vector specifying the fill color for the columns of the bar plot. One of the palettes in \code{\link[RColorBrewer]{display.brewer.all}} (e.g. "Reds" or "Blues", ...) can be applied.
#' @param xlab  the title of the x axis
#' @param ylab  the title of the y axis
#' @param errorbar Type of error bar, can be \code{se} or \code{ci}.
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
#' res <- qpcrANOVARE(data_3factor, numberOfrefGenes = 1)$Result
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
#'    errorbar = "se",
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
#'    y.axis.adjust = 1,
#'    y.axis.by = 2,
#'    letter.position.adjust = 0.6,
#'    legend.title = "Genotype",
#'    fontsize = 12,
#'    legend.position = c(0.2, 0.8),
#'    show.letters = TRUE)
#'
#'
#'


threeFACTORplot <- function(res,
                         arrangement = c(1, 2, 3),
                         bar.width = 0.5,
                         fill = "Reds",
                         xlab = "none",
                         ylab = "Relative Expression",
                         errorbar = "se",
                         y.axis.adjust = 0.5,
                         y.axis.by = 2,
                         letter.position.adjust = 0.3,
                         legend.title = "Legend Title",
                         legend.position = c(0.4, 0.8),
                         fontsize = 12,
                         fontsizePvalue = 5,
                         show.letters = TRUE,
                         axis.text.x.angle = 0,
                         axis.text.x.hjust = 0.5){
  
  x <- res
  x <- x[, c(arrangement, 4:ncol(x))]
  se <- x$se
  LCL <- x$LCL
  UCL <- x$UCL
  
  
  if (any(grepl("letters", names(x)))) {
    x$letters <- gsub(" ", "", x$letters)
  }
  
  
  if (any(grepl("RE", names(x)))) {
    RE <- x$RE
    if(errorbar == "se") {
  pp1 <- ggplot(x, aes(x[,1], y = RE, fill = x[,2])) +
    geom_bar(stat = "identity", position = "dodge", width =  bar.width, col = "black") +
    geom_errorbar(aes(ymax = RE + se, ymin = RE),
                  position = position_dodge(bar.width), width = 0.15, color = "black") +
    facet_grid( ~ x[,3])+
    ylab(ylab) +
    theme_bw() +
    theme(legend.position = legend.position) +
    scale_fill_brewer(palette = fill) +
    scale_y_continuous(breaks = seq(0, max(x$RE) + max(x$se) +
                                      y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(x$RE) + max(x$se) + y.axis.adjust),
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
      geom_text(data = x, aes(label=letters, y = RE + se + letter.position.adjust), color = "black",
                show.legend = FALSE, position = position_dodge(bar.width), size = fontsizePvalue)
  }
  
  if(xlab == "none"){
    pp1 <- pp1 + 
      labs(x = NULL)
  }else{
    pp1 <- pp1 +
      xlab(xlab)
  }
  } else if(errorbar == "ci") {
  pp1 <- ggplot(x, aes(x[,1], y = RE, fill = x[,2])) +
    geom_bar(stat = "identity", position = "dodge", width =  bar.width, col = "black") +
    geom_errorbar(aes(ymin = LCL, ymax = UCL),
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
    pp1 <- pp1 +
      geom_text(data = x, aes(label=letters, y = UCL + letter.position.adjust), color = "black",
                show.legend = FALSE, position = position_dodge(bar.width), size = fontsizePvalue)
  }
  
  if(xlab == "none"){
    pp1 <- pp1 + 
      labs(x = NULL)
  }else{
    pp1 <- pp1 +
      xlab(xlab)
  }
  } 
  
  # if(errorbar == "se") {
  #   out <- list(plot = pp1)
  #   
  # } else if(errorbar == "ci") {
  #   out <- list(plot = pp1)
  # }
  }

  
  
  
  
  
  
  
  
  
  
  
  if (any(grepl("FC", names(x)))) {
    letters <- x$sig
    FC <- x$FC
    if(errorbar == "se") {
    pp1 <- ggplot(x, aes(x[,1], y = FC, fill = x[,2])) +
      geom_bar(stat = "identity", position = "dodge", width =  bar.width, col = "black") +
      geom_errorbar(aes(ymax = FC + se, ymin = FC),
                    position = position_dodge(bar.width), width = 0.15, color = "black") +
      facet_grid( ~ x[,3])+
      ylab(ylab) +
      theme_bw() +
      theme(legend.position = legend.position) +
      scale_fill_brewer(palette = fill) +
      scale_y_continuous(breaks = seq(0, max(x$FC) + max(x$se) +
                                        y.axis.adjust, by = y.axis.by),
                         limits = c(0, max(x$FC) + max(x$se) + y.axis.adjust),
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
        geom_text(data = x, aes(label=letters, y = FC + se + letter.position.adjust), color = "black",
                  show.legend = FALSE, position = position_dodge(bar.width), size = fontsizePvalue)
    }
    
    if(xlab == "none"){
      pp1 <- pp1 + 
        labs(x = NULL)
    }else{
      pp1 <- pp1 +
        xlab(xlab)
    }
    } else if(errorbar == "ci") {
    pp1 <- ggplot(x, aes(x[,1], y = FC, fill = x[,2])) +
      geom_bar(stat = "identity", position = "dodge", width =  bar.width, col = "black") +
      geom_errorbar(aes(ymin = LCL, ymax = UCL),
                    position = position_dodge(bar.width), width = 0.15, color = "black") +
      facet_grid( ~ x[,3])+
      ylab(ylab) +
      theme_bw() +
      theme(legend.position = legend.position) +
      scale_fill_brewer(palette = fill) +
      scale_y_continuous(breaks = seq(0, max(x$FC) + max(x$LCL) +
                                        y.axis.adjust, by = y.axis.by),
                         limits = c(0, max(x$FC) + max(x$LCL) + y.axis.adjust),
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
      pp1 <- pp1 +
        geom_text(data = x, aes(label=letters, y = UCL + letter.position.adjust), color = "black",
                  show.legend = FALSE, position = position_dodge(bar.width), size = fontsizePvalue)
    }
    
    if(xlab == "none"){
      pp1 <- pp1 + 
        labs(x = NULL)
    }else{
      pp1 <- pp1 +
        xlab(xlab)
    }
    } 
    
    # if(errorbar == "se") {
    #   out <- list(plot = pp1)
    #   
    # } else if(errorbar == "ci") {
    #   out <- list(plot = pp1)
    # }
  }
  
  
  
  return(pp1)
}
