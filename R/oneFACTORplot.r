#' @title Bar plot of the relative gene expression from a one-factor experiment
#' @description Bar plot of the relative expression of a gene along with the 95\% confidence interval and significance. There is another function, \code{qpcrTTESTplot}, that is used when the factor has two levels and represents Fold Change. \code{oneFACTORplot} is mainly useful for a one-factor experiment with more than two levels.
#' @details The \code{oneFACTORplot} function generates the bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @author Ghader Mirzaghaderi
#' @export oneFACTORplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x a data frame created by \link{qpcrANOVA} function via running \code{res <- qpcrANOVA(x)$Result}.
#' @param width a positive number determining bar width.
#' @param fill  specify a fill color.
#' @param y.axis.adjust  a negative or positive number for reducing or increasing the length of the y axis.
#' @param y.axis.by determines y axis step length.
#' @param errorbar Type of error bar, can be \code{std} or \code{ci}.
#' @param letter.position.adjust adjust the distance between the grouping letters to the error bars.
#' @param xlab  the title of the x axis.
#' @param ylab  the title of the y axis.
#' @param fontsize size of all fonts  of the plot.
#' @param fontsizePvalue font size of the pvalue labels
#' @param show.letters a logical variable. If TRUE, mean grouping letters are added to the bars.
#' @return Bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @examples
#'
#' # Before plotting, the result needs to be extracted as below:
#' out <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)$Result
#'
#' # Bar plot
#' oneFACTORplot(out,
#'          width = 0.2,
#'          fill = "skyblue",
#'          y.axis.adjust = 0,
#'          y.axis.by = 0.2,
#'          errorbar = "std",
#'          show.letters = TRUE,
#'          letter.position.adjust = 0.05,
#'          ylab = "Relative Expression",
#'          xlab = "Factor Levels",
#'          fontsize = 12)
#'
#'



oneFACTORplot <- function(x,
                          width = 0.2,
                          fill = "skyblue",
                          y.axis.adjust = 0.5,
                          y.axis.by = 2,
                          errorbar = "std",
                          show.letters = TRUE,
                          letter.position.adjust = 0.1,
                          ylab = "Relative Expression",
                          xlab = "Factor",
                          fontsize = 12,
                          fontsizePvalue = 7){

  RE <- x$RE
  std <- x$std
  LCL <- x$LCL
  UCL <- x$UCL

  q1f1 <- ggplot(x, aes(rownames(x), y = RE, group = rownames(x))) +
    geom_col(color = "black", fill = fill, width = width) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = RE, ymax = RE + std), width = 0.1) +
    ylab(ylab) +
    xlab(xlab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize)) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(std) + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(RE) + max(std) + y.axis.adjust), expand = c(0, 0))

  if (show.letters) {
    q1f1 <-q1f1 +
      geom_text(data = x, aes(label = letters, x = rownames(x), y = RE + std + letter.position.adjust),
                vjust = -0.5, size = fontsizePvalue)
  }




  q1f2 <- ggplot(x, aes(rownames(x), y = RE, group = rownames(x))) +
    geom_col(color = "black", fill = fill, width = width) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = RE - UCL, ymax = RE + LCL), width = 0.1) +
    ylab(ylab) +
    xlab(xlab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize)) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(LCL) + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(RE) + max(LCL) + y.axis.adjust), expand = c(0, 0))

  if (show.letters) {
    q1f2 <- q1f2 +
      geom_text(data = x, aes(label = letters, x = rownames(x), y = RE + LCL + letter.position.adjust),
                vjust = -0.5, size = fontsizePvalue)
  }


  
  if(errorbar == "std") {
    out1 <- list(plot = q1f1)
    
  } else if(errorbar == "ci") {
    out1 <- list(plot = q1f2)
  }
  
  return(out1)
}


