#' @title Bar plot of the relative gene expression from a two-factor experiment
#' @description Bar plot of the relative expression of a gene along with the 95\% confidence interval and significance
#' @details The \code{twoFACTORplot} function generates the bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @author Ghader Mirzaghaderi
#' @export twoFACTORplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x a data frame created by \link{qpcrANOVA} function via running \code{res <- qpcrANOVA(x)$Result}.
#' @param x.axis.factor x-axis factor.
#' @param group.factor grouping factor.
#' @param width a positive number determining bar width.
#' @param fill  specify the fill color vector for the columns of the bar plot. One of the palettes in \code{\link[RColorBrewer]{display.brewer.all}} (e.g. "Reds" or "Blues", ...) can be applied.
#' @param y.axis.adjust  a negative or positive number for reducing or increasing the length of the y axis.
#' @param y.axis.by determines y axis step length.
#' @param errorbar Type of error bar, can be \code{std} or \code{ci}.
#' @param letter.position.adjust adjust the distance between the grouping letters to the error bars.
#' @param xlab  the title of the x axis.
#' @param ylab  the title of the y axis.
#' @param legend.position a two digit vector specifying the legend position.
#' @param fontsize size of all fonts  of the plot.
#' @param show.letters a logical variable. If TRUE, mean grouping letters are added to the bars. 
#' @return Bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @examples
#' 
#' # See a sample data frame
#' data_2factor
#'
#' # Before generating plot, the result table needs to be extracted as below:
#' res <- qpcrANOVA(data_2factor, numberOfrefGenes = 1)$Result
#' res
#'
#' # Plot of the 'res' data with 'Genotype' as grouping factor
#' twoFACTORplot(res,
#'    x.axis.factor = Drought,
#'    group.factor = Genotype,
#'    width = 0.5,
#'    fill = "Greens",
#'    y.axis.adjust = 1,
#'    y.axis.by = 2,
#'    ylab = "Relative Expression",
#'    xlab = "Drought Levels",
#'    show.letters = TRUE, 
#'    letter.position.adjust = 0.3,
#'    legend.position = c(0.09, 0.8),
#'    errorbar = "ci")
#'
#' # Plotting the same data with 'Drought' as grouping factor
#' twoFACTORplot(res,
#'    x.axis.factor = Genotype,
#'    group.factor = Drought,
#'    xlab = "Genotype",
#'    fill = "Blues",
#'    fontsize = 12,
#'    show.letters = FALSE)
#'
#'




twoFACTORplot <- function(x,
                          x.axis.factor,
                          group.factor,
                          width = 0.5,
                          fill = "Blues",
                          y.axis.adjust = 0.5,
                          y.axis.by = 2,
                          errorbar = "std",
                          show.letters = TRUE,
                          letter.position.adjust = 0.1,
                          ylab = "Relative Expression",
                          xlab = "Gene",
                          legend.position = c(0.09, 0.8),
                          fontsize = 12){
  RE <- x$RE
  std <- x$std
  LCL <- x$LCL
  UCL <- x$UCL
  
  qq1 <- ggplot(x, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                       y = RE, group = {{ group.factor }}, fill = {{ group.factor }})) +
    geom_col(color = "black", position = "dodge", width = width) +
    scale_fill_brewer(palette = fill) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = RE, ymax = RE + std), width = 0.2, position = position_dodge(width = 0.5)) +
    ylab(ylab) +
    xlab(xlab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize)) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(std) + y.axis.adjust, by = y.axis.by), 
                       limits = c(0, max(RE) + max(std) + y.axis.adjust), expand = c(0, 0)) +
    theme(legend.position  = legend.position) +
    theme(legend.title = element_text(size = fontsize, color = "black"))
  
  if (show.letters) {
    qq1 <-qq1 + 
      geom_text(data = x, aes(label = letters, x = {{ x.axis.factor }}, y = RE + std + letter.position.adjust), 
                vjust = -0.5, size = 4, position = position_dodge(width = 0.5))
  }
  
  
  
  
  qq2 <- ggplot(x, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                       y = RE, group = {{ group.factor }}, fill = {{ group.factor }})) +
    geom_col(color = "black", position = "dodge", width = width) +
    scale_fill_brewer(palette = fill) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = RE - UCL, ymax = RE + LCL), width = 0.2, position = position_dodge(width = 0.5)) +
    ylab(ylab) +
    xlab(xlab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize)) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(LCL) + y.axis.adjust, by = y.axis.by), 
                       limits = c(0, max(RE) + max(LCL) + y.axis.adjust), expand = c(0, 0)) +
    theme(legend.position  = legend.position) +
    theme(legend.title = element_text(size = fontsize, color = "black"))
  
  if (show.letters) {
    qq2 <- qq2 + 
      geom_text(data = x, aes(label = letters, x = {{ x.axis.factor }}, y = RE + LCL + letter.position.adjust), 
                vjust = -0.5, size = 4, position = position_dodge(width = 0.5))
  }
  
  
  if(errorbar == "std") {
    print(qq1)
  } else if(errorbar == "ci") {
    print(qq2)
  }
}
