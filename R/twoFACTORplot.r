#' @title Bar plot of the relative gene expression from a two-factor experiment
#' @description Bar plot of the relative expression (RE) of a gene along with the 95\% confidence interval and significance
#' @details The \code{twoFACTORplot} function generates the bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @author Ghader Mirzaghaderi
#' @export twoFACTORplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x an object created by \link{qpcrANOVA} function via running \code{res <- qpcrANOVA(x)}.
#' @param x.axis.factor x-axis factor.
#' @param group.factor grouping factor.
#' @param width a positive number determining bar width.
#' @param fill  specify the fill color vector for the columns of the bar plot. One of the palettes in \code{\link[RColorBrewer]{display.brewer.all}} (e.g. "Reds" or "Blues", ...) can be applied.
#' @param y.axis.adjust  a negative or positive number for reducing or increasing the length of the y axis.
#' @param y.axis.by determines y axis step length.
#' @param show.errorbars show errorbars
#' @param show.points show points
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
#' res <- qpcrANOVA(data_2factor, numberOfrefGenes = 1)
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
#'    show.letters = FALSE,
#'    show.points = TRUE,
#'    show.errorbars = FALSE)
#'
#'




twoFACTORplot <- function(x,
                          x.axis.factor,
                          group.factor,
                          width = 0.5,
                          fill = "Blues",
                          y.axis.adjust = 0.5,
                          y.axis.by = 2,
                          show.errorbars = TRUE,
                          errorbar = "std",
                          show.letters = TRUE,
                          show.points = FALSE,
                          letter.position.adjust = 0.1,
                          ylab = "Relative Expression",
                          xlab = "Gene",
                          legend.position = c(0.09, 0.8),
                          fontsize = 12){
  b <- x$Result
  a <- x$Final_data
  RE <- b$RE
  std <- b$std
  LCL <- b$LCL
  UCL <- b$UCL
  wDCt <- a$wDCt
  
  qq1 <- ggplot(b, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                       y = RE, group = {{ group.factor }}, fill = {{ group.factor }})) +
    geom_col(color = "black", position = "dodge", width = width) +
    scale_fill_brewer(palette = fill) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE)  +
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
  
  
  if (show.errorbars) {
    qq1 <-qq1 +
    geom_errorbar(aes(ymin = RE, ymax = RE + std), width = 0.2, position = position_dodge(width = 0.5))
  }
  
  
  if (show.points) {
    qq1 <-qq1 + 
    geom_point(data = a, aes(x = factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))), 
                             y = (10^(-wDCt)), group = as.factor({{ group.factor }}), fill = as.factor({{ group.factor }})), 
               position = position_dodge(width = width),color = "grey30",  fill = "grey30", shape = 21, size = 2)
  }
  if (show.letters) {
    qq1 <-qq1 + 
      geom_text(data = b, aes(label = letters, x = {{ x.axis.factor }}, y = RE + std + letter.position.adjust), 
                vjust = -0.5, size = 4, position = position_dodge(width = 0.5))
  }
  
  
  
  
  qq2 <- ggplot(b, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                       y = RE, group = {{ group.factor }}, fill = {{ group.factor }})) +
    geom_col(color = "black", position = "dodge", width = width) +
    scale_fill_brewer(palette = fill) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
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
  
  
  if (show.errorbars) {
    qq2 <-qq2 +
      geom_errorbar(aes(ymin = RE - UCL, ymax = RE + LCL), width = 0.2, position = position_dodge(width = 0.5))
  }
  
  if (show.points) {
    qq2 <-qq2 + 
      geom_point(data = a, aes(x = factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))), 
                               y = (10^(-wDCt)), group = as.factor({{ group.factor }}), fill = as.factor({{ group.factor }})), 
                 position = position_dodge(width = width),color = "grey30",  fill = "grey30", shape = 21, size = 3)
  }
  if (show.letters) {
    qq2 <- qq2 + 
      geom_text(data = b, aes(label = letters, x = {{ x.axis.factor }}, y = RE + LCL + letter.position.adjust), 
                vjust = -0.5, size = 4, position = position_dodge(width = 0.5))
  }
  
  
  if(errorbar == "std") {
    out2 <- list(plot = qq1)
    
  } else if(errorbar == "ci") {
    out2 <- list(plot = qq2)
  }
  
  return(out2)
}


