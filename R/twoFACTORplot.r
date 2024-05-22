#' @title Bar plot of the relative gene expression (\eqn{\Delta C_T} method) from the \code{qpcrANOVARE} output of a two-factorial experiment data
#' 
#' @description Bar plot of the relative expression (\eqn{\Delta C_T} method) of a gene along with the standard error (se), 95\% confidence interval (ci) and significance
#' 
#' @details The \code{twoFACTORplot} function generates the bar plot of the average fold change for target genes along with the significance, standard error (se) and the 95\% confidence interval (ci) as error bars.
#' 
#' @author Ghader Mirzaghaderi
#' @export twoFACTORplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @param res the FC data frame created by \code{qpcrANOVARE(x)$Result} function on a two factor data such as \code{data_2factor}.
#' @param x.axis.factor x-axis factor.
#' @param group.factor grouping factor.
#' @param width a positive number determining bar width.
#' @param fill  specify the fill color vector for the columns of the bar plot. One of the palettes in \code{\link[RColorBrewer]{display.brewer.all}} (e.g. "Reds" or "Blues", ...) can be applied.
#' @param y.axis.adjust  a negative or positive number for reducing or increasing the length of the y axis.
#' @param y.axis.by determines y axis step length.
#' @param show.errorbars show errorbars
#' @param errorbar Type of error bar, can be \code{se} or \code{ci}.
#' @param letter.position.adjust adjust the distance between the grouping letters to the error bars.
#' @param xlab  the title of the x axis.
#' @param ylab  the title of the y axis.
#' @param legend.position a two digit vector specifying the legend position.
#' @param fontsize size of all fonts  of the plot.
#' @param fontsizePvalue font size of the pvalue labels
#' @param show.letters a logical variable. If TRUE, mean grouping letters are added to the bars. 
#' @param axis.text.x.angle angle of x axis text
#' @param axis.text.x.hjust horizontal justification of x axis text
#' @return Bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @examples
#' 
#' # See a sample data frame
#' data_2factor
#'
#' # Before generating plot, the result table needs to be extracted as below:
#' res <- qpcrANOVARE(data_2factor, numberOfrefGenes = 1, block = NULL)$Result
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
#'    letter.position.adjust = 0.2,
#'    legend.position = c(0.2, 0.8),
#'    errorbar = "se")
#'    
#'
#' # Plotting the same data with 'Drought' as grouping factor
#' twoFACTORplot(res,
#'               x.axis.factor = Genotype,
#'               group.factor = Drought,
#'               xlab = "Genotype",
#'               fill = "Blues",
#'               fontsizePvalue = 5,
#'               errorbar = "ci")
#'               
#'               
#'               
#' # Combining FC results of two different genes:
#' a <- qpcrREPEATED(data_repeated_measure_1,
#'                   numberOfrefGenes = 1,
#'                   factor = "time", block = NULL)
#' 
#' b <- qpcrREPEATED(data_repeated_measure_2,
#'                   factor = "time",
#'                   numberOfrefGenes = 1, block = NULL)
#' 
#' a1 <- a$FC_statistics_of_the_main_factor
#' b1 <- b$FC_statistics_of_the_main_factor
#' c <- rbind(a1, b1)
#' c$gene <- factor(c(1,1,1,2,2,2))
#' c
#' 
#' twoFACTORplot(c, x.axis.factor = contrast, 
#'               group.factor = gene, fill = 'Reds',
#'               ylab = "FC", axis.text.x.angle = 45,
#'               errorbar = "se",
#'               axis.text.x.hjust = 1)
#'
#'




twoFACTORplot <- function(res, x.axis.factor, group.factor, width = 0.5, fill = "Blues", 
                          y.axis.adjust = 0.5, y.axis.by = 2, show.errorbars = TRUE,
                          errorbar, show.letters = TRUE, letter.position.adjust = 0.1,
                          ylab = "Relative Expression", xlab = "none", legend.position = c(0.09, 0.8),
                          fontsize = 12, fontsizePvalue = 5, axis.text.x.angle = 0, axis.text.x.hjust = 0.5){
  
  
  if (!("data.frame" %in% class(res))) {
    stop(deparse(substitute(res)), " is not a data frame")
  }
  
  if (missing(x.axis.factor)) {
    stop("argument 'x.axis.factor' is missing, with no default")
  }
  if (missing(group.factor)) {
    stop("argument 'group.factor' is missing, with no default")
  }
  if (missing(errorbar)) {
    stop("argument 'errorbar' is missing, with no default")
  }  
  
  
  
  
  b <- res
  se <- b$se
  LCL <- b$LCL
  UCL <- b$UCL
  
  
  if (any(grepl("letters", names(b)))) {
    b$letters <- gsub(" ", "", b$letters)
  }
  
  
  if (any(grepl("RE", names(b)))) {
    RE <- b$RE
    if(errorbar == "se") { 
  qq1 <- ggplot(b, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                       y = RE, group = {{ group.factor }}, fill = {{ group.factor }})) +
    geom_col(color = "black", position = "dodge", width = width) +
    scale_fill_brewer(palette = fill) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE)  +
    ylab(ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize),
          legend.background = element_rect(fill = "transparent")) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(se) + y.axis.adjust, by = y.axis.by), 
                       limits = c(0, max(RE) + max(se) + y.axis.adjust), expand = c(0, 0)) +
    theme(legend.position  = legend.position) +
    theme(legend.title = element_text(size = fontsize, color = "black")) 
  
  
  if (show.errorbars) {
    qq1 <-qq1 +
    geom_errorbar(aes(ymin = RE, ymax = RE + se), width = 0.2, position = position_dodge(width = 0.5))
  }
  
  
  if (show.letters) {
    qq1 <-qq1 + 
      geom_text(data = b, aes(label = letters, x = {{ x.axis.factor }}, y = RE + se + letter.position.adjust), 
                vjust = -0.5, size = fontsizePvalue, position = position_dodge(width = 0.5))
  }
  
  if(xlab == "none"){
    qq1 <- qq1 + 
      labs(x = NULL)
  }else{
    qq1 <- qq1 +
      xlab(xlab)
  }
    } else if(errorbar == "ci") {
  qq1 <- ggplot(b, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                       y = RE, group = {{ group.factor }}, fill = {{ group.factor }})) +
    geom_col(color = "black", position = "dodge", width = width) +
    scale_fill_brewer(palette = fill) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
    ylab(ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize),
          legend.background = element_rect(fill = "transparent")) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(LCL) + y.axis.adjust, by = y.axis.by), 
                       limits = c(0, max(RE) + max(LCL) + y.axis.adjust), expand = c(0, 0)) +
    theme(legend.position  = legend.position) +
    theme(legend.title = element_text(size = fontsize, color = "black")) 
  
  
  if (show.errorbars) {
    qq1 <-qq1 +
      geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, position = position_dodge(width = 0.5))
  }
  

  if (show.letters) {
    qq1 <- qq1 + 
      geom_text(data = b, aes(label = letters, x = {{ x.axis.factor }}, y = UCL + letter.position.adjust), 
                vjust = -0.5, size = fontsizePvalue, position = position_dodge(width = 0.5))
  }
  
  if(xlab == "none"){
    qq1 <- qq1 + 
      labs(x = NULL)
  }else{
    qq1 <- qq1 +
      xlab(xlab)
  }
    }
  
  # if(errorbar == "se") {
  #   out2 <- list(plot = qq1)
  #   
  # } else if(errorbar == "ci") {
  #   out2 <- list(plot = qq1)
  # }
  }
  
  
  
  
  
  
  
  if (any(grepl("FC", names(b)))) {
    letters <- b$sig
    FC <- b$FC
    if(errorbar == "se") {
    qq1 <- ggplot(b, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                         y = FC, group = {{ group.factor }}, fill = {{ group.factor }})) +
      geom_col(color = "black", position = "dodge", width = width) +
      scale_fill_brewer(palette = fill) +
      #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE)  +
      ylab(ylab) +
      theme_bw() +
      theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
            axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
            axis.title  = element_text(size = fontsize),
            legend.text = element_text(size = fontsize)) +
      scale_y_continuous(breaks = seq(0, max(FC) + max(se) + y.axis.adjust, by = y.axis.by), 
                         limits = c(0, max(FC) + max(se) + y.axis.adjust), expand = c(0, 0)) +
      theme(legend.position  = legend.position) +
      theme(legend.title = element_text(size = fontsize, color = "black"),
            legend.background = element_rect(fill = "transparent")) 
    
    
    if (show.errorbars) {
      qq1 <-qq1 +
        geom_errorbar(aes(ymin = FC, ymax = FC + se), width = 0.2, position = position_dodge(width = 0.5))
    }
    
    
    if (show.letters) {
      qq1 <-qq1 + 
        geom_text(data = b, aes(label = letters, x = {{ x.axis.factor }}, y = FC + se + letter.position.adjust), 
                  vjust = -0.5, size = fontsizePvalue, position = position_dodge(width = 0.5))
    }
    
    if(xlab == "none"){
      qq1 <- qq1 + 
        labs(x = NULL)
    }else{
      qq1 <- qq1 +
        xlab(xlab)
    }
    } else if(errorbar == "ci") {
    qq1 <- ggplot(b, aes(factor({{x.axis.factor}}, levels = as.character(unique({{x.axis.factor}}))),
                         y = FC, group = {{ group.factor }}, fill = {{ group.factor }})) +
      geom_col(color = "black", position = "dodge", width = width) +
      scale_fill_brewer(palette = fill) +
      #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
      ylab(ylab) +
      theme_bw() +
      theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
            axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
            axis.title  = element_text(size = fontsize),
            legend.text = element_text(size = fontsize),
            legend.background = element_rect(fill = "transparent")) +
      scale_y_continuous(breaks = seq(0, max(FC) + max(LCL) + y.axis.adjust, by = y.axis.by), 
                         limits = c(0, max(FC) + max(LCL) + y.axis.adjust), expand = c(0, 0)) +
      theme(legend.position  = legend.position) +
      theme(legend.title = element_text(size = fontsize, color = "black")) 
    
    
    if (show.errorbars) {
      qq1 <-qq1 +
        geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, position = position_dodge(width = 0.5))
    }
    
    
    if (show.letters) {
      qq1 <- qq1 + 
        geom_text(data = b, aes(label = letters, x = {{ x.axis.factor }}, y = UCL + letter.position.adjust), 
                  vjust = -0.5, size = fontsizePvalue, position = position_dodge(width = 0.5))
    }
    
    if(xlab == "none"){
      qq1 <- qq1 + 
        labs(x = NULL)
    }else{
      qq1 <- qq1 +
        xlab(xlab)
    }
    }  
    

  }
  
  
  return(qq1)
}



