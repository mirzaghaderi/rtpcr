#' @title Bar plot of the relative gene expression (RE) from the \code{qpcrANOVA} output of a one-factor experiment data
#' @description Bar plot of the relative expression of a gene along with the standard error (se), 95\% confidence interval (ci) and significance. \code{oneFACTORplot} is mainly useful for a one-factor experiment with more than two levels.
#' @details The \code{oneFACTORplot} function generates the bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @author Ghader Mirzaghaderi
#' @export oneFACTORplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import agricolae
#' @param res an FC data frame object created by \code{qpcrANOVA(x)$Result} function on a one factor data such as \code{data_1factor}.
#' @param width a positive number determining bar width.
#' @param fill  specify a fill color.
#' @param y.axis.adjust  a negative or positive number for reducing or increasing the length of the y axis.
#' @param y.axis.by determines y axis step length.
#' @param errorbar Type of error bar, can be \code{se} or \code{ci}.
#' @param letter.position.adjust adjust the distance between the grouping letters to the error bars.
#' @param xlab  the title of the x axis.
#' @param ylab  the title of the y axis.
#' @param fontsize size of all fonts  of the plot.
#' @param fontsizePvalue font size of the pvalue labels
#' @param show.letters a logical variable. If TRUE, mean grouping letters are added to the bars.
#' @param axis.text.x.angle angle of x axis text
#' @param axis.text.x.hjust horizontal justification of x axis text
#' @return Bar plot of the average fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' @examples
#'
#' # Before plotting, the result needs to be extracted as below:
#' res <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)$Result
#'
#' # Bar plot
#' oneFACTORplot(res,
#'          width = 0.2,
#'          fill = "skyblue",
#'          y.axis.adjust = 0,
#'          y.axis.by = 0.2,
#'          errorbar = "se",
#'          show.letters = TRUE,
#'          letter.position.adjust = 0.05,
#'          ylab = "Relative Expression",
#'          xlab = "Factor Levels",
#'          fontsize = 12)
#'
#'



oneFACTORplot <- function(res,
                          width = 0.2,
                          fill = "skyblue",
                          y.axis.adjust = 0.5,
                          y.axis.by = 2,
                          errorbar = "se",
                          show.letters = TRUE,
                          letter.position.adjust = 0.1,
                          ylab = "Relative Expression",
                          xlab = "none",
                          fontsize = 12,
                          fontsizePvalue = 7,
                          axis.text.x.angle = 0,
                          axis.text.x.hjust = 0.5){

  
  x <- res
  LCL <- x$LCL
  UCL <- x$UCL
  se <- x$se
  
  if (any(grepl("RE", names(x)))) {
  RE <- x$RE
  
  
  q1f1 <- ggplot(x, aes(rownames(x), y = RE, group = rownames(x))) +
    geom_col(color = "black", fill = fill, width = width) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = RE, ymax = RE + se), width = 0.1) +
    ylab(ylab) +
    xlab(xlab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize)) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(se) + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(RE) + max(se) + y.axis.adjust), expand = c(0, 0))

  if (show.letters) {
    q1f1 <-q1f1 +
      geom_text(data = x, aes(label = letters, x = rownames(x), y = RE + se + letter.position.adjust),
                vjust = -0.5, size = fontsizePvalue)
  }
  
  if(xlab == "none"){
    q1f1 <- q1f1 + 
      labs(x = NULL)
  }else{
    q1f1 <- q1f1 +
      xlab(xlab)
  }



  q1f2 <- ggplot(x, aes(rownames(x), y = RE, group = rownames(x))) +
    geom_col(color = "black", fill = fill, width = width) +
    #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.1) +
    ylab(ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize),
          legend.text = element_text(size = fontsize)) +
    scale_y_continuous(breaks = seq(0, max(RE) + max(LCL) + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(RE) + max(LCL) + y.axis.adjust), expand = c(0, 0))

  if (show.letters) {
    q1f2 <- q1f2 +
      geom_text(data = x, aes(label = letters, x = rownames(x), y = UCL + letter.position.adjust),
                vjust = -0.5, size = fontsizePvalue)
  }

  if(xlab == "none"){
    q1f2 <- q1f2 + 
      labs(x = NULL)
  }else{
    q1f2 <- q1f2 +
      xlab(xlab)
  }
  

  
  if(errorbar == "se") {
    out1 <- list(plot = q1f1)
    
  } else if(errorbar == "ci") {
    out1 <- list(plot = q1f2)
  }
  }
  

  
  
  
  
  
  if (any(grepl("FC", names(x)))) {
    x$FC <- as.numeric(x$FC)
    letters <- x$sig
    FC <- x$FC
    
    
    q1f1 <- ggplot(x, aes(contrast, y = FC, group = contrast)) +
      geom_col(color = "black", fill = fill, width = width) +
      #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
      geom_errorbar(aes(ymin = FC, ymax = FC + se), width = 0.1) +
      ylab(ylab) +
      xlab(xlab) +
      theme_bw() +
      theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
            axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
            axis.title  = element_text(size = fontsize),
            legend.text = element_text(size = fontsize)) +
      scale_y_continuous(breaks = seq(0, max(FC) + max(se) + y.axis.adjust, by = y.axis.by),
                         limits = c(0, max(FC) + max(se) + y.axis.adjust), expand = c(0, 0))
    
    if (show.letters) {
      q1f1 <-q1f1 +
        geom_text(data = x, aes(label = letters, x = contrast, y = FC + se + letter.position.adjust),
                  vjust = -0.5, size = fontsizePvalue)
    }
    
    if(xlab == "none"){
      q1f1 <- q1f1 + 
        labs(x = NULL)
    }else{
      q1f1 <- q1f1 +
        xlab(xlab)
    }
    
    
    
    q1f2 <- ggplot(x, aes(contrast, y = FC, group = contrast)) +
      geom_col(color = "black", fill = fill, width = width) +
      #geom_hline(aes(yintercept = 1), col = "red", linetype = 2, show.legend = FALSE) +
      geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.1) +
      ylab(ylab) +
      theme_bw() +
      theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
            axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
            axis.title  = element_text(size = fontsize),
            legend.text = element_text(size = fontsize)) +
      scale_y_continuous(breaks = seq(0, max(FC) + max(LCL) + y.axis.adjust, by = y.axis.by),
                         limits = c(0, max(FC) + max(LCL) + y.axis.adjust), expand = c(0, 0))
    
    if (show.letters) {
      q1f2 <- q1f2 +
        geom_text(data = x, aes(label = letters, x = contrast, y = UCL + letter.position.adjust),
                  vjust = -0.5, size = fontsizePvalue)
    }
    
    if(xlab == "none"){
      q1f2 <- q1f2 + 
        labs(x = NULL)
    }else{
      q1f2 <- q1f2 +
        xlab(xlab)
    }
    
    
    
    if(errorbar == "se") {
      out1 <- list(plot = q1f1)
      
    } else if(errorbar == "ci") {
      out1 <- list(plot = q1f2)
    }
  }
  
  return(out1)
}


