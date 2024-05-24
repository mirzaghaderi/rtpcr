#' @title Bar plot of the average fold change (\eqn{\Delta \Delta C_T} method) of target genes
#' 
#' @description Bar plot of the fold change (\eqn{\Delta \Delta C_T} method) values for for any number of target 
#' genes under a two-level conditional experimental (e.g. control and treatment).
#' 
#' @details The \code{qpcrTTESTplot} function applies a t.test based analysis to any number of target genes 
#' along with one or two reference gene(s), that have been evaluated under control and treatment conditions. 
#' It returns the bar plot of the fold change (FC) values for target genes along with the 95\% CI and significance. 
#' Sampling may be unpaired or paired. Unpaired and paired samples are commonly analyzed using unpaired and 
#' paired t-test, respectively.Paired samples in quantitative PCR refer to two sample data that are collected 
#' from one set of individuals 
#' at two different conditions, for example before and after a treatment or at two different time points. While 
#' for unpaired samples, two sets of individuals are used: one under untreated and the other set under treated 
#' condition. Paired samples allow to compare gene expression changes within the same individual, reducing 
#' inter-individual variability. 
#' 
#' @author Ghader Mirzaghaderi
#' @export qpcrTTESTplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' 
#' @param x a data frame. The data frame consists of 4 columns belonging to condition levels, 
#' E (efficiency), genes and Ct values, respectively. Each Ct in the following data frame is 
#' the mean of technical replicates. Complete amplification efficiencies of 2 is assumed here 
#' for all wells but the calculated efficienies can be used we well. We use this data set for 
#' fold change expression analysis of the target genes in treatment condition compared to 
#' normal condition. See \href{../doc/vignette.html}{\code{vignette}} for details.
#' @param numberOfrefGenes number of reference genes. Up to two reference genes can be handled.
#' @param order a vector determining genes order on the output graph.
#' @param paired  a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param width a positive number determining bar width.
#' @param fill  specify the fill color for the columns of the bar plot.
#' @param y.axis.adjust  a negative or positive value for reducing or increasing the length of the y axis.
#' @param letter.position.adjust adjust the distance between the signs and the error bars.
#' @param y.axis.by determines y axis step length
#' @param xlab  the title of the x axis
#' @param ylab  the title of the y axis
#' @param fontsize fonts size of the plot
#' @param fontsizePvalue font size of the pvalue labels
#' @param axis.text.x.angle angle of x axis text
#' @param axis.text.x.hjust horizontal justification of x axis text
#' @param errorbar Type of error bar, can be \code{se} or \code{ci}.
#' @return Bar  plot of the average fold change for target genes along with the significance and the 95 percent CI as error bars.
#' @examples
#'
#' # See a sample data frame
#' data_ttest
#'
#'
#' qpcrTTESTplot(data_ttest, 
#'               numberOfrefGenes = 1,
#'               errorbar = "ci")
#'
#'
#' # Producing the plot
#' qpcrTTESTplot(data_ttest,
#'               numberOfrefGenes = 1,
#'               order = c("C2H2-01", "C2H2-12", "C2H2-26"),
#'               paired = FALSE,
#'               var.equal = TRUE,
#'               width = 0.5,
#'               fill = "skyblue",
#'               y.axis.adjust = 0,
#'               y.axis.by = 2,
#'               letter.position.adjust = 0.3,
#'               ylab = "Fold Change in Treatment vs Control",
#'               xlab = "Gene",
#'               errorbar = "se")
#'
#'


qpcrTTESTplot <- function(x, order = "none", numberOfrefGenes, paired = FALSE,
                          var.equal = TRUE, width = 0.5, fill = "skyblue", y.axis.adjust = 0,
                          y.axis.by = 2, letter.position.adjust = 0.3, ylab = "Average Fold Change",
                          xlab = "none", fontsize = 12, fontsizePvalue = 7, axis.text.x.angle = 0,
                          axis.text.x.hjust = 0.5, errorbar = "se"){

  default.order <- unique(x[,2])[-length(unique(x[,2]))]

  # convert_to_character function
  convert_to_character <- function(numbers) {
    characters <- character(length(numbers))  # Initialize a character vector to store the results

    for (i in seq_along(numbers)) {
      if (numbers[i] < 0.01) {
        characters[i] <- "**"
      } else if (numbers[i] < 0.05) {
        characters[i] <- "*"
      } else if (numbers[i] < 0.1) {
        characters[i] <- "."
      } else {
        characters[i] <- "ns"
      }
    }

    return(characters)
  }


  
  

  
  if(numberOfrefGenes == 1) {
    TTESTRES <- qpcrTTEST(x, numberOfrefGenes = 1, paired = paired, var.equal = var.equal)
  } else {
    TTESTRES <- qpcrTTEST(x, numberOfrefGenes = 2, paired = paired, var.equal = var.equal)
  }
  
  
  
  
  
  
  # Getting barplot:
  df2 <- as.data.frame(TTESTRES[[2]])

  # Convert the column to a factor with specified levels
  if(any(order == "none")){
    df2$Gene <- factor(df2$Gene, levels = default.order)
  } else {
    df2$Gene <- factor(df2$Gene, levels = order)
  }
 

  # Order the data frame based on the specified column
  df2 <- df2[order(df2$Gene), ]

  Gene <- df2$Gene
  Gene <- factor(Gene, levels = Gene)
  Fold_Change <- df2$FC
  Lower.Er <- df2$LCL
  Upper.Er <- df2$UCL
  pvalue <- df2$pvalue
  se <- df2$se

  p <- ggplot(df2, aes(Gene, as.numeric(Fold_Change))) +
    geom_col(col = "black", fill = fill, width = width)
    
    
    if(errorbar == "ci") {
      p <- p +
        geom_errorbar(aes(Gene, ymin=as.numeric(Lower.Er), ymax = as.numeric(Upper.Er)), width=0.1)
    } else if(errorbar == "se") {
      p <- p +
        geom_errorbar(aes(Gene, ymin = 2^(log2(as.numeric(Fold_Change)) - as.numeric(se)), ymax =  2^(log2(as.numeric(Fold_Change)) + as.numeric(se))), width=0.1)
    }
    

  
  if(errorbar == "ci") {
    p <- p +
      geom_text(aes(label = convert_to_character(pvalue),
                    x = Gene,
                    y = as.numeric(Upper.Er)),
                vjust = -0.5, size = fontsizePvalue)
  } else if(errorbar == "se") {
    p <- p +
      geom_text(aes(label = convert_to_character(pvalue),
                    x = Gene,
                    y = 2^(log2(as.numeric(Fold_Change)) + as.numeric(se))),
                vjust = -0.5, size = fontsizePvalue)
  }

    
    
    p <- p + ylab(ylab) + 
    theme_bw()+
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize)) 
    
    
    if(errorbar == "ci") {
      p <- p +
        scale_y_continuous(breaks=seq(0, max(as.numeric(Fold_Change)) + max(as.numeric(Upper.Er) + y.axis.adjust) + 1, by = y.axis.by),
                           limits = c(0, max(as.numeric(Fold_Change)) + max(as.numeric(Upper.Er) + y.axis.adjust) + 1), expand = c(0, 0))
    } else if(errorbar == "se") {
      p <- p +
        scale_y_continuous(breaks = seq(0, max(as.numeric(Fold_Change)) + max(as.numeric(se) + y.axis.adjust) + 1, by = y.axis.by),
                           limits = c(0, max(as.numeric(Fold_Change)) + max(as.numeric(se) + y.axis.adjust) + 1), expand = c(0, 0))
      }
    
    
    
    
      p <- p + theme(legend.text = element_text(colour = "black", size = fontsize),
          legend.background = element_rect(fill = "transparent"))
  
  if(xlab == "none"){
    p <- p + 
      labs(x = NULL)
  }else{
    p <- p +
      xlab(xlab)
  }

  out <- list(plot = p)
  return(out)
}

