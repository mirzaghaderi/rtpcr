#' @title Bar plot of the average fold change of one target gene with two or more levels.
#' @description Bar plot of the relative expression of a gene along with the 95\% confidence interval and significance. There is another function, \code{qpcrTTESTplot}, that is used when the factor has two levels and represents Fold Change. \code{oneFACTORplot} is mainly useful for a one-factor experiment with more than two levels.
#' @details The \code{oneFACTORfcplot} function applies ANOVA based analysis where one target and one reference gene, that have been evaluated under two or more than two levels of a factor. It returns the bar plot of the average fold change for the target gene along with the 95\% CI and significance.
#' @author Ghader Mirzaghaderi
#' @export oneFACTORfcplot
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x a data frame. The data frame consists of 4 columns belonging to condition levels, E (efficiency), genes and Ct values, respectively. Each Ct in the following data frame is the mean of technical replicates. Complete amplification efficiencies of 2 is assumed here for all wells but the calculated efficienies can be used we well. We use this data set for Fold Change expression analysis of the target genes in treatment condition compared to normal condition.
#' @param numberOfrefGenes number of reference genes. Up to two reference genes can be handled.
#' @param levels a numeric vector with the length equal to the factor levels. First number indicates Control.
#' @param level.names  a vector determining level names in the x axis on the plot.
#' @param showCheckLevel a logical variable determining if the check level column be show in the plot or not.
#' @param width a positive number determining bar width.
#' @param fill  specify the fill color for the columns of the bar plot.
#' @param y.axis.adjust  a negative or positive value for reducing or increasing the length of the y axis.
#' @param letter.position.adjust adjust the distance between the signs and the error bars.
#' @param y.axis.by determines y axis step length
#' @param xlab  the title of the x axis
#' @param ylab  the title of the y axis
#' @param fontsize all fonts size of the plot
#' @return A list with 2 elements:
#' \describe{
#'   \item{Table}{Bar plot of the average fold change for target genes.}
#'   \item{plot}{Table of FC values and significance and the 95 percent CI as error bars.}
#' }
#' @examples
#'
#' # See sample data
#' data_1factor
#'
#'
#' # Producing the plot
#' oneFACTORfcplot(data_1factor,
#'                 numberOfrefGenes = 1,
#'                 levels = c(3, 2, 1),
#'                 level.names = "none",
#'                 showCheckLevel = TRUE,
#'                 width = 0.5,
#'                 fill = "skyblue",
#'                 y.axis.adjust = 1,
#'                 y.axis.by = 1,
#'                 letter.position.adjust = 0.3,
#'                 ylab = "Average Fold Change",
#'                 xlab = "Pairs",
#'                 fontsize = 12)
#'
#'


oneFACTORfcplot <- function(
    x,
    numberOfrefGenes,
    levels,
    showCheckLevel = TRUE,
    level.names = "none",
    width = 0.5,
    fill = "skyblue",
    y.axis.adjust = 1,
    y.axis.by = 1,
    letter.position.adjust = 0.1,
    ylab = "Average Fold Change",
    xlab = "Pairs",
    fontsize = 12){
  
  
  xfl <- x[,1]
  levels <- rev(levels)
  colnames(x)[1] <- "condition"
  x$condition <- levels[as.factor(xfl)]
  
  if(numberOfrefGenes == 1) {
    FINALDATA  <- qpcrANOVA(x, numberOfrefGenes = 1)$Final_data
    POSTHUC <- qpcrANOVA(x, numberOfrefGenes = 1)$Post_hoc_Test
  } else {
    FINALDATA  <- qpcrANOVA(x, numberOfrefGenes = 2)$Final_data
    POSTHUC <- qpcrANOVA(x, numberOfrefGenes = 2)$Post_hoc_Test
  }
  
  
  Nrows <- length(unique(FINALDATA[,1])[-1])
  withControl  <- POSTHUC[1:Nrows,]
  withControl
  
  # if check level be shown in the plot or not
  if(showCheckLevel == TRUE) {
    tableC <- rbind(data.frame(row.names = "1 - 1", FC = 1, pvalue=1, signif.=" ", LCL=0, UCL=0), withControl)
    # default level names of add level.names
    if (any(level.names == "none")) {
      rownames(tableC) <- rownames(tableC)
    } else {
      rownames(tableC) <- level.names
    }
  } else {
    tableC <- withControl
    # default level names of add level.names
    if (any(level.names == "none")) {
      rownames(tableC) <- rownames(tableC)
    } else {
      rownames(tableC) <- level.names[-1]
    }
  }
  
  
  pairs <- rownames(tableC)
  UCLp <- tableC$UCL
  LCLp <- tableC$LCL
  FCp <- tableC$FC
  significance <- tableC$signif.
  
  
  pfc <- ggplot(tableC, aes(factor(pairs, levels = rownames(tableC)), FCp)) +
    geom_col(col = "black", fill = fill, width = width) +
    geom_errorbar(aes(pairs, ymin = FCp - LCLp, ymax =  FCp + UCLp),
                  width=0.1) +
    geom_text(aes(label = significance,
                  x = pairs,
                  y = FCp + UCLp + letter.position.adjust),
              vjust = -0.5, size = 8) +
    ylab(ylab) + xlab(xlab) +
    theme_bw()+
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize)) +
    scale_y_continuous(breaks=seq(0, max(UCLp) + max(FCp) + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(UCLp) + max(FCp) + y.axis.adjust + y.axis.adjust), expand = c(0, 0)) +
    theme(legend.text = element_text(colour = "black", size = fontsize),
          legend.background = element_rect(fill = "transparent"))
  
  
  outlist <- list(plot = pfc,
                  Table = tableC)
  return(outlist)
}
