#' @title ANOVA based on a factorial design or analysis of covariance plus bar plot of FC
#' @description Analysis of variance of a factorial design or analysis of covariance can be performed using `qpcrANCOVA` function, if more than a factor exist. Bar plot of fold changes (FCs) along with the 95\% confidence interval is also returned by the `qpcrANCOVA` function. There is also a function called \code{oneFACTORplot} which returns FC values and related plot for a one-factor-experiment with more than two levels. 
#' @details The \code{qpcrANCOVA} function applies ANOVA based analysis where one target and one reference gene, that have been evaluated under two or more than two levels of a factor. It returns the bar plot of the average fold change for the target gene along with the 95\% CI and significance.
#' @author Ghader Mirzaghaderi
#' @export qpcrANCOVA
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x a data frame. The data frame consists of 4 columns belonging to condition levels, E (efficiency), genes and Ct values, respectively. Each Ct in the following data frame is the mean of technical replicates. Complete amplification efficiencies of 2 is assumed here for all wells but the calculated efficienies can be used we well. We use this data set for Fold Change expression analysis of the target genes in treatment condition compared to normal condition.
#' @param numberOfrefGenes number of reference genes. Up to two reference genes can be handled.
#' @param factorial a logical variable showing if the data is analysed based on factorial arrangement. 
#' @param factorDesired the factor for which the FC is calculated for its level pairs.
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
#' @param block column name of the blocking factor (for correct column arrangement see example data.)
#' @param p.adj Method for adjusting p values (see p.adjust)
#' @return A list with 2 elements:
#' \describe{
#'   \item{plot}{Bar plot of the average fold change for desired factor levels.}
#'   \item{Table}{Table of FC values and significance and the 95 percent CI as error bars.}
#' }
#' @examples
#'
#' # See sample data
#' data_2factor
#'
#'
#' qpcrANCOVA(data_2factor, 
#'            numberOfrefGenes = 1, 
#'            factorial = TRUE, 
#'            factorDesired = 2,  
#'            levels = c(3, 2, 1), 
#'            showCheckLevel = TRUE)
#'            
#'            
#' 

qpcrANCOVA <- function(x,
                      numberOfrefGenes,
                      block = NULL,
                      factorial = TRUE,
                      factorDesired = 1,
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
                      fontsize = 12,
                      p.adj = c("none","holm","hommel", "hochberg", "bonferroni", "BH", "BY", "fdr")){
  
  
  
  xfl <- x[,factorDesired]
  levels <- rev(levels)
  colnames(x)[factorDesired] <- "condition"
  x$condition <- levels[as.factor(xfl)]
  
  
  
  # Check if there is block
  if (is.null(block)) {
    
    
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-5)]
      CRDfactors <- paste(factors, collapse = "*")
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log10(x$Etarget)*x$Cttarget)-(log10(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      
      factors <- colnames(x)[1:(ncol(x)-7)]
      CRDfactors <- paste(factors, collapse = "*")
      colnames(x)[ncol(x)-6] <- "rep"
      colnames(x)[ncol(x)-5] <- "Etarget"
      colnames(x)[ncol(x)-4] <- "Cttarget"
      colnames(x)[ncol(x)-3] <- "Eref"
      colnames(x)[ncol(x)-2] <- "Ctref"
      colnames(x)[ncol(x)-1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x[1:(ncol(x)-2)], wDCt = (log10(x$Etarget)*x$Cttarget)-
                        ((log10(x$Eref)*x$Ctref) + (log10(x$Eref2)*x$Ctref2))/2)
    }
    
  } else {    # if there is Block
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-6)]
      CRDfactors <- paste(factors, collapse = "*")
      colnames(x)[ncol(x)-5] <- "block"
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log10(x$Etarget)*x$Cttarget)-(log10(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      
      factors <- colnames(x)[1:(ncol(x)-8)]
      CRDfactors <- paste(factors, collapse = "*")
      colnames(x)[ncol(x)-7] <- "block"
      colnames(x)[ncol(x)-6] <- "rep"
      colnames(x)[ncol(x)-5] <- "Etarget"
      colnames(x)[ncol(x)-4] <- "Cttarget"
      colnames(x)[ncol(x)-3] <- "Eref"
      colnames(x)[ncol(x)-2] <- "Ctref"
      colnames(x)[ncol(x)-1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x[1:(ncol(x)-2)], wDCt = (log10(x$Etarget)*x$Cttarget)-
                        ((log10(x$Eref)*x$Ctref) + (log10(x$Eref2)*x$Ctref2))/2)
    }
  }
  
  
  
  
  
  # Check if there is block
  if (is.null(block)) {
    
    # If ANOVA based on factorial design was desired:
    lm0 <- stats::lm(stats::as.formula(paste("wDCt ~", CRDfactors)), x)
    factorialANOVA <- stats::anova(lm0)

  } else {
    # If ANOVA based on factorial design was desired with blocking factor:
    lm0 <- stats::lm(stats::as.formula(paste("wDCt ~",  "block +", CRDfactors)), x)
    factorialANOVA <- stats::anova(lm0)
  }
  
  
  
  # Reverse ordering of the grouping letters
  invOrder <- function(invg){
    collapsed <- paste(invg,sep="", collapse = "")
    u <- unique(strsplit(collapsed, "")[[1]])
    if(length(u) < 2){
      return( invg)
    }
    u <- u[order(u)]
    m <- matrix(nrow = NROW(invg), ncol=length(u))
    m[] <-F
    for(i in 1:length( invg)){
      s <- strsplit( invg[i],"")[[1]]
      index <- match(s, u)
      m[i, index] <- T
    }
    for(i in 1:(length(u) - 1)){
      firstColT <- match(T, m[, i])[1]
      firstT <- match(T, rowSums(m[, i:length(u)] > 0))[1]
      if(firstT < firstColT){
        colT <- match(T, m[firstT, i:length(u)])[1]
        colT <- colT + i - 1
        tmp <- m[, colT]
        m[, colT] <- m[,i]
        m[, i] <- tmp
      }
    }
    res <- vector(mode = "character", length = length("trt"))
    for(i in 1:length(invg)){
      l <- u[m[i, ]]
      res[i] <- paste(l, sep = "",collapse = "")
    }
    return(res)
  }
  
  
  
  
  # If the analysis and post-hoc testing is based on factorial
  lmf <- lm0

  FACTOR <- base::attr(stats::terms(lmf), "term.labels")[factorDesired]

  
  # Preparing final result table including letter grouping of the means
  g <-  LSD.test(lmf, FACTOR, group = T, console = F, alpha = 0.05, p.adj = p.adj)$groups
  g <- g[rev(rownames(g)),] #order the result the way you want
  g$groups <- invOrder(as.character(g$groups))
  mean <- LSD.test(lmf, FACTOR, group = T, console = F, alpha = 0.05, p.adj = p.adj)$means
  
  
  # Comparing mean pairs that also returns CI
  # Preparing final result table including letter grouping of the means
  meanPP <- LSD.test(lmf, FACTOR, group = F, console = F, alpha = 0.05, p.adj = p.adj)
  meanPairs <- meanPP$comparison
  ROWS <- rownames(meanPairs)
  diffs <- meanPairs$difference
  pval <- meanPairs$pvalue
  signif <- meanPairs$signif.
  ucl <- meanPairs$UCL
  lcl <- meanPairs$LCL
  Post_hoc_Testing <- data.frame(row.names = ROWS,
                                 FC = round(10^(-diffs), 4),
                                 pvalue = pval,
                                 signif. = signif,
                                 LCL = round(10^(-ucl), 4),
                                 UCL = round(10^(-lcl), 4))
  
  
  
  

  
  FINALDATA <- x
  POSTHUC <- Post_hoc_Testing
  
  Nrows <- length(unique(FINALDATA[,factorDesired])[-1])
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
  
  
  pfc2 <- ggplot(tableC, aes(factor(pairs, levels = rownames(tableC)), FCp)) +
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
  
  
  
  
  
  outlist2 <- list(Final_data = x,
                   lm = lmf,
                   ANOVA_factorial = factorialANOVA,
                   Table = tableC,
                   plot = pfc2)
  
  
  
  return(outlist2)
}

