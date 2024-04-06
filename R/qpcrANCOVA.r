#' @title ANCOVA and ANOVA on a factorial design plus bar plot of FC
#' @description ANCOVA (analysis of covariance) and ANOVA (analysis of variance) can be performed using 
#' \code{qpcrANCOVA} function, for uni- or multi-factorial experiment data. Bar plot of the fold changes (FC) 
#' values along with the 95\% confidence interval is also returned by the \code{qpcrANCOVA} function. 
#' @details The \code{qpcrANCOVA} function applies both ANCOVA and ANOVA analysis to the data of a uni- or 
#' multi-factorial experiment. ANCOVA is basically suitable when the levels of a factor are 
#' also affected by uncontrolled quantitative covariate(s). 
#' For example, suppose that wDCt of a target gene in a plant is affected by temperature. The gene may 
#' also be affected by drought. Since we already know that temperature affects the target gene, we are 
#' interested to now if the gene expression is also altered by the drought levels. We can design an 
#' experiment to understand the gene behavior at both temperature and drought levels at the same time. 
#' The drought is another factor (the covariate) that may affect the expression of our gene under the 
#' levels of the first factor i.e. temperature. The data of such an experiment can be analyzed by ANCOVA 
#' or even ANOVA based on a factorial experiment using \code{qpcrANCOVA} function. Bar plot of fold changes (FC) values along with the 95\% confidence interval is also returned by the 
#' \code{qpcrANCOVA} function. There is also a function called \code{oneFACTORplot} which returns RE values 
#' and related plot for a one-factor-experiment with more than two levels.
#' Along with the ANCOVA, the \code{qpcrANCOVA} also performs a full model factorial analysis of variance. 
#' If there is covariate variable(s), before ANCOVA analysis, it is better to run ANOVA based on a 
#' factorial design to see if the main factor and covariate(s) interaction is significant or not. 
#' If the pvalue of the interation effect is smaller than 0.05, then the interaction between the main factor and covariate 
#' is significant, suggesting that ANCOVA is not appropriate in this case.
#' @author Ghader Mirzaghaderi
#' @export qpcrANCOVA
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x a data frame of condition (or conditions) levels, E (efficiency), genes and Ct values. Each Ct value in the data frame is the mean of technical replicates. Please refer to the vignette for preparing your data frame correctly.
#' @param numberOfrefGenes number of reference genes. Up to two reference genes can be handled.
#' @param analysisType should be one of "ancova" or "anova".
#' @param main.factor.column main factor for which the levels FC is compared. The remaining factors are considered as covariates.
#' @param level.names  a vector determining level names in the x axis on the plot.
#' @param levels a numeric vector corresponding to the main factor levels. The first number is used as control level in computing the FC value of the other levels.
#' @param width a positive number determining bar width.
#' @param fill  specify the fill color for the columns of the bar plot.
#' @param y.axis.adjust  a negative or positive value for reducing or increasing the length of the y axis.
#' @param letter.position.adjust adjust the distance between the signs and the error bars.
#' @param y.axis.by determines y axis step length
#' @param xlab  the title of the x axis
#' @param ylab  the title of the y axis
#' @param fontsize font size of the plot
#' @param fontsizePvalue font size of the pvalue labels
#' @param axis.text.x.angle angle of x axis text
#' @param axis.text.x.hjust horizontal justification of x axis text
#' @param block column name of the block if there is a blocking factor (for correct column arrangement see example data.)
#' @param p.adj method for adjusting p values (see \code{p.adjust})
#' @return A list with 2 elements:
#' \describe{
#'   \item{Final_data}{}
#'   \item{lm_ANOVA}{lm of factorial analysis-tyle}
#'   \item{lm_ANCOVA}{lm of ANCOVA analysis-type}
#'   \item{ANOVA_table}{ANOVA table}
#'   \item{ANCOVA_table}{ANCOVA table}
#'   \item{FC Table}{Table of FC values, significance and confidence limits for the main factor levels.}
#'   \item{Bar plot of FC values}{Bar plot of the fold change values for the main factor levels.}
#' }
#' 
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method. Methods 25 (4). doi:10.1006/meth.2001.1262.
#'
#' Ganger, MT, Dietz GD, and Ewing SJ. 2017. A common base method for analysis of qPCR data
#' and the application of simple blocking in qPCR experiments. BMC bioinformatics 18, 1-11.
#'
#' Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006.
#' Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85). doi:10.1186/1471-2105-7-85.
#' 
#' 
#' 
#' @examples
#'
#' # See sample data
#' data_2factor
#'
#'
#' qpcrANCOVA(data_2factor, 
#'            numberOfrefGenes = 1, 
#'            analysisType = "ancova", 
#'            main.factor.column = 1,
#'            levels = c(1, 2))
#'            
#'
#' qpcrANCOVA(data_2factorBlock, 
#'            numberOfrefGenes = 1,
#'            block = "block",  
#'            main.factor.column = 2,
#'            levels = c(3, 2, 1))
#'            
#'            
#'  qpcrANCOVA(data_1factor, 
#'            numberOfrefGenes = 1,
#'            analysisType = "ancova", 
#'            main.factor.column = 1,
#'            levels = 3:1)
#'            
#'            
#' 

qpcrANCOVA <- function(x,
                       numberOfrefGenes,
                       block = NULL,
                       analysisType = "ancova",
                       main.factor.column,
                       levels,
                       level.names = "none",
                       width = 0.5,
                       fill = "skyblue",
                       y.axis.adjust = 1,
                       y.axis.by = 1,
                       letter.position.adjust = 0.1,
                       ylab = "Average Fold Change",
                       xlab = "Pairs",
                       fontsize = 12,
                       fontsizePvalue = 7,
                       axis.text.x.angle = 0,
                       axis.text.x.hjust = 0.5,
                       p.adj = c("none","holm","hommel", "hochberg", "bonferroni", "BH", "BY", "fdr")){
  
  
  x <- x[, c(main.factor.column, (1:ncol(x))[-main.factor.column])]      
  x <- x[order(x[,1]),]
  xfl <- x[,1]
  levels <- rev(levels)
  #colnames(x)[1] <- "condition"
  x[,1] <- levels[as.factor(xfl)]
  
  
  # Check if there is block
  if (is.null(block)) {
    
    
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-5)]
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log10(x$Etarget)*x$Cttarget)-(log10(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      
      factors <- colnames(x)[1:(ncol(x)-7)]
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
      colnames(x)[ncol(x)-5] <- "block"
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log10(x$Etarget)*x$Cttarget)-(log10(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      
      factors <- colnames(x)[1:(ncol(x)-8)]
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
    
    # ANOVA based on factorial design
    formula_ANOVA <- paste("wDCt ~", paste("as.factor(", factors, ")", collapse = " * "))
    lmf <- lm(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmf)
    formula_ANCOVA <- paste("wDCt ~", paste("as.factor(", rev(factors), ")", collapse = " + "))
    lmc <- lm(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmc)
    
  } else {
    # If ANOVA based on factorial design was desired with blocking factor:
    formula_ANOVA <- paste("wDCt ~", paste("as.factor(", "block",") +"), paste("as.factor(", factors, ")", collapse = " * "))
    lmf <- lm(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmf)
    formula_ANCOVA <- paste("wDCt ~", paste("as.factor(", "block",") +"), paste("as.factor(", rev(factors), ")", collapse = " + "))
    lmc <- lm(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmc)
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
  
  
  
  
  # Type of analysis: ancova or anova
  if(analysisType == "ancova") {
    FACTOR <- rev(base::attr(stats::terms(lmc), "term.labels"))[1]
    lm <- lmc
  } 
  else{
    FACTOR <- base::attr(stats::terms(lmf), "term.labels")[1]
    lm <- lmf
  }
  
  
  
  
  # Preparing final result table including letter grouping of the means
  g <-  LSD.test(lm, FACTOR, group = T, console = F, alpha = 0.05, p.adj = p.adj)$groups
  g <- g[rev(rownames(g)),] #order the result the way you want
  g$groups <- invOrder(as.character(g$groups))
  mean <- LSD.test(lm, FACTOR, group = T, console = F, alpha = 0.05, p.adj = p.adj)$means
  
  
  # Comparing mean pairs that also returns CI
  # Preparing final result table including letter grouping of the means
  meanPP <- LSD.test(lm, FACTOR, group = F, console = F, alpha = 0.05, p.adj = p.adj)
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
  
  #Nrows <- length(unique(FINALDATA[,1])[-1])
  #withControl  <- POSTHUC[1:Nrows,]
  #withControl
  
  
  
  tableC <- POSTHUC
  # default level names of add level.names
  if (any(level.names == "none")) {
    rownames(tableC) <- rownames(tableC)
  } else {
    rownames(tableC) <- level.names
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
              vjust = -0.5, size = fontsizePvalue) +
    ylab(ylab) + xlab(xlab) +
    theme_bw()+
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize)) +
    scale_y_continuous(breaks=seq(0, max(UCLp) + max(FCp) + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(UCLp) + max(FCp) + y.axis.adjust + y.axis.adjust), expand = c(0, 0)) +
    theme(legend.text = element_text(colour = "black", size = fontsize),
          legend.background = element_rect(fill = "transparent"))
  
  
  
  
  
  outlist2 <- list(Final_data = x,
                   lm_ANOVA = lmf,
                   lm_ANCOVA = lmc,
                   ANOVA_table = ANOVA,
                   ANCOVA_table = ANCOVA,
                   Table  = tableC,
                   Plot = pfc2)
  
  names(outlist2)[6] <- "Fold change statistics for the main factor:"
  names(outlist2)[7] <- "Bar plot of the fold change values for the main factor levels:"
  
  return(outlist2)
}
