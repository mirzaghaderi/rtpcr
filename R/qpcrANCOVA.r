#' @title Fold change (\eqn{\Delta \Delta C_T} method) analysis using ANCOVA
#' 
#' @description ANCOVA (analysis of covariance) and ANOVA (analysis of variance) can be performed using 
#' 
#' \code{qpcrANCOVA} function, for uni- or multi-factorial experiment data. This function performs fold 
#' change or \eqn{\Delta \Delta C_T} method analysis for the target gene even
#' if there is only one factor (without covariate variable), although, for the data with 
#' only one factor, the analysis turns into ANOVA. The bar plot of the fold changes (FC) 
#' values along with the standard error (se) and confidence interval (ci) is also returned by 
#' the \code{qpcrANCOVA} function. 
#' 
#' @details The \code{qpcrANCOVA} function applies both ANCOVA and ANOVA analysis to the data of a uni- or 
#' multi-factorial experiment, although for the data with 
#' only one factor, the analysis turns to ANOVA. ANCOVA is basically appropriate when the 
#' levels of a factor are 
#' also affected by uncontrolled quantitative covariate(s). 
#' For example, suppose that wDCt of a target gene in a plant is affected by temperature. The gene may 
#' also be affected by drought. Since we already know that temperature affects the target gene, we are 
#' interested to know if the gene expression is also altered by the drought levels. We can design an 
#' experiment to understand the gene behavior at both temperature and drought levels at the same time. 
#' The drought is another factor (the covariate) that may affect the expression of our gene under the 
#' levels of the first factor i.e. temperature. The data of such an experiment can be analyzed by ANCOVA 
#' or even ANOVA based on a factorial experiment using \code{qpcrANCOVA}. This function performs FC 
#' analysis even there is only one factor (without covariate or factor  variable). Bar plot of fold changes 
#' (FC) values along with the pair-wise errors (square roots of pooled variances of each pair of samples) are also returned by the 
#' \code{qpcrANCOVA} function. There is also a function called \code{oneFACTORplot} which returns RE values 
#' and related plot for a one-factor-experiment with more than two levels.
#' Along with the ANCOVA, the \code{qpcrANCOVA} also performs a full model factorial analysis of variance. 
#' If there is covariate variable(s), before ANCOVA analysis, it is better to run ANOVA based on a 
#' factorial design to see if the main factor and covariate(s) interaction is significant or not. 
#' If the pvalue of the interaction effect is smaller than 0.05, then the interaction between the main factor and covariate 
#' is significant, suggesting that ANCOVA is not appropriate in this case.
#' @author Ghader Mirzaghaderi
#' @export qpcrANCOVA
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lmerTest
#' @import emmeans
#' @param x a data frame of condition(s), biological replicates, efficiency (E) and Ct values of target and reference genes. Each Ct value in the data frame is the mean of technical replicates. NOTE: Each line belongs to a separate individual reflecting a non-repeated measure experiment). Please refer to the vignette for preparing your data frame correctly.
#' @param numberOfrefGenes number of reference genes which is 1 or 2 (Up to two reference genes can be handled).
#' @param analysisType should be one of "ancova" or "anova". Default is "ancova".
#' @param mainFactor.column the factor for which FC is calculated for its levels. The remaining factors (if any) are considered as covariate(s).
#' @param mainFactor.level.order  NULL or a vector of main factor level names. If \code{NULL}, the first level of the \code{mainFactor.column} is used 
#' as reference or calibrator. If a vector of main factor levels (in any order) is specified, the first level in the vector is used as calibrator. Calibrator is the reference level or sample that all others are compared to. Examples are untreated 
#' of time 0. The FC value of the reference or calibrator level is 1 because it is not changed compared to itself.
#' If NULL, the first level of the main factor column is used as calibrator.
#' @param width a positive number determining bar width. 
#' @param fill  specify the fill color for the columns in the bar plot. If a vector of two colors is specified, the reference level is differentialy colored.
#' @param y.axis.adjust  a negative or positive value for reducing or increasing the length of the y axis.
#' @param letter.position.adjust adjust the distance between the signs and the error bars.
#' @param y.axis.by determines y axis step length
#' @param xlab  the title of the x axis
#' @param ylab  the title of the y axis
#' @param fontsize font size of the plot
#' @param fontsizePvalue font size of the pvalue labels
#' @param axis.text.x.angle angle of x axis text
#' @param axis.text.x.hjust horizontal justification of x axis text
#' @param x.axis.labels.rename a vector replacing the x axis labels in the bar plot
#' @param block column name of the block if there is a blocking factor (for correct column arrangement see example data.). When a qPCR experiment is done in multiple qPCR plates, variation resulting from the plates may interfere with the actual amount of gene expression. One solution is to conduct each plate as a complete randomized block so that at least one replicate of each treatment and control is present on a plate. Block effect is usually considered as random and its interaction with any main effect is not considered.
#' @param p.adj Method for adjusting p values
#' @param errorbar Type of error bar, can be \code{se} or \code{ci}.
#' @return A list with 7 elements:
#' \describe{
#'   \item{Final_data}{Input data frame plus the weighted Delat Ct values (wDCt)}
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
#'
#'  qpcrANCOVA(data_1factor, 
#'             numberOfrefGenes = 1,
#'             mainFactor.column = 1,
#'             fill = c("#CDC673", "#EEDD82"),
#'             fontsizePvalue = 5,
#'             y.axis.adjust = 0.1)
#'
#' qpcrANCOVA(data_2factor, 
#'            numberOfrefGenes = 1,
#'            mainFactor.column = 2,
#'            mainFactor.level.order = c("D0", "D1", "D2"),
#'            fill = c("#79CDCD", "#B4EEB4"),
#'            analysisType = "ancova",
#'            fontsizePvalue = 5,
#'            y.axis.adjust = 0.4)
#'
#'
#' # Data from Lee et al., 2020 
#' # Here, the data set contains technical replicates so 
#' # we get mean of technical reps first:
#' df <- meanTech(Lee_etal2020qPCR, groups = 1:3)
#' order <- rev(unique(df$DS))
#' qpcrANCOVA(df, 
#'            numberOfrefGenes = 1, 
#'            analysisType = "ancova", 
#'            mainFactor.column = 2,
#'            mainFactor.level.order = order,
#'            fill = c("skyblue", "#BFEFFF"),
#'            y.axis.adjust = 0.05)
#' 
#'
#' df <- meanTech(Lee_etal2020qPCR, groups = 1:3) 
#' df2 <- df[df$factor1 == "DSWHi",][-1]
#' qpcrANCOVA(df2, 
#'            mainFactor.column = 1,
#'            mainFactor.level.order = c("D7", "D12", "D15","D18"),
#'            numberOfrefGenes = 1,
#'            analysisType = "ancova",
#'            fontsizePvalue = 5,
#'            y.axis.adjust = 0.1)
#'
#'
#' qpcrANCOVA(data_2factorBlock,
#'            numberOfrefGenes = 1,
#'            mainFactor.column = 1, 
#'            mainFactor.level.order = c("S", "R"),
#'            block = "block", 
#'            fill = c("#CDC673", "#EEDD82"),
#'            analysisType = "ancova",
#'            fontsizePvalue = 7,
#'            y.axis.adjust = 0.1, 
#'            width = 0.35)
#'
#' addline_format <- function(x,...){gsub('\\s','\n',x)}
#' order <- unique(data_2factor$Drought)
#' qpcrANCOVA(data_1factor,
#'    numberOfrefGenes = 1,
#'    mainFactor.column = 1,
#'    mainFactor.level.order = c("L1","L2","L3"),
#'    width = 0.5,
#'    fill = c("skyblue","#79CDCD"),
#'    y.axis.by = 1,
#'    letter.position.adjust = 0,
#'    y.axis.adjust = 1,
#'    ylab = "Fold Change",
#'    fontsize = 12,
#'    x.axis.labels.rename = addline_format(c("Control", 
#'                                          "Treatment_1 vs Control", 
#'                                          "Treatment_2 vs Control")))
#'                                                        
#'                                                        



qpcrANCOVA <- function(x,
                       numberOfrefGenes,
                       analysisType = "ancova",
                       mainFactor.column,
                       mainFactor.level.order = NULL,
                       block = NULL,
                       width = 0.5,
                       fill = "#BFEFFF",
                       y.axis.adjust = 1,
                       y.axis.by = 1,
                       letter.position.adjust = 0.1,
                       ylab = "Fold Change",
                       xlab = "none",
                       fontsize = 12,
                       fontsizePvalue = 7,
                       axis.text.x.angle = 0,
                       axis.text.x.hjust = 0.5,
                       x.axis.labels.rename = "none",
                       p.adj = "none",
                       errorbar = "se"){


  
  x <- x[, c(mainFactor.column, (1:ncol(x))[-mainFactor.column])] 
  
  
  if (is.null(mainFactor.level.order)) {
    mainFactor.level.order <- unique(x[,1])
    calibrartor <- x[,1][1]
    warning(paste("The", calibrartor, "level was used as calibrator."))
  } else if (any(is.na(match(unique(x[,1]), mainFactor.level.order))) == TRUE){
    stop("The `mainFactor.level.order` doesn't match main factor levels.")
  } else {
    x <- x[order(match(x[,1], mainFactor.level.order)), ]
    x[,1] <- factor(x[,1], levels = mainFactor.level.order)
  }
  
  
  # The data frame doesn't have ? columns. Please refer to the vignette to ensure that our data is properly structured.
  
  
  if (is.null(block)) {
    
    
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-5)]
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      
      factors <- colnames(x)[1:(ncol(x)-7)]
      colnames(x)[ncol(x)-6] <- "rep"
      colnames(x)[ncol(x)-5] <- "Etarget"
      colnames(x)[ncol(x)-4] <- "Cttarget"
      colnames(x)[ncol(x)-3] <- "Eref"
      colnames(x)[ncol(x)-2] <- "Ctref"
      colnames(x)[ncol(x)-1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-
                        ((log2(x$Eref)*x$Ctref) + (log2(x$Eref2)*x$Ctref2))/2)
    }
    
  } else {
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-6)]
      colnames(x)[ncol(x)-5] <- "block"
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
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
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-
                        ((log2(x$Eref)*x$Ctref) + (log2(x$Eref2)*x$Ctref2))/2)
    }
  }
  
  
  
  # Check if there is block
  if (is.null(block)) {
    
    # ANOVA based on factorial design
    formula_ANOVA <- paste("wDCt ~", paste("as.factor(", factors, ")", collapse = " * "))
    lmf <- lm(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmf)
    # ANCOVA 
    formula_ANCOVA <- paste("wDCt ~", paste("as.factor(", rev(factors), ")", collapse = " + "))
    lmc <- lm(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmc)
    #rownames(ANCOVA) <- as.vector(cat(paste0('"', rev(factors), '"'), '"Residuals"'))
    
  } else {
    # If ANOVA based on factorial design was desired with blocking factor:
    formula_ANOVA <- paste("wDCt ~", paste("as.factor(", "block",") +"), paste("as.factor(", factors, ")", collapse = " * "))
    lmfb <- lm(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmfb)
    # ANCOVA 
    formula_ANCOVA <- paste("wDCt ~", paste("as.factor(", "block",") +"), paste("as.factor(", rev(factors), ")", collapse = " + "))
    lmcb <- lm(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmcb)
  }
  
  
  
  
  
  # Type of analysis: ancova or anova
  if (is.null(block)) {
    if(analysisType == "ancova") {
      lm <- lmc
    } 
    else{
      lm <- lmf
    }
  } else {
    if(analysisType == "ancova") {
      lm <- lmcb
    } 
    else{
      lm <- lmfb
    } 
  }
  
  
  
  
  pp1 <- emmeans(lm, colnames(x)[1], data = x, adjust = p.adj)
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  pp3 <- pp2[1:length(mainFactor.level.order)-1,]
  ci <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(unique(x[,1]))-1,]
  pp <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)

  

  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(factor = x[,1], bwDCt = bwDCt), x[,1]),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))  
  
  
  sig <- .convert_to_character(pp$p.value)
  contrast <- pp$contrast
  post_hoc_test <- data.frame(contrast, 
                              FC = round(1/(2^-(pp$estimate)), 7),
                              pvalue = pp$p.value,
                              sig = sig,
                              LCL = 1/(2^-pp$lower.CL),
                              UCL = 1/(2^-pp$upper.CL),
                              se = se$se[-1])
  
  reference <- data.frame(contrast = mainFactor.level.order[1],
                          FC = 1,
                          pvalue = 1, 
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])
  
  tableC <- rbind(reference, post_hoc_test)
  
  FINALDATA <- x
  
  tableC$contrast <- sapply(strsplit(tableC$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))
  
  if(any(x.axis.labels.rename == "none")){
    tableC
  }else{
    tableC$contrast <- x.axis.labels.rename
  }
  
  
  
  
  tableC$contrast <- factor(tableC$contrast, levels = unique(tableC$contrast))
  contrast <- tableC$contrast
  LCL <- tableC$LCL
  UCL <- tableC$UCL
  FCp <- as.numeric(tableC$FC)
  significance <- tableC$sig
  se <- tableC$se
  
  

  
  pfc2 <- ggplot(tableC, aes(contrast, FCp, fill = contrast)) +
    geom_col(col = "black", width = width)
  
  
  
  if(errorbar == "ci") {
    pfc2 <- pfc2 +
      geom_errorbar(aes(contrast, ymin = LCL, ymax =  UCL), width=0.1) +
      geom_text(aes(label = significance, x = contrast,
                           y = UCL + letter.position.adjust),
                       vjust = -0.5, size = fontsizePvalue)
  } else if(errorbar == "se") {
    pfc2 <- pfc2 +
      geom_errorbar(aes(contrast, ymin = FCp, ymax =  FCp + se), width=0.1) +
      geom_text(aes(label = significance, x = contrast,
                           y = FCp + se + letter.position.adjust),
                       vjust = -0.5, size = fontsizePvalue)
    }
    
    
  pfc2 <- pfc2 +
    ylab(ylab) +
    theme_bw()+
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize)) +
    scale_y_continuous(breaks=seq(0, max(FCp) + max(se)  + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(FCp) + max(se) + y.axis.adjust), expand = c(0, 0)) +
    theme(legend.text = element_text(colour = "black", size = fontsize),
          legend.background = element_rect(fill = "transparent"))
  
  
  if(length(fill) == 2) {
    pfc2 <- pfc2 +
             scale_fill_manual(values = c(fill[1], rep(fill[2], nrow(tableC)-1)))
  } 
  if (length(fill) == 1) {
    pfc2 <- pfc2 +
      scale_fill_manual(values = rep(fill, nrow(tableC)))
  }
  
  pfc2 <- pfc2 + guides(fill = "none") 
  
  
  if(xlab == "none"){
    pfc2 <- pfc2 + 
      labs(x = NULL)
  }else{
    pfc2 <- pfc2 +
    xlab(xlab)
  }
  
  
  
  #changing as.factor(x) to x

  if (is.null(block)) {
    lmf$coefficients <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", lmf$coefficients)
    lmf$coefficients <- gsub(":as factor", ":", lmf$coefficients)
    
    rownames(ANOVA) <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", rownames(ANOVA))
    rownames(ANOVA) <- gsub(":as factor", ":", rownames(ANOVA))
    
    lmc$coefficients <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", lmc$coefficients)
    lmc$coefficients <- gsub(":as factor", ":", lmc$coefficients)
    
    rownames(ANCOVA) <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", rownames(ANCOVA))
    rownames(ANCOVA) <- gsub(":as factor", ":", rownames(ANCOVA))
    lm_ANOVA <- lmf
    lm_ANCOVA <- lmc
  } else {
    lmfb$coefficients <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", lmfb$coefficients)
    lmfb$coefficients <- gsub(":as factor", ":", lmfb$coefficients)
    
    rownames(ANOVA) <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", rownames(ANOVA))
    rownames(ANOVA) <- gsub(":as factor", ":", rownames(ANOVA))
    
    lmcb$coefficients <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", lmcb$coefficients)
    lmcb$coefficients <- gsub(":as factor", ":", lmcb$coefficients)
    
    rownames(ANCOVA) <- gsub("as\\.factor\\(([^)]+)\\)", "\\1", rownames(ANCOVA))
    rownames(ANCOVA) <- gsub(":as factor", ":", rownames(ANCOVA))
    lm_ANOVA <- lmfb
    lm_ANCOVA <- lmcb
  }
  
  outlist2 <- list(Final_data = x,
                   lm_ANOVA = lm_ANOVA,
                   lm_ANCOVA = lm_ANCOVA,
                   ANOVA_table = ANOVA,
                   ANCOVA_table = ANCOVA,
                   FC_statistics_of_the_main_factor  = tableC,
                   FC_Plot_of_the_main_factor_levels = pfc2)

  return(outlist2)
}
