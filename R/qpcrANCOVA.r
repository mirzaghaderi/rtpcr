#' @title Fold change (FC) analysis using ANCOVA
#' @description ANCOVA (analysis of covariance) and ANOVA (analysis of variance) can be performed using 
#' \code{qpcrANCOVA} function, for uni- or multi-factorial experiment data. This function performs FC analysis even
#' if there is only one factor (without covariate variable), although, for the data with 
#' only one factor, the analysis turns into ANOVA. The bar plot of the fold changes (FC) 
#' values along with the confidence interval is also returned by the \code{qpcrANCOVA} function. 
#' @details The \code{qpcrANCOVA} function applies both ANCOVA and ANOVA analysis to the data of a uni- or 
#' multi-factorial experiment, although for the data with 
#' only one factor, the analysis turns to ANOVA. ANCOVA is basically appropriate when the 
#' levels of a factor are 
#' also affected by uncontrolled quantitative covariate(s). 
#' For example, suppose that wDCt of a target gene in a plant is affected by temperature. The gene may 
#' also be affected by drought. Since we already know that temperature affects the target gene, we are 
#' interested to now if the gene expression is also altered by the drought levels. We can design an 
#' experiment to understand the gene behavior at both temperature and drought levels at the same time. 
#' The drought is another factor (the covariate) that may affect the expression of our gene under the 
#' levels of the first factor i.e. temperature. The data of such an experiment can be analyzed by ANCOVA 
#' or even ANOVA based on a factorial experiment using \code{qpcrANCOVA}. This function performs FC 
#' analysis even there is only one factor (without covariate or factor  variable). Bar plot of fold changes 
#' (FC) values along with the confidence interval is also returned by the 
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
#' @import lme4
#' @import emmeans
#' @param x a data frame of condition (or conditions) levels, E (efficiency), genes and Ct values. Each Ct value in the data frame is the mean of technical replicates. Please refer to the vignette for preparing your data frame correctly.
#' @param numberOfrefGenes number of reference genes. Up to two reference genes can be handled.
#' @param analysisType should be one of "ancova" or "anova".
#' @param mainFactor.column main factor for which the levels FC is compared. The remaining factors are considered as covariates.
#' @param mainFactor.level.order  a vector of main factor level names. The first level in the vector is used as reference.
#' @param mainFactor.level.edit a vector of names replacing the column names of the output plot.
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
#' @param block column name of the block if there is a blocking factor (for correct column arrangement see example data.). When a qPCR experiment is done in multiple qPCR plates, variation resulting from the plates may interfere with the actual amount of gene expression. One solution is to conduct each plate as a complete randomized block so that at least one replicate of each treatment and control is present on a plate. Block effect is usually considered as random and its interaction with any main effect is not considered.
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
#' # Data from Lee et al., 2020 
#'
#'df <- meanTech(Lee_etal2020qPCR, groups = 1:3)
#'order2 <- unique(df$DS)
#'qpcrANCOVA(df, 
#'            numberOfrefGenes = 1, 
#'            analysisType = "ancova", 
#'            mainFactor.column = 2,
#'            mainFactor.level.order = c("D7", "D12", "D15","D18"),
#'            fill = c("skyblue", "#BFEFFF"),
#'            y.axis.adjust = 0.05)
#' 
#'
#' df <- meanTech(Lee_etal2020qPCR, groups = 1:3) 
#' df2 <- df[df$factor1 == "DSWHi",][-1]
#' qpcrANCOVA(df2, 
#'           mainFactor.column = 1,
#'           mainFactor.level.order = c("D7", "D12", "D15","D18"),
#'           numberOfrefGenes = 1,
#'           analysisType = "ancova",
#'           fontsizePvalue = 5,
#'           y.axis.adjust = 0)
#'
#'

qpcrANCOVA <- function(x,
                       numberOfrefGenes,
                       analysisType = "ancova",
                       mainFactor.column,
                       mainFactor.level.order,
                       mainFactor.level.edit = "none",
                       block = NULL,
                       width = 0.5,
                       fill = "#BFEFFF",
                       y.axis.adjust = 1,
                       y.axis.by = 1,
                       letter.position.adjust = 0.1,
                       ylab = "Fold Change",
                       xlab = "Pairs",
                       fontsize = 12,
                       fontsizePvalue = 7,
                       axis.text.x.angle = 0,
                       axis.text.x.hjust = 0.5){


  
  x <- x[, c(mainFactor.column, (1:ncol(x))[-mainFactor.column])] 
  x <- x[order(match(x[,1], mainFactor.level.order)), ]
  x[,1] <- factor(x[,1], levels = mainFactor.level.order)
  
  
  
  resultx <- .addwDCt(x)
  x <- resultx$x
  factors <- resultx$factors
  
  
  
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
    lmf <- lm(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmf)
    # ANCOVA 
    formula_ANCOVA <- paste("wDCt ~", paste("as.factor(", "block",") +"), paste("as.factor(", rev(factors), ")", collapse = " + "))
    lmc <- lm(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmc)
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
      lm <- lmc
    } 
    else{
      lm <- lmf
    } 
  }
  
  
  
  
  pp1 <- emmeans(lm, colnames(x)[1], data = x)
  pp <- as.data.frame(pairs(pp1))
  pp <- pp[1:length(mainFactor.level.order)-1,]
  
  
  
  resid <- lm$residuals
  x <- data.frame(x, resid = resid)
  # calculatinf sd for pairs
  wDCt <- x$wDCt
  result <- x %>%
    group_by(x[,1]) %>%
    summarize(variance = stats::var(wDCt))
  
  meanVar = c()
  for(i in 1:nrow(result)){
    meanVar[i] = sqrt((result[1,2] + result[i,2])/2)
  }
  sd <- unlist(meanVar)
  
  
  sig <- .convert_to_character(pp$p.value)
  
  
  pp <- data.frame(pp, sd = sd[-1])
  contrast <- pp[,1]
  post_hoc_test <- data.frame(contrast, 
                              FC = 1/(10^-(pp$estimate)),
                              sd = 10^-(sd[-1]),
                              pvalue = pp$p.value,
                              sig = sig)
  
  reference <- data.frame(contrast = mainFactor.level.order[1],
                          FC = "1",
                          sd = sd[1], 
                          pvalue = 1, 
                          sig = " ")
  
  post_hoc_test <- rbind(reference, post_hoc_test)
  
  
  
  FINALDATA <- x
  tableC <- post_hoc_test
  
  
  pairs <- tableC$contrast
  sd <- tableC$sd
  FCp <- as.numeric(tableC$FC)
  significance <- tableC$sig
  
  
  if (mainFactor.level.edit == "none") {
    tableC <- tableC
  } else {
    colnames(tableC) <- mainFactor.level.edit
  }
  
  
  
  
  
  
  pfc2 <- ggplot(tableC, aes(factor(pairs, levels = contrast), FCp, fill = pairs)) +
    geom_col(col = "black", width = width) +
    geom_errorbar(aes(pairs, ymin = FCp, ymax =  FCp + sd),
                  width=0.1) +
    geom_text(aes(label = significance,
                  x = pairs,
                  y = FCp + sd + letter.position.adjust),
              vjust = -0.5, size = fontsizePvalue) +
    ylab(ylab) + xlab(xlab) +
    theme_bw()+
    theme(axis.text.x = element_text(size = fontsize, color = "black", angle = axis.text.x.angle, hjust = axis.text.x.hjust),
          axis.text.y = element_text(size = fontsize, color = "black", angle = 0, hjust = 0.5),
          axis.title  = element_text(size = fontsize)) +
    scale_y_continuous(breaks=seq(0, max(sd) + max(FCp) + y.axis.adjust, by = y.axis.by),
                       limits = c(0, max(sd) + max(FCp) + y.axis.adjust + y.axis.adjust), expand = c(0, 0)) +
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
  pfc2 <- pfc2 + guides(fill = FALSE) 
  
  
  outlist2 <- list(Final_data = x,
                   lm_ANOVA = lmf,
                   lm_ANCOVA = lmc,
                   ANOVA_table = ANOVA,
                   ANCOVA_table = ANCOVA,
                   FC_statistics_of_the_main_factor  = tableC,
                   FC_Plot_of_the_main_factor_levels = pfc2)

  return(outlist2)
}
