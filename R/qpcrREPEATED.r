#' @title Repeated measure analysis
#' @description Fold change (FC) analysis of observations repeatedly taken over time
#' \code{qpcrREPEATED} function, for Repeated measure analysis of uni- or multi-factorial experiment data. The bar plot of the fold changes (FC) 
#' values along with the standard error (se) or confidence interval (ci) is also returned by the \code{qpcrREPEATED} function. 
#' @details The \code{qpcrREPEATED} function prtforms Fold change (FC) analysis of observations repeatedly taken over time. It 
#' performs a full model factorial analysis of variance to the data and returns FC values
#' along with confidence interval and standard error for the FC values.
#' @author Ghader Mirzaghaderi
#' @export qpcrREPEATED
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import emmeans
#' @param x a data frame of condition(s), biological replicates, efficiency (E) and Ct values of target and reference genes. Each Ct value in the data frame is the mean of technical replicates. NOTE: Each line belongs to a separate individual reflecting a non-repeated measure experiment). Please refer to the vignette for preparing your data frame correctly.
#' @param numberOfrefGenes number of reference genes which is 1 or 2 (Up to two reference genes can be handled).
#' as reference or calibrator which is the reference level or sample that all others are compared to. Examples are untreated 
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
#' @param block column name of the block if there is a blocking factor (for correct column arrangement see example data.). 
#' Block effect is usually considered as random and its interaction with any main effect is not considered.
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
#' 
#' 
#' @examples
#' 
#' a <- qpcrREPEATED(data_repeated_measure_1,
#'             numberOfrefGenes = 1,
#'             block = NULL)
#'
#' qpcrREPEATED(data_repeated_measure_2,
#'              numberOfrefGenes = 1,
#'              block = NULL)
#'                                                        
#'                                                        




qpcrREPEATED <- function(x,
                         numberOfrefGenes,
                         block = NULL,
                         width = 0.5,
                         fill = "#BFEFFF",
                         y.axis.adjust = 1,
                         y.axis.by = 1,
                         ylab = "Fold Change",
                         xlab = "none",
                         fontsize = 12,
                         fontsizePvalue = 7,
                         axis.text.x.angle = 0,
                         axis.text.x.hjust = 0.5,
                         x.axis.labels.rename = "none",
                         letter.position.adjust = 0,
                         p.adj = "none",
                         errorbar = "se"){
  

  
  
    if (is.null(block)) {
    
    
    if(numberOfrefGenes == 1) {
      id <- colnames(x)[1]
      if((ncol(x)-5) <= 1){
        factors = NULL
      } else {
        factors <- colnames(x)[2:(ncol(x)-5)]
      }
      colnames(x)[ncol(x)-4] <- "time"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      id <- colnames(x)[1]
      factors <- colnames(x)[2:(ncol(x)-7)]
      colnames(x)[ncol(x)-6] <- "time"
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
      id <- colnames(x)[1]
      factors <- colnames(x)[2:(ncol(x)-6)]
      colnames(x)[ncol(x)-5] <- "block"
      colnames(x)[ncol(x)-4] <- "time"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      id <- colnames(x)[1]
      factors <- colnames(x)[2:(ncol(x)-8)]
      colnames(x)[ncol(x)-7] <- "block"
      colnames(x)[ncol(x)-6] <- "time"
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
  
  
  # converting columns 1 to time as factor
  x[, 1:which(names(x) == "time")] <- lapply(x[, 1:which(names(x) == "time")], factor)
  
  
  # Check if there is block
  if (is.null(block)) {
    if((ncol(x)-5) <= 2){
      formula <- wDCt ~ time + (1 | id)
    } else {
      formula <- paste("wDCt ~", paste("time"," *"), paste(factors, collapse = " * "), "+ (1 | id)")
      #formula <- paste("wDCt ~", paste("as.factor(","time",") *"), paste("as.factor(", factors, ")", collapse = " * "), "+ (1 | id)")
      }
    } else {
      #x <- x[, c(match("block", names(x)), (1:ncol(x))[-match("block", names(x))])]
      if((ncol(x)-6) <= 2){
      formula <- wDCt ~ time + (1|id) + (1|block/id)
      } else {
      formula <- paste("wDCt ~ ", paste("as.factor(", "time",") *"), paste("as.factor(", factors, ")", collapse = " * "), "+ (1 | id) + (1|block/id)")
      }
    }
  lm <- lmer(formula, data = x)
  ANOVA <- stats::anova(lm) 
  ANOVA
  
  
  
  pp1 <- emmeans(lm, "time", data = x, adjust = p.adj)
  pp <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  pp <- pp[1:length(unique(x$time))-1,]
  
  
  
  # Preparing t-test results
  t_test_results <- list()
  
  for (i in 1:length(x$time)) {
    level_data <- subset(x, x$time == x$time[i])$wDCt
    t_test_result <- stats::t.test(level_data, subset(x, x$time == x$time[1])$wDCt)
    t_test_results[[paste("t_test_result_", x$time[i], "_vs_", x$time[1])]] <- t_test_result
  }
  
  confidence_intervals <- data.frame(
    Comparison = sapply(names(t_test_results), function(x) gsub("t_test_result_", "", x)),
    CI_lower = sapply(t_test_results, function(x) x$conf.int[1]),
    CI_upper = sapply(t_test_results, function(x) x$conf.int[2]),
    df = sapply(t_test_results, function(x) x$parameter),
    p.value = sapply(t_test_results, function(x) x$p.value))
  
  CI <- data.frame(Comparison = confidence_intervals$Comparison,
                   LCL = 2^-confidence_intervals$CI_upper,
                   UCL = 2^-confidence_intervals$CI_lower,
                   df = confidence_intervals$df,
                   p.value = confidence_intervals$p.value)
  
  
  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(factor = x$time, bwDCt = bwDCt), x$time),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))  
  
  
  sig <- .convert_to_character(pp$p.value)
  
  
  
  
  
  contrast <- pp[,1]
  post_hoc_test <- data.frame(contrast, 
                              FC = round(1/(2^-(pp$estimate)), 4),
                              pvalue = round(pp$p.value, 4),
                              sig = sig,
                              LCL = CI[-1,]$LCL,
                              UCL = CI[-1,]$UCL,
                              se = se$se[-1])
  
  reference <- data.frame(contrast = "time1",
                          FC = "1",
                          pvalue = 1, 
                          sig = " ",
                          LCL = CI[1,2],
                          UCL = CI[1,3],
                          se = se$se[1])
  
  post_hoc_test <- rbind(reference, post_hoc_test)
  
  
  
  FINALDATA <- x
  tableC <- post_hoc_test
  
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
  
  

  
  outlist2 <- list(Final_data = x,
                   lm = lm,
                   ANOVA_table = ANOVA,
                   FC_statistics_of_the_main_factor  = tableC,
                   FC_Plot = pfc2)
  
  return(outlist2)
  
  
  
}
