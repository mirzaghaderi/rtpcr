#' @title Relative efficiency (\eqn{\Delta C_T} method) analysis using ANOVA 
#' 
#' @description Analysis of Variance of relative efficiency (\eqn{\Delta C_T} method) values based on a completely randomized design (CRD). Even there are more than a factor in the experiment, it is still possible to apply CRD analysis on the factor-level combinations as treatments. Analysis of variance based on factorial design or analysis of covariance can be performed using \code{qpcrANCOVA} function.  
#' @details The \code{qpcrANOVA} function performs analysis of variance (ANOVA) of relative efficiency (RE) values based on a completely randomized design (CRD). 
#' It is suitable when relative expression (RE) analysis between different treatment combinations 
#' (in a Uni- or multi-factorial experiment) is desired. If there are more than a factor in the experiment, 
#' it is still possible to apply CRD analysis on the factor-level combinations as treatments. 
#' For this, a column of treatment combinations is made first as a grouping factor Fold change analysis based 
#' on factorial design or analysis of covariance for the can be performed using \link{qpcrANCOVA}.
#' @author Ghader Mirzaghaderi
#' @export qpcrANOVA
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import ggplot2
#' @import lmerTest
#' @import agricolae
#' @param x A data frame consisting of condition columns, target gene efficiency (E), target Gene Ct, reference gene efficiency and reference gene Ct values, respectively. Each Ct in the data frame is the mean of technical replicates. Complete amplification efficiencies of 2 was assumed in the example data for all wells but the calculated efficienies can be used instead. NOTE: Each line belongs to a separate individual reflecting a non-repeated measure experiment).
#' @param numberOfrefGenes number of reference genes (1 or 2). Up to two reference genes can be handled.
#' @param block column name of the blocking factor (for correct column arrangement see example data.). When a qPCR experiment is done in multiple qPCR plates, variation resulting from the plates may interfere with the actual amount of gene expression. One solution is to conduct each plate as a complete randomized block so that at least one replicate of each treatment and control is present on a plate. Block effect is usually considered as random and its interaction with any main effect is not considered.
#' @param p.adj Method for adjusting p values (see p.adjust)
#' @return A list with 5 elements:
#' \describe{
#'   \item{Final_data}{The row data plus weighed delta Ct (wDCt) values.}
#'   \item{lm}{The output of linear model analysis including ANOVA tables based on factorial experiment and completely randomized design (CRD).}
#'   \item{ANOVA_factorial}{ANOVA table based on factorial arrangement}
#'   \item{ANOVA_CRD}{ANOVA table based on CRD}
#'   \item{Result}{The result table including treatments and factors, RE (Relative Expression), LCL, UCL, letter grouping and standard deviation of relative expression.}
#'   \item{Post_hoc_Test}{Post hoc test of FC (Fold Change), pvalue, significance and confidence interval (LCL, UCL).}
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
#' @examples
#'
#' # If the data include technical replicates, means of technical replicates
#' # should be calculated first using meanTech function.
#'
#' # Applying ANOVA analysis
#' qpcrANOVA(
#'      data_3factor,
#'      numberOfrefGenes = 1,
#'      p.adj = "none")
#'
#'
#' qpcrANOVA(
#'     data_2factorBlock,
#'     block = "Block",
#'     numberOfrefGenes = 1)
#'
#'



qpcrANOVA <- function(x,
                      numberOfrefGenes,
                      block = NULL,
                      p.adj = c("none","holm","hommel", 
                                "hochberg", "bonferroni", "BH", "BY", "fdr")){
  
  

  
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
  if(numberOfrefGenes == 1) {
    if (is.null(block)) {
    # Concatenate the columns using paste0
    x$T <- do.call(paste, c(x[1:(ncol(x)-6)], sep = ":"))

    lm <- stats::lm(wDCt ~ T, x)
    anovaCRD <- stats::anova(lm)
    
  } else {
    # Concatenate the columns using paste0
    x$T <- do.call(paste, c(x[1:(ncol(x)-7)], sep = ":"))
    lm <- stats::lm(wDCt ~ block + T, x)
    anovaCRD <- stats::anova(lm)
  }
  } 
  if(numberOfrefGenes == 2) {
    if (is.null(block)) {
      # Concatenate the columns using paste0
      x$T <- do.call(paste, c(x[1:(ncol(x)-8)], sep = ":"))
      
      lm <- stats::lm(wDCt ~ T, x)
      anovaCRD <- stats::anova(lm)
      
    } else {
      # Concatenate the columns using paste0
      x$T <- do.call(paste, c(x[1:(ncol(x)-9)], sep = ":"))
      lm <- stats::lm(wDCt ~ block + T, x)
      anovaCRD <- stats::anova(lm)
    }
  }
  
  
  # Preparing final result table including letter grouping of the means for T
  g <- LSD.test(lm, "T", group = T, console = F, alpha = 0.05, p.adj = p.adj)$groups
  g <- g[rev(rownames(g)),] #order the result the way you want
  g$groups <- .invOrder(as.character(g$groups))
  mean <- LSD.test(lm, "T", group = T, console = F, alpha = 0.05, p.adj = p.adj)$means
  
  
  # Comparing mean pairs that also returns CI for T
  # Preparing final result table including letter grouping of the means
  meanPP <- LSD.test(lm, "T", group = F, console = F, alpha = 0.05, p.adj = p.adj)
  meanPairs <- meanPP$comparison
  ROWS <- rownames(meanPairs)
  diffs <- meanPairs$difference
  pval <- meanPairs$pvalue
  signif <- meanPairs$signif.
  ucl <- meanPairs$UCL
  lcl <- meanPairs$LCL
  Post_hoc_Testing <- data.frame(row.names = ROWS,
                                 FC = round(2^(-diffs), 4),
                                 pvalue = pval,
                                 signif. = signif,
                                 LCL = round(2^(-ucl), 4),
                                 UCL = round(2^(-lcl), 4))
  
  
  RowNames <- rownames(mean)
  mean$RowNames <- RowNames
  

  mean <- separate(mean, RowNames, into = factors, sep = ":", remove = T)

  
  mean <- mean[order(rownames(mean)),]
  g <- g[order(rownames(g)),]
  
  bwDCt <- x$wDCt    #bwDCt <- 2^(-x$wDCt)
  sdRow <- summarise(
    group_by(data.frame(T = x$T, bwDCt = bwDCt), T),
    se = stats::sd(bwDCt/sqrt(length(bwDCt)), na.rm = TRUE))   #sd = sd(bwDCt, na.rm = TRUE))
  se <- sdRow[order(sdRow$T),]      #sd <- sdRow[order(sdRow$T),]
  
  Results <- data.frame(mean[,(ncol(mean)-2):ncol(mean)],
                        RE = round(2^(-mean$wDCt), 5),
                        LCL = round(2^(-mean$UCL), 5),
                        UCL = round(2^(-mean$LCL), 5),
                        letters = g$groups,
                        se = round(se$se, 5))     #std = round(sd$sd, 5))
  
  
  # removing additional columns!
  if(length(factors) == 1) {
    Results <- Results[, -(1:2)]
    
  } else if(length(factors) == 2) {
    Results <- Results[, -1]
    
  } else if(length(factors) == 3) {
    Results <- Results
  }
  
  
  
  xx <- x[, -(ncol(x))] # Removing the last column of T
  # rownames(Results) <- NULL # Removing rownames 
  
  outlist <- list(Final_data = xx,
                  lmCRD = lm,
                  ANOVA = anovaCRD,
                  Result = Results)  # Post_hoc_Test = Post_hoc_Testing
  
  
  return(outlist)
}
