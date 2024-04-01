#' @title qPCR Analysis of Variance
#' @description Analysis of Variance of qPCR data
#' @details The qpcrANOVA performs ANOVA (analysis of variance) of qPCR data. It is suitable when there is a factor with more than two levels or when the are more that a condition factor in the experiment or when there is a blocking factor.
#' @author Ghader Mirzaghaderi
#' @export qpcrANOVA
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x A data frame consisting of condition columns, target gene efficiency (E), target Gene Ct, reference gene efficiency and reference gene Ct values, respectively. Each Ct in the data frame is the mean of technical replicates. Complete amplification efficiencies of 2 was assumed in the example data for all wells but the calculated efficienies can be used instead.
#' @param numberOfrefGenes number of reference genes (1 or 2). Up to two reference genes can be handled.
#' @param block column name of the blocking factor (for correct column arrangement see example data.)
#' @param p.adj Method for adjusting p values (see p.adjust)
#' @return A list with 5 elements:
#' \describe{
#'   \item{Final_data}{The final row data including weighed delta Ct (wDCt) values.}
#'   \item{lm}{The output of linear model analysis including ANOVA tables based on factorial experiment and completely randomized design (CRD).}
#'   \item{ANOVA_factorial}{ANOVA table based on factorial arrangement}
#'   \item{ANOVA_CRD}{ANOVA table based on CRD}
#'   \item{Result}{The main result table including treatments and factors, wDCt, LCL, UCL, letters and standard deviation of average relative expression data.}
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
#'      data_3factor_a,
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
                      p.adj = c("none","holm","hommel", "hochberg", "bonferroni", "BH", "BY", "fdr")){




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


    # Concatenate the columns using paste0
    x$T <- do.call(paste, c(x[1:(ncol(x)-6)], sep = ":"))
    x
    lm <- stats::lm(wDCt ~ T, x)
    anovaCRD <- stats::anova(lm)


  } else {
    # If ANOVA based on factorial design was desired with blocking factor:
    lm0 <- stats::lm(stats::as.formula(paste("wDCt ~",  "block +", CRDfactors)), x)
    factorialANOVA <- stats::anova(lm0)


    # Concatenate the columns using paste0
    x$T <- do.call(paste, c(x[1:(ncol(x)-7)], sep = ":"))
    x
    lm <- stats::lm(wDCt ~ block + T, x)
    anovaCRD <- stats::anova(lm)
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


  # Preparing final result table including letter grouping of the means
  g <- LSD.test(lm, "T", group = T, console = F, alpha = 0.05, p.adj = p.adj)$groups
  g <- g[rev(rownames(g)),] #order the result the way you want
  g$groups <- invOrder(as.character(g$groups))
  mean <- LSD.test(lm, "T", group = T, console = F, alpha = 0.05, p.adj = p.adj)$means


  # Comparing mean pairs that also returns CI
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
                                 FC = round(10^(-diffs), 4),
                                 pvalue = pval,
                                 signif. = signif,
                                 LCL = round(10^(-ucl), 4),
                                 UCL = round(10^(-lcl), 4))


  RowNames <- rownames(mean)
  mean$RowNames <- RowNames
  mean <- separate(mean, RowNames, into = factors, sep = ":", remove = T)
  mean <- mean[order(rownames(mean)),]
  g <- g[order(rownames(g)),]

  bwDCt <- 10^(-x$wDCt)
  sdRow <- summarise(
    group_by(data.frame(T = x$T, bwDCt = bwDCt), T),
    sd = sd(bwDCt, na.rm = TRUE))
  sd <- sdRow[order(sdRow$T),]

  Results <- data.frame(mean[,(ncol(mean)-2):ncol(mean)],
                        RE = round(10^(-mean$wDCt), 4),
                        LCL = round(10^(-mean$LCL), 4),
                        UCL = round(10^(-mean$UCL), 4),
                        letters = g$groups,
                        std = round(sd$sd, 4))


  # removing additional columns!
  if(length(factors) == 1) {
    Results <- Results[, -(1:2)]

  } else if(length(factors) == 2) {
    Results <- Results[, -1]

  } else if(length(factors) == 3) {
    Results <- Results
  }


  xx <- x[, -(ncol(x))] # Removing the last column of T


  outlist <- list(Final_data = xx,
                  lmCRD = lm,
                  ANOVA_factorial = factorialANOVA,
                  ANOVA_CRD = anovaCRD,
                  Result = Results,
                  Post_hoc_Test = Post_hoc_Testing)
  return(outlist)
}
