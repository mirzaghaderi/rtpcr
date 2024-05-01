#' @title Fold change (\eqn{\Delta \Delta C_T} method) analysis of target genes using t-test
#' 
#' @description t.test based analysis of the fold change expression for any number of target genes.
#' 
#' @details The \code{qpcrTTEST} function applies a t.test based analysis to calculate fold change (\eqn{\Delta \Delta C_T} method) expression and returns related statistics for any number of target genes along with one or two 
#' reference gene(s), that have been evaluated under control and treatment conditions. Sampling may be paired or 
#' unpaired. Paired samples in quantitative PCR refer to two sample data that are collected from one set of individuals 
#' at two different conditions, for example before and after a treatment or at two different time points. While 
#' for unpaired samples, two sets of individuals are used: one under untreated and the other set under treated 
#' condition.  Paired samples allow to compare gene expression changes within the same individual, reducing 
#' inter-individual variability. Unpaired and paired samples are commonly analyzed using unpaired and paired t-test, respectively.
#' @author Ghader Mirzaghaderi
#' @export qpcrTTEST
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import agricolae
#' @param x a data frame of 4 columns including Conditions, E (efficiency), Gene and Ct values (see example below). Biological replicates needs to be equal for all Genes. Each Ct value is the mean of technical replicates. Complete amplification efficiencies of 2 is assumed here for all wells but the calculated efficienies can be used instead.
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param numberOfrefGenes number of reference genes. Up to two reference genes can be handled.
#' @return A list of two elements:
#' \describe{
#'   \item{Row_data}{The row data including Genes and weighed delta Ct (wDCt) values.}
#'   \item{Result}{Output table including the Fold Change values, lower and upper confidence interval and the pvalues from compairing fold change between treated and non-treated conditions}
#' }
#' For more information about the test procedure and its arguments,
#' refer \code{\link[stats]{t.test}}, and \code{\link[stats]{lm}}.
#' If the residuals of the model do not follow normal distribution and variances between the two groups are not homoGene, \code{\link[stats]{wilcox.test}} procedure may be concidered
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method. Methods 25 (4). doi:10.1006/meth.2001.1262.
#'
#' Ganger, MT, Dietz GD, and Ewing SJ. 2017. A common base method for analysis
#' of qPCR data and the application of simple blocking in qPCR experiments.
#' BMC bioinformatics 18, 1-11.
#'
#' Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006.
#' Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85).
#' doi:10.1186/1471-2105-7-85.
#'
#' @examples
#'
#' # See the sample data structure
#' data_ttest
#'
#' # Getting t.test results
#' qpcrTTEST(data_ttest,
#'    paired = FALSE,
#'    var.equal = TRUE,
#'    numberOfrefGenes = 1)
#'
#'
#'
#' qpcrTTEST(Taylor_etal2019, 
#'           numberOfrefGenes = 2, 
#'           var.equal = TRUE)
#'  
#'           
#' qpcrTTESTplot(Taylor_etal2019, 
#'               numberOfrefGenes = 2, 
#'               var.equal = TRUE,
#'               y.axis.adjust = -0.7,
#'               y.axis.by = 0.5)
#'
#'




qpcrTTEST <- function(x,
                      numberOfrefGenes,
                      paired = FALSE,
                      var.equal = TRUE) {

  colnames(x)[1] <- "Condition"
  colnames(x)[2] <- "Gene"
  colnames(x)[3] <- "E"
  colnames(x)[4] <- "Ct"




  r <- nrow(x)/(2 * length(unique(x$Gene)))

  if(!all(r %% 1 == 0)) {
    stop("Error: Replicates are not equal for all Genes!")
  } else {



    x <- data.frame(x, wCt = log2(x$E) * x$Ct)

    if(numberOfrefGenes == 1) {
      x <- x
    } else {
      a <- (((2 * r) * (length(unique(x$Gene)) - 2)) + 1)
      b <- ((length(unique(x$Gene)) - 1) * 2 * r)
      mwCT <- (x$wCt[a:b] + x$wCt[(a+(2*r)):(b+(2*r))])/2
      x$wCt[a:b] <- mwCT
      x <- x[-((a+(2*r)):(b+(2*r))),]
    }

    
    
    

    
    
    GENE <- x$Gene


    levels_to_compare <- unique(GENE)[-length(unique(GENE))]
    res <- matrix(nrow = length(levels_to_compare), ncol=7)
    colnames(res) <- c("Gene", "dif", "FC", "LCL", "UCL", "pvalue", "se")
    subset <- matrix(NA, nrow = 2 * r, ncol=length(levels_to_compare))
    ttest_result <- vector("list", length(levels_to_compare))
    


    for (i in 1:length(levels_to_compare)) {
      subset[,i] <- x[GENE == levels_to_compare[i], "wCt"] - x[GENE == utils::tail(unique(GENE), 1), "wCt"]
      ttest_result[[i]] <- stats::t.test(subset[(r + 1):(2 * r), i], subset[1:r, i], paired = paired, var.equal = var.equal)
      
      res[i, ] <- c(levels_to_compare[i],
                    round(mean(subset[(r+1):(2*r),i]) - mean(subset[1:r, i]), 4),
                    round(2^-((mean(subset[(r+1):(2*r), i]) - mean(subset[1:r,i]))), 4),
                    round(2^(-ttest_result[[i]]$conf.int[2]), 4), # Lower error bar point
                    round(2^(-ttest_result[[i]]$conf.int[1]), 4), # Upper error bar point
                    round(ttest_result[[i]]$p.value, 4),
                    round(stats::sd(subset[(r+1):(2*r),i])/sqrt(r), 4))

    }
    Raw_df <- melt(subset, value.name = "wDCt")[-1]
    res <- list(Raw_data = Raw_df, Result = data.frame(res))
    return(res)
  }
}
