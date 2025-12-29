#' @title ΔΔCt method t-test analysis
#'
#' @description
#' The \code{TTEST_DDCt} function performs fold change expression analysis based on
#' the ΔΔCt method using Student's t-test. It supports analysis
#' of one or more target genes evaluated under two experimental conditions
#' (e.g. control vs treatment).
#'
#' @details
#' Relative expression values are computed using one or more reference genes for normalization.
#' Both paired and unpaired experimental designs are supported.
#'
#' Paired samples in quantitative PCR refer to measurements collected from the same
#' individuals under two different conditions (e.g. before vs after treatment),
#' whereas unpaired samples originate from different individuals in each condition.
#' Paired designs allow within-individual comparisons and typically reduce
#' inter-individual variability.
#'
#' The function returns numerical summaries as well as bar plots based on either
#' relative expression (RE) or log2 fold change (log2FC).
#'
#' @author Ghader Mirzaghaderi
#'
#' @export
#'
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#'
#' @param x
#' A data frame containing experimental conditions, biological replicates, and
#' amplification efficiency and Ct values for target and reference genes.
#' The number of biological replicates must be equal across genes. If this 
#' is not true, or there are \code{NA} values use \code{ANODA_DDCt} function 
#' for independent samples or \code{REPEATED_DDCt} for paired samples.
#' See the package vignette for details on the required data structure.
#'
#' @param paired
#' Logical; if \code{TRUE}, a paired t-test is performed.
#'
#' @param var.equal
#' Logical; if \code{TRUE}, equal variances are assumed and a pooled variance
#' estimate is used. Otherwise, Welch's t-test is applied.
#'
#' @param numberOfrefGenes
#' Integer specifying the number of reference genes used for normalization.
#'
#' @param Factor.level.order
#' Optional character vector specifying the order of factor levels.
#' If \code{NULL}, the first level of the factor column is used as the calibrator.
#'
#' @param p.adj
#' Method for p-value adjustment. One of
#' \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#' \code{"BH"}, \code{"BY"}, \code{"fdr"}, or \code{"none"}. See \code{\link[stats]{p.adjust}}.
#'
#' @param order
#' Optional character vector specifying the order of genes in the output plot.
#'
#' @param plotType
#' Plot scale to use: \code{"RE"} for relative expression or
#' \code{"log2FC"} for log2 fold change.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{Result}{Table containing RE values, log2FC, p-values, significance codes,
#'   confidence intervals, standard errors, and lower/upper SE limits.}
#'   \item{RE_Plot}{Bar plot of relative expression values.}
#'   \item{log2FC_Plot}{Bar plot of log2 fold change values.}
#' }
#'
#' @references
#' Livak, K. J. and Schmittgen, T. D. (2001).
#' Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR
#' and the Double Delta CT Method.
#' \emph{Methods}, 25(4), 402–408.
#' doi:10.1006/meth.2001.1262
#'
#' Ganger, M. T., Dietz, G. D., and Ewing, S. J. (2017).
#' A common base method for analysis of qPCR data and the application of
#' simple blocking in qPCR experiments.
#' \emph{BMC Bioinformatics}, 18, 1–11.
#'
#' Yuan, J. S., Reed, A., Chen, F., and Stewart, N. (2006).
#' Statistical Analysis of Real-Time PCR Data.
#' \emph{BMC Bioinformatics}, 7, 85.
#'
#' @examples
#' # Example data structure
#' data1 <- read.csv(system.file("extdata", "data_1factor_one_ref.csv", package = "rtpcr"))
#'
#' # Unpaired t-test
#' TTEST_DDCt(
#'   data1,
#'   paired = FALSE,
#'   var.equal = TRUE,
#'   numberOfrefGenes = 1)
#'
#' # With amplification efficiencies
#' data2 <- read.csv(system.file("extdata", "data_1factor_one_ref_Eff.csv", package = "rtpcr"))
#'
#' TTEST_DDCt(
#'   data2,
#'   paired = FALSE,
#'   var.equal = TRUE,
#'   numberOfrefGenes = 1)
#'
#' # Two reference genes
#' data3 <- read.csv(system.file("extdata", "data_1factor_Two_ref.csv", package = "rtpcr"))
#' TTEST_DDCt(
#'   data3,
#'   numberOfrefGenes = 2,
#'   var.equal = TRUE,
#'   p.adj = "BH")



TTEST_DDCt <- function(x,
                       numberOfrefGenes,
                       Factor.level.order = NULL,
                       paired = FALSE, 
                       var.equal = TRUE, 
                       p.adj = "none",
                       order = "none", 
                       plotType = "RE") {
  

  ## Basic checks
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame")
  }
  if (!is.numeric(numberOfrefGenes) || length(numberOfrefGenes) != 1 || numberOfrefGenes < 1) {
    stop("`numberOfrefGenes` must be a single numeric value >= 1")
  }
  

  
  ## Factor level handling
  if (is.null(Factor.level.order)) {
    x[,1] <- factor(x[,1], levels = unique(x[,1]))
    Factor.level.order <- unique(x[,1])
    calibrartor <- x[,1][1]
    on.exit(cat(structure(
      paste("*** The", calibrartor, "level was used as calibrator.\n")
    )))
  } else if (any(is.na(match(unique(x[,1]), Factor.level.order)))) {
    stop("The `Factor.level.order` doesn't match factor levels.")
  } else {
    x <- x[order(match(x[,1], Factor.level.order)), ]
    x[,1] <- factor(x[,1], levels = Factor.level.order)
    calibrartor <- x[,1][1]
    on.exit(cat(structure(
      paste("*** The", calibrartor, "level was used as calibrator.\n")
    )))
  }
  

  ## Number of targets
  y <- (ncol(x) - ((numberOfrefGenes * 2) + 2)) / 2
  
  cat(sprintf(
    "*** %s target(s) using %s reference gene(s) was analysed!\n",
    y, numberOfrefGenes
  ))
  

  ## Wide to long
  x <- .wide_to_long(x)
  colnames(x)[1:4] <- c("Condition", "Gene", "E", "Ct")
  
  default.order <- unique(x$Gene)[-length(unique(x$Gene))]
  
  r <- nrow(x) / (2 * length(unique(x$Gene)))
  if (!all(r %% 1 == 0)) {
    stop("Error: Replicates are not equal for all Genes, or there are more than two conditions!")
  }
  

  ## Weighted Ct
  x$wCt <- log2(x$E) * x$Ct
  

  
  ## GENERALIZED reference gene handling (ANY number ≥ 1)
  genes <- unique(x$Gene)
  nGenes <- length(genes)
  
  if (numberOfrefGenes >= nGenes) {
    stop("Number of reference genes must be smaller than total number of genes")
  }
  
  # Reference genes = last numberOfrefGenes
  ref_genes <- tail(genes, numberOfrefGenes)
  
  ref_data <- x[x$Gene %in% ref_genes, ]
  
  # Replicates per condition
  r_ref <- nrow(ref_data) / (length(ref_genes) * length(unique(x$Condition)))
  if (r_ref %% 1 != 0) {
    stop("Unequal replicates among reference genes")
  }
  
  # Average wCt across reference genes per condition & replicate
  ref_mean <- do.call(
    rbind,
    lapply(split(ref_data, ref_data$Condition), function(d) {
      wCt_mat <- matrix(d$wCt, nrow = r_ref, byrow = FALSE)
      data.frame(
        Condition = unique(d$Condition),
        Gene = ref_genes[1],   # pseudo reference gene
        E = NA,
        Ct = NA,
        wCt = rowMeans(wCt_mat)
      )
    })
  )
  
  # Remove original reference genes
  x <- x[!x$Gene %in% ref_genes, ]
  
  # Append averaged reference gene
  x <- rbind(x, ref_mean)
  

  ## DDCt and t-tests
  GENE <- x$Gene
  levels_to_compare <- unique(GENE)[-length(unique(GENE))]
  
  res <- matrix(nrow = length(levels_to_compare), ncol = 7)
  colnames(res) <- c("Gene", "dif", "RE", "LCL", "UCL", "pvalue", "se")
  
  subset <- matrix(NA, nrow = 2 * r_ref, ncol = length(levels_to_compare))
  ttest_result <- vector("list", length(levels_to_compare))
  
  for (i in seq_along(levels_to_compare)) {
    
    subset[, i] <-
      x[GENE == levels_to_compare[i], "wCt"] -
      x[GENE == utils::tail(unique(GENE), 1), "wCt"]
    
    ttest_result[[i]] <-
      stats::t.test(
        subset[(r_ref + 1):(2 * r_ref), i],
        subset[1:r_ref, i],
        paired = paired,
        var.equal = var.equal
      )
    
    res[i, ] <- c(
      levels_to_compare[i],
      mean(subset[(r_ref + 1):(2 * r_ref), i]) - mean(subset[1:r_ref, i]),
      2^-(mean(subset[(r_ref + 1):(2 * r_ref), i]) - mean(subset[1:r_ref, i])),
      round(2^(-ttest_result[[i]]$conf.int[2]), 4),
      round(2^(-ttest_result[[i]]$conf.int[1]), 4),
      ttest_result[[i]]$p.value,
      stats::sd(subset[(r_ref + 1):(2 * r_ref), i]) / sqrt(r_ref)
      )
  }


  
  ## Result table
  res <- as.data.frame(res)
  res$RE <- as.numeric(res$RE)
  res$se <- as.numeric(res$se)
  res$dif <- NULL
  
  res <- data.frame(
    res,
    log2FC = log2(res$RE),
    Lower.se.RE = 2^(log2(res$RE) - res$se),
    Upper.se.RE = 2^(log2(res$RE) + res$se),
    p.adj = stats::p.adjust(res$pvalue, method = p.adj)
  )
  
  res$pvalue <- res$p.adj
  res$p.adj <- NULL
  
  res$sig <- cut(
    res$pvalue,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", " ")
  )
  

  #-------------------- 
  a <- data.frame(res, d = 0, Lower.se.log2FC = 0, Upper.se.log2FC = 0)
  res$Lower.se
  for (i in 1:length(res$RE)) {
    if (res$RE[i] < 1) {
      a$Lower.se.log2FC[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i]
      a$Upper.se.log2FC[i] <- (res$Lower.se.RE[i]*log2(res$RE[i]))/res$RE[i]
      a$d[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i] - 0.2
    } else {
      a$Lower.se.log2FC[i] <- (res$Lower.se.RE[i]*log2(res$RE[i]))/res$RE[i]
      a$Upper.se.log2FC[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i]
      a$d[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i] + 0.2
    }
  }
  
  res <- data.frame(res, 
                    Lower.se.log2FC = a$Lower.se.log2FC,
                    Upper.se.log2FC = a$Upper.se.log2FC)
  #-----------------
  
  
  ## Plotting
  df2 <- res
  
  if (any(order == "none")) {
    df2$Gene <- factor(df2$Gene, levels = default.order)
  } else {
    df2$Gene <- factor(df2$Gene, levels = order)
  }
  
  df2 <- df2[order(df2$Gene), ]
  
  
  df2 <- data.frame(df2, d = 0)
  for (i in 1:length(df2$RE)) {
    if (df2$RE[i] < 1) {
      df2$d[i] <- (df2$Upper.se.RE[i]*log2(df2$RE[i]))/df2$RE[i] - 0.2
    } else {
      df2$d[i] <- (df2$Upper.se.RE[i]*log2(df2$RE[i]))/df2$RE[i] + 0.2
    }
  }
  
  if (plotType == "RE") {
    p <- ggplot(df2, aes(Gene, RE)) +
      geom_col() +
      geom_errorbar(aes(ymin = Lower.se.RE, ymax = Upper.se.RE), width = 0.1) +
      geom_text(aes(label = sig, y = Upper.se.RE + 0.2)) +
      ylab("Relative expression (DDCt method)")
  }
  
  if (plotType == "log2FC") {
    p <- ggplot(df2, aes(Gene, log2FC)) +
      geom_col() +
      geom_errorbar(
        aes(ymin = log2(Lower.se.RE), ymax = log2(Upper.se.RE)),
        width = 0.1
      ) +
      geom_text(aes(label = sig, y = d)) + 
      ylab("log2 fold change")
  }
  
  res <- res %>%
    dplyr::mutate_if(is.numeric, ~ round(., 4))
  
  return(list(Result = res, plot = p))
}