#' @title Delta Delta Ct method wilcox.test analysis
#'
#' @description
#' The \code{WILCOX_DDCt} function performs fold change expression analysis based on
#' the \eqn{\Delta \Delta C_T} method using wilcox.test. It supports analysis
#' of one or more target genes evaluated under two experimental conditions
#' (e.g. control vs treatment).
#'
#' @details
#' Relative expression values are computed using reference gene(s) for normalization.
#' Both paired and unpaired experimental designs are supported.
#'
#' Paired samples in quantitative PCR refer to measurements collected from the same
#' individuals under two different conditions (e.g. before vs after treatment),
#' whereas unpaired samples originate from different individuals in each condition.
#' Paired designs allow within-individual comparisons and typically reduce
#' inter-individual variability.
#'
#' The function returns expression table. The expression table returned by `TTEST_DDCt()`, 
#' `WILCOX_DDCt()`, `ANOVA_DDCt()`, `ANCOVA_DDCt()`, and `REPEATED_DDCt()` functions 
#' include these columns: gene (name of target genes), 
#' contrast (calibrator level and contrasts for which the relative expression is computed), 
#' RE (relative expression or fold change),  log2FC (log(2) of relative expression or fold change), 
#' pvalue, sig (per-gene significance), LCL (95\% lower confidence level), UCL (95\% upper confidence level),
#' se (standard error of mean calculated from the weighted delta Ct values of each of the main factor levels),
#' Lower.se.RE (The lower limit error bar for RE which is 2^(log2(RE) - se)), 
#' Upper.se.RE (The upper limit error bar for RE which is 2^(log2(RE) + se)),
#' Lower.se.log2FC (The lower limit error bar for log2 RE), and 
#' Upper.se.log2FC (The upper limit error bar for log2 RE)
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
#' Logical; if \code{TRUE}, a paired wilcox.test is performed.
#'
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
#' @param set_missing_target_Ct_to_40 If \code{TRUE}, missing target gene Ct 
#'   values become 40; if \code{FALSE} (default), they become NA.
#'
#'
#' @return
#' A table containing RE values, log2FC, p-values, significance, 
#' confidence intervals, standard errors, and lower/upper SE limits.
#'
#' @references
#' Yuan, J. S., Reed, A., Chen, F., and Stewart, N. (2006).
#' Statistical Analysis of Real-Time PCR Data.
#' \emph{BMC Bioinformatics}, 7, 85.
#'
#' @examples
#' # Example data structure
#' data <- read.csv(system.file("extdata", "data_Yuan2006PMCBioinf.csv", package = "rtpcr"))
#'
#' # Unpaired t-test
#' WILCOX_DDCt(
#'   data,
#'   paired = FALSE,
#'   numberOfrefGenes = 1)
#'
#'
#' # Two reference genes
#' data2 <- read.csv(system.file("extdata", "data_1factor_Two_ref.csv", package = "rtpcr"))
#' WILCOX_DDCt(
#'   data2,
#'   numberOfrefGenes = 2,
#'   p.adj = "none")

WILCOX_DDCt <- function(x,
                       numberOfrefGenes,
                       Factor.level.order = NULL,
                       paired = FALSE,
                       p.adj = "none",
                       set_missing_target_Ct_to_40 = FALSE) {
  
  stopifnot(is.data.frame(x))
  stopifnot(numberOfrefGenes >= 1)

  
  # Factor handling
  if (is.null(Factor.level.order)) {
    x[, 1] <- factor(x[, 1], levels = unique(x[, 1]))
    calibrator <- x[, 1][1]
  } else {
    x[, 1] <- factor(x[, 1], levels = Factor.level.order)
    calibrator <- Factor.level.order[1]
  }
  
  on.exit(cat(sprintf(
    "*** The %s level was used as calibrator.\n",
    calibrator
  )))
  
  
  # Column structure
  nc <- ncol(x)
  
  nTargets <- (nc - 2 - numberOfrefGenes * 2) / 2
  if (nTargets < 1) stop("No target genes detected")
  
  cat(sprintf(
    "*** %s target(s) using %s reference gene(s) was analysed!\n",
    nTargets, numberOfrefGenes
  ))
  
  target_start <- 3
  target_end   <- target_start + nTargets * 2 - 1
  ref_start    <- target_end + 1
  
  gene_names <- sub("_E$", "", colnames(x)[seq(target_start, target_end, by = 2)])
  
  # Replicates per condition
  r <- table(x[, 1])[1]
  
  
  ## Result container
  res <- matrix(NA, nrow = nTargets, ncol = 7)
  colnames(res) <- c("gene", "ddCt", "RE", "LCL", "UCL", "pvalue", "se")
  
  
  # Loop over targets
  for (i in seq_len(nTargets)) {
    
    E_col  <- target_start + (i - 1) * 2
    Ct_col <- E_col + 1
    
    # Build per-target wide table
    tmp <- x[, c(1, 2, E_col, Ct_col, ref_start:nc)]
    
    # Compute wDCt using helper
    tmp <- compute_wDCt(
      x = tmp,
      numOfFactors = 1, 
      numberOfrefGenes = numberOfrefGenes, 
      set_missing_target_Ct_to_40 = set_missing_target_Ct_to_40,
      block = NULL
    )
    
    wDCt <- tmp$wDCt
    
    cal_vals <- wDCt[tmp[, 1] == calibrator]
    trt_vals <- wDCt[tmp[, 1] != calibrator]
    

    tt <- suppressWarnings(stats::wilcox.test(
      trt_vals,
      cal_vals,
      paired = paired,
      conf.int = TRUE
    ))
    
    
    set <- stats::t.test(
      trt_vals,
      cal_vals,
      paired = paired, var.equal = F)
    se <- set$stderr
    
    res[i, ] <- c(
      gene_names[i],
      tt$estimate,
      2^(-tt$estimate),
      2^(-tt$conf.int[2]),
      2^(-tt$conf.int[1]),
      tt$p.value,
      se
    )
  }
  
  # Result table
  res <- as.data.frame(res)
  num_cols <- c("ddCt", "RE", "LCL", "UCL", "pvalue", "se")
  res[num_cols] <- lapply(res[num_cols], as.numeric)
  
  # res$estimate <- NULL
  
  res <- transform(
    res,
    log2FC = log2(RE),
    Lower.se.RE = 2^(log2(RE) - se),
    Upper.se.RE = 2^(log2(RE) + se),
    pvalue = stats::p.adjust(pvalue, method = p.adj)
  )
  
  res$sig <- cut(
    res$pvalue,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", " ")
  )
  default.order <- unique(res$Gene)      
  #res$sig[res$RE < 0.001] <- "ND" 
  

  a <- data.frame(res, d = 0, Lower.se.log2FC = 0, Upper.se.log2FC = 0)
  
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
  
  desired_order <- c("gene", "ddCt", "RE", "log2FC", "LCL", "UCL", "se",
                     "Lower.se.RE", "Upper.se.RE", "Lower.se.log2FC", "Upper.se.log2FC", 
                     "pvalue", "sig")
  res <- res[, desired_order]
  

  
  Result <- res %>%
    dplyr::mutate_if(is.numeric, ~ round(., 5))
  
  Result
}
