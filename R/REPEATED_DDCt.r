#' @title Fold change (\eqn{\Delta\Delta C_T}) analysis of repeated-measure qPCR data
#'
#' @description
#' The \code{REPEATED_DDCt} function performs fold change (FC) analysis using the
#' \eqn{\Delta\Delta C_T} method for qPCR data obtained from repeated measurements
#' over time. Data may originate from uni- or multi-factorial experimental designs.
#'
#' In addition to numerical results, bar plots of relative expression (RE) or log2
#' fold change values with associated uncertainty are optionally produced.
#'
#' @details
#' The analysis is carried out using a linear mixed-effects model in which repeated
#' measurements are accounted for by a random effect of individual (\code{id}).
#' The factor of interest (e.g. time or treatment) is specified via the
#' \code{factor} argument. The first level of this factor (or the level specified
#' by \code{calibratorLevel}) is used as the calibrator.
#'
#' The function supports one or more reference genes. When multiple reference genes
#' are supplied, their contributions are averaged when computing weighted
#' \eqn{\Delta C_T} values.
#'
#' @author Ghader Mirzaghaderi
#'
#' @export
#'
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import emmeans
#' @import lmerTest
#'
#' @param x
#' A data frame in which the first column is the individual identifier (\code{id}),
#' followed by one or more factor columns (including \code{time}).
#' Expression-related columns (time, target gene, reference gene(s)) must appear
#' at the end of the data frame in the required order.
#'
#' @param numberOfrefGenes
#' Integer specifying the number of reference genes (must be \eqn{\ge 1}).
#'
#' @param factor
#' Character string specifying the factor for which fold changes are analysed
#' (commonly \code{"time"}).
#'
#' @param calibratorLevel
#' A level of \code{factor} to be used as the calibrator (reference level).
#'
#' @param block
#' Optional blocking factor column name. If supplied, block effects are treated
#' as random effects.
#'
#' @param x.axis.labels.rename
#' Optional character vector used to replace x-axis labels in the bar plot.
#'
#' @param p.adj
#' Method for p-value adjustment (passed to \code{emmeans}).
#'
#' @param plot
#' Logical; if \code{FALSE}, plots are not produced.
#'
#' @param plotType
#' Either \code{"RE"} (relative expression) or \code{"log2FC"} (log2 fold change).
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{Final_data}{Input data frame augmented with weighted \eqn{\Delta C_T} values.}
#'   \item{lm}{Fitted linear mixed-effects model object.}
#'   \item{ANOVA_table}{ANOVA table for fixed effects.}
#'   \item{Relative_Expression_table}{Table containing RE values, log2FC, p-values,
#'   significance codes, confidence intervals, and standard errors.}
#'   \item{RE_Plot}{Bar plot of relative expression values (if requested).}
#'   \item{log2FC_Plot}{Bar plot of log2 fold change values (if requested).}
#' }
#'
#' @examples
#' REPEATED_DDCt(
#'   data_repeated_measure_1,
#'   numberOfrefGenes = 1,
#'   factor = "time",
#'   calibratorLevel = "1",
#'   block = NULL
#' )
#'
#' REPEATED_DDCt(
#'   data_repeated_measure_2,
#'   numberOfrefGenes = 1,
#'   factor = "time",
#'   calibratorLevel = "1",
#'   block = NULL
#' )



REPEATED_DDCt <- function(x, 
                          numberOfrefGenes,
                          factor, 
                          calibratorLevel,
                          block,
                          x.axis.labels.rename = "none",
                          p.adj = "none",
                          plot = TRUE,
                          plotType = "RE"){
  
  ## ---- basic checks ----
  if (!is.data.frame(x)) stop("`x` must be a data.frame")
  if (missing(factor)) stop("argument 'factor' is missing")
  if (missing(calibratorLevel)) stop("argument 'calibratorLevel' is missing")
  if (!is.numeric(numberOfrefGenes) || numberOfrefGenes < 1)
    stop("`numberOfrefGenes` must be >= 1")
  if (missing(block)) stop("argument 'block' is missing")
  
  # rearrange_repeatedMeasureData
  x <- .rearrange_repeatedMeasureData(x, column_name = factor, level = calibratorLevel)  
  
  
  ## ---- validate number of target genes ----
  ## ---- validate that only ONE target gene exists ----
  expr_cols_expected <- if (is.null(block)) {
    3 + 2 * numberOfrefGenes   # time + target + refs
  } else {
    4 + 2 * numberOfrefGenes   # block + time + target + refs
  }
  
  non_expr_cols <- ncol(x) - expr_cols_expected
  
  if (non_expr_cols < 1) {
    stop(
      "Input data structure error:\n",
      "At least one non-expression column (id) must exist before expression columns.",
      call. = FALSE
    )
  }
  
  ## if expression columns are MORE than expected → extra target genes
  actual_expr_cols <- ncol(x) - non_expr_cols
  
  if (actual_expr_cols != expr_cols_expected) {
    stop(
      sprintf(
        paste0(
          "Exactly ONE target gene is allowed.\n\n",
          "Expected expression columns:\n",
          "  %d  (= time + 1 target + %d reference gene(s)%s)\n\n",
          "But detected:\n",
          "  %d expression-related columns\n\n",
          "This usually means:\n",
          "  more than one target gene is present, or\n",
          "  numberOfrefGenes is incorrect, or\n",
          "  expression columns are not at the end of the data frame."
        ),
        expr_cols_expected,
        numberOfrefGenes,
        if (is.null(block)) "" else " + block",
        actual_expr_cols
      ),
      call. = FALSE
    )
  }
  
  
  
  
  
  
  
  id <- colnames(x)[1]
  
  ## ---- column parsing ----
  if (is.null(block)) {
    
    n_expr <- 3 + 2 * numberOfrefGenes
    factors <- if ((ncol(x) - n_expr) <= 1) NULL else colnames(x)[2:(ncol(x) - n_expr)]
    
    colnames(x)[(ncol(x) - n_expr + 1)] <- "time"
    colnames(x)[(ncol(x) - n_expr + 2)] <- "Etarget"
    colnames(x)[(ncol(x) - n_expr + 3)] <- "Cttarget"
    
    ref_start <- ncol(x) - (2 * numberOfrefGenes) + 1
    ref_cols <- ref_start:ncol(x)
    
  } else {
    
    n_expr <- 4 + 2 * numberOfrefGenes
    factors <- if ((ncol(x) - n_expr) <= 1) NULL else colnames(x)[2:(ncol(x) - n_expr)]
    
    colnames(x)[(ncol(x) - n_expr + 1)] <- "block"
    colnames(x)[(ncol(x) - n_expr + 2)] <- "time"
    colnames(x)[(ncol(x) - n_expr + 3)] <- "Etarget"
    colnames(x)[(ncol(x) - n_expr + 4)] <- "Cttarget"
    
    ref_start <- ncol(x) - (2 * numberOfrefGenes) + 1
    ref_cols <- ref_start:ncol(x)
  }
  
  ## ---- compute wDCt (GENERALIZED) ----
  target_part <- log2(x$Etarget) * x$Cttarget
  
  ref_matrix <- matrix(
    mapply(
      function(E, Ct) log2(E) * Ct,
      x[, ref_cols[seq(1, length(ref_cols), 2)]],
      x[, ref_cols[seq(2, length(ref_cols), 2)]]
    ),
    ncol = numberOfrefGenes
  )
  
  ref_part <- rowMeans(ref_matrix)
  x <- data.frame(x, wDCt = target_part - ref_part)
  
  ## ---- convert factors ----
  for (i in 2:which(names(x) == "time")) {
    x[[i]] <- factor(x[[i]], levels = unique(x[[i]]))
  }
  
  ## ---- model formula ----
  if (is.null(block)) {
    if (is.null(factors)) {
      formula <- wDCt ~ time + (1 | id)
    } else {
      formula <- as.formula(
        paste("wDCt ~ time *", paste(factors, collapse = " * "), "+ (1 | id)")
      )
    }
  } else {
    if (is.null(factors)) {
      formula <- wDCt ~ time + (1 | id) + (1 | block/id)
    } else {
      formula <- as.formula(
        paste("wDCt ~ time *", paste(factors, collapse = " * "),
              "+ (1 | id) + (1 | block/id)")
      )
    }
  }
  
  lm <- lmerTest::lmer(formula, data = x)
  ANOVA <- stats::anova(lm)
  
  #post hoc
  v <- match(colnames(x), factor)
  n <- which(!is.na(v))
  factor <- colnames(x)[n]
  lvls <- unique(x[,n])
  calibrartor <- lvls[1]
  
  on.exit(cat(paste("The level", calibrartor, " of the selected factor was used as calibrator.\n")))
  pp1 <- emmeans(lm, factor, data = x, adjust = p.adj, mode = "satterthwaite")
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  if (length(lvls) >= 3){
    pp3 <- pp2[1:length(lvls) - 1,] 
  } else {
    pp3 <- pp2
  }
  ci <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(lvls)-1,]
  pp <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
  
  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(factor = x[n], bwDCt = bwDCt), x[n]),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))  
  
  
  sig <- .convert_to_character(pp$p.value)
  contrast <- pp$contrast
  post_hoc_test <- data.frame(contrast, 
                              RE = 1/(2^-(pp$estimate)),
                              log2FC = log2(1/(2^-(pp$estimate))),
                              pvalue = pp$p.value,
                              sig = sig,
                              LCL = 1/(2^-pp$lower.CL),
                              UCL = 1/(2^-pp$upper.CL),
                              se = se$se[-1])
  
  
  words <- strsplit(as.character(contrast[1]), " ")[[1]]
  referencelevel <- words[1]
  
  
  reference <- data.frame(contrast = as.character(referencelevel),
                          RE = 1,
                          log2FC = 0,
                          pvalue = 1, 
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])
  
  tableC  <- rbind(reference, post_hoc_test)
  
  #round tableC to 4 decimal places
  #tableC[, sapply(tableC, is.numeric)] <- lapply(tableC[, sapply(tableC, is.numeric)], function(x) round(x, 4))
  
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
  
  
  tableC <- data.frame(tableC, 
                       Lower.se.RE = 2^(log2(tableC$RE) - tableC$se), 
                       Upper.se.RE = 2^(log2(tableC$RE) + tableC$se))  
  ##################################################
  a <- data.frame(tableC, d = 0)
  
  for (i in 1:length(tableC$RE)) {
    if (tableC$RE[i] < 1) {
      a$Lower.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] - 0.2
    } else {
      a$Lower.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] + 0.2
    }
  }
  pfc1 <- ggplot(a, aes(contrast,RE)) + 
    geom_col() +
    geom_errorbar(aes(ymin = tableC$Lower.se.RE, ymax=tableC$Upper.se.RE), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = tableC$Upper.se.RE + 0.2)) +
    ylab("Relative Expression (DDCt)")
  pfc2 <- ggplot(a, aes(contrast,log2FC)) +
    geom_col() +
    geom_errorbar(aes(ymin = Upper.se, ymax=Lower.se), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = d)) +
    ylab("log2FC")
  
  tableC <- data.frame(tableC, Lower.se.log2FC = a$Lower.se, Upper.se.log2FC = a$Upper.se)
  ##################################################    
  tableC <- tableC %>%
    mutate_if(is.numeric, ~ round(., 4))
  
  outlist2 <- structure(list(Final_data = x,
                             lm = lm,
                             ANOVA_table = ANOVA,
                             Relative_Expression_table  = tableC,
                             RE_Plot = pfc1,
                             log2FC_Plot = pfc2), class = "XX")
  
  print.XX <- function(outlist2){
    print(outlist2$ANOVA_table)
    cat("\n", sep = '',"Expression table", "\n")
    print(outlist2$Relative_Expression_table)
    
    if (plot == TRUE){
      if(plotType == "RE"){
        cat("\n", sep = '', "Expression plot", "\n")
        print(outlist2$RE_Plot)
      }else{
        cat("\n", sep = '', "Expression plot", "\n")
        print(outlist2$log2FC_Plot)
      }
    }
    
    invisible(outlist2)
  }
  print.XX(outlist2)
}
