#' @title Efficiency, standard curves and related statistics
#' @description The \code{efficiency} function calculates amplification efficiency and returns standard curves and related statistics.
#' @details The \code{efficiency} function calculates amplification efficiency of a target and a reference gene, and present the related standard curves along with the Slope, Efficiency, and R2 statistics. The function also compares the slopes of the two standard curves.
#' @author Ghader Mirzaghaderi
#' @export efficiency
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x a data frame
#' @return A list including standard curves along with the Slope, Efficiency, and R2 statistics
#' @examples
#'
#' 
#' # locate and read the sample data
#' df <- system.file('extdata', 'data_efficiency.csv', package = 'rtpcr')
#' df <- read.csv(df)
#'
#' # Applying the efficiency function
#' efficiency(df)
#'
#'



efficiency <- function(x) {

  # renaming the first column
  colnames(x)[1] <- "dilutions"

  # making a long format data
  e <- melt(x, id = "dilutions")

  dilutions <- e$dilutions
  value <- e$value

  colnames(e)[colnames(e) == "variable"] <- "Gene"
  # Fitting list of regression models using lmList function of lme4 package separately for each gene.
  fits <- lmList(value ~ log10(dilutions) | Gene, data = e)

  # Generating a table of Slope, E (efficiency), and R2 statistics for each of the reference and target genes.
  Gene  <- colnames(x)[-1]
  Slope <- numeric(length(Gene))
  E     <- numeric(length(Gene))
  R2    <- numeric(length(Gene))

  for (i in seq_along(Gene)) {
    Slope[i] <- round(as.numeric(fits[[Gene[i]]]$coefficients[2]), digits = 3)
    E[i] <- round(as.numeric(10^(-1/fits[[Gene[i]]]$coefficients[2])), digits = 3)
    R2[i] <- round(summary(fits[[Gene[i]]])$r.squared, digits = 3)
  }

  res <- data.frame(Gene = Gene, Slope = Slope, E = E, R2 = R2)


  # Producing efficiency graphs
  p <- ggplot(e, aes(y = value, x = log10(dilutions))) +
    geom_point() +
    xlab("Log(dilution)") +
    ylab("Ct") +
    theme_bw() +
    geom_smooth(method = "lm", color = "red", fill = "lightyellow") +
    theme(strip.background = element_rect(fill = "#F0FFFF")) +
    geom_text(data = res, aes(x = Inf, y = Inf,
                              label = paste0("Slope = ", Slope, "\n", "E = ", E, "\n", " R2 = ", R2)),
              hjust = 1.0, vjust = 1.0) +
    facet_wrap(~ Gene)

  # Comparing slopes. For this, a regression line is fitted using the DeltaCt values. If 2^-DeltaDelta Ct method is intended, the slope should not exceeds 0.2!
  Dilutions0 <- x$dilutions
  D1 <- data.frame(Dilutions = Dilutions0, diffDilutions = x[,2]-x[,3])
  Dilutions <- D1$Dilutions
  diffDilutions <- D1$diffDilutions
  D <- summarise(group_by(D1, Dilutions), meandd = mean(diffDilutions))
  meandd <- D$meandd
  lm1 <- stats::lm(meandd ~ log10(Dilutions), data = D)
  summary_model <- summary(lm1)
  Slope_of_meand <- summary_model$coefficients[2]

  
  
  outlist1 <- list(plot = p,
                  Efficiency_Analysis_Results = res, 
                  Slope_of_differences = Slope_of_meand)
  
  
  return(outlist1)
  
  return("\nLow values of Slope_of_differences (e.g. less than 0.2) indicates equal slopes for target and reference genes.")
}
