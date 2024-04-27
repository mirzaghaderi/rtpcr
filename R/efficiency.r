#' @title Slope, R2 and Efficiency (E) statistics
#' @description The \code{efficiency} function calculates amplification efficiency and returns related statistics.
#' @details The \code{efficiency} function calculates amplification efficiency of genes, and present the Slope, Efficiency, and R2 statistics. 
#' @author Ghader Mirzaghaderi
#' @export efficiency
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import purrr
#' @param df a data frame of dilutions and Ct of genes. First column is dilutions and other columns are Ct values for different genes.
#' @return A data frame including  Slope, R2 and Efficiency (E) statistics for each gene.
#' @examples
#'
#' 
#' # locate and read the sample data
#' data_efficiency
#'
#' # Applying the efficiency function
#' efficiency(data_efficiency)
#'
#'



efficiency <- function(df) {

  
  # renaming the first column
  colnames(df)[1] <- "dilutions"
  dilutions <- df$dilutions
  # Fit the linear regressions and extract the slope and R-squared
  results <- df %>%
    select(-dilutions) %>%
    map_df(~{
      model <- lm(. ~ log10(dilutions))
      Slope <- coef(model)[2]
      R2 <- summary(model)$r.squared
      E <- 10^(-1/coef(model)[2])
      data.frame(Slope, R2, E)
    })
  
  # Add the column names to the results data frame
  results <- cbind(gene = colnames(df)[2:ncol(df)], results)
  colnames(results) <- c("Gene", "Slope", "R2", "E")
  rownames(results) <- NULL

   return(results)
}
