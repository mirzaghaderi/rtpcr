#' @title Efficiency, standard curves and related statistics
#' @description The \code{efficiency} function calculates amplification efficiency and returns standard curves and related statistics.
#' @details The \code{efficiency} function calculates amplification efficiency of a target and a reference gene, and present the related standard curves along with the Slope, Efficiency, and R2 statistics. The function also compares the slopes of the two standard curves.
#' @author Ghader Mirzaghaderi
#' @export efficiency
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import purrr
#' @param df a data frame
#' @return A list including standard curves along with the Slope, Efficiency, and R2 statistics
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
