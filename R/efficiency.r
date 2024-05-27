#' @title Slope, R2 and Efficiency (E) statistics and standard curves
#' @description The \code{efficiency} function calculates amplification efficiency and returns related statistics and standard curves.
#' @details The \code{efficiency} function calculates amplification efficiency of genes, and present the Slope, Efficiency, and R2 statistics. 
#' @author Ghader Mirzaghaderi
#' @export efficiency
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import purrr
#' @import emmeans
#' @param df a data frame of dilutions and Ct of genes. First column is dilutions and other columns are Ct values for different genes.
#' @return A list 3 elements.
#' \describe{
#'   \item{efficiency}{Slope, R2 and Efficiency (E) statistics}
#'   \item{Slope_compare}{slope comparison table}
#'   \item{plot}{standard curves}
#'   }
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

  

  
  # COMPAIRING SLOPES
  # making a long format data
  e <- melt(df, id = "dilutions")
  
  dilutions <- e$dilutions
  value <- e$value
  
  lm <- lm(value ~ log10(dilutions) * variable, data = e)
  slopes <- emtrends(lm, pairwise ~ variable, var = "log10(dilutions)")
  
  
  
  fits <- lapply(df[,-1], function(x) lm(x ~ log10(df[,1])))
  mdat <- melt(df, id="dilutions")
  Ct <- mdat$value
  Gene <- mdat$variable
  variable <- mdat$Gene
  p <- ggplot(data = mdat) + 
    geom_point(aes(y = Ct, x = log10(dilutions), color = Gene)) +
    geom_smooth(data = mdat,aes(x = log10(dilutions), y = Ct, color = Gene), formula = y ~ x,
                method = "lm", se = F)
  
 
  
  # Compairing intercepts
  # emm = emmeans(lm, "variable", at = list(log10(dilutions) == 0))
  # emmIntercepts <- pairs(emm)


  
  res <- list(Efficiency = results, Slope_compare = slopes, plot = p)
   return(res)
}
