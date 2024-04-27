#' @title Slope, R2 and Efficiency (E) statistics
#' @description The \code{efficiency} function calculates amplification efficiency and returns related statistics.
#' @details The \code{efficiency} function calculates amplification efficiency of genes, and present the Slope, Efficiency, and R2 statistics. 
#' @author Ghader Mirzaghaderi
#' @export efficiency
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
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

  

  
  # COMPAIRING SLOPES
  # making a long format data
  e <- melt(df, id = "dilutions")
  
  dilutions <- e$dilutions
  value <- e$value
  
  lm <- lm(value ~ log10(dilutions) * variable, data = e)
  slopes <- emtrends(lm, pairwise ~ variable, var = "log10(dilutions)")
  
  
  
  
  
  fits <- lapply(df[,-1], function(x) lm(x ~ log10(df[,1])))
  mdat <- melt(df,id="dilutions")
  variable <- mdat$variable
  p <- ggplot(data = mdat) + 
    geom_point(aes(y = value, x = log10(dilutions), color=variable)) +
    geom_smooth(data=mdat,aes(x = log10(dilutions),y=value,color = variable),formula = y ~ x,
                method = "lm", se = F)
  
  
  
  res <- list(Efficiency = results, Slope_compare = slopes, plot = p)
   return(res)
}
