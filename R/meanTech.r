#' @title Calculating mean of technical replicates
#' @description Calculating arithmetic mean of technical replicates for subsequent ANOVA analysis
#' @details The meanTech calculates mean of technical replicates. Arithmetic mean of technical replicates can be calculated in order to 
#' simplify the statistical comparison between sample groups.
#' @author Ghader Mirzaghaderi
#' @export meanTech
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @param x  A raw data frame including technical replicates.
#' @param groups grouping columns based on which the mean technical replicates are calculated.
#' @return A data frame with the mean of technical replicates.
#'
#' @examples
#'
#' # See example input data frame:
#' data_withTechRep
#'
#' # Calculating mean of technical replicates
#' meanTech(data_withTechRep, groups = 1:4)
#' 
#' # Calculating mean of technical replicates
#' meanTech(Lee_etal2020qPCR, groups = 1:3)
#'
#'

meanTech <- function(x, groups){

  g <- colnames(x)[groups]

  # Concatenate the columns using paste0
  x$T <- do.call(paste, c(x[groups], sep = ":"))

  data_meanBiolRep <- x %>%
    group_by(across(all_of(g))) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))

  y <- data_meanBiolRep
  m <- ncol(y)-4
  databiol <- y[,-m]
  return(as.data.frame(databiol))
}
