#' @title Calculate mean of technical replicates
#'
#' @description
#' Computes the arithmetic mean of technical replicates for each sample or group.
#' This is often performed before ANOVA or other statistical analyses to simplify
#' comparisons between experimental groups.
#'
#' @details
#' The \code{meanTech} function calculates the mean of technical replicates
#' based on one or more grouping columns. This reduces the dataset to a single
#' representative value per group, facilitating downstream analysis such as
#' fold change calculation or ANOVA.
#'
#' @author
#' Ghader Mirzaghaderi
#'
#' @export
#'
#' @import tidyr
#' @import dplyr
#' @import reshape2
#'
#' @param x
#' A raw data frame containing technical replicates.
#'
#' @param groups
#' An integer vector or character vector specifying the column(s) to group
#' by before calculating the mean of technical replicates.
#'
#' @return
#' A data frame with the mean of technical replicates for each group.
#'
#' @examples
#'
#' # Example input data frame with technical replicates
#' data1 <- read.csv(system.file("extdata", "data_withTechRep.csv", package = "rtpcr"))
#'
#' # Calculate mean of technical replicates using first four columns as groups
#' meanTech(data1, groups = 1:4)
#' 
#' # Another example using different dataset and grouping columns
#' data2 <- read.csv(system.file("extdata", "Lee_etal2020qPCR.csv", package = "rtpcr"))
#' meanTech(data2, groups = 1:3)


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
