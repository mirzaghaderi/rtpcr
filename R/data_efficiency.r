#' A sample dataset of amplification efficiency qPCR data.
#'
#' A sample qPCR dataset for demonstrating efficiency calculation.
#'
#' @format A data frame with 21 observations and 3 variables:
#' \describe{
#'   \item{dilutions}{Dilution factor}
#'   \item{C2H2.26}{Target gene}
#'   \item{GAPDH}{Reference gene}
#' }
#'
#' @source Where the data comes from (if applicable)
#'
#' @usage data(data_efficiency)
#'
#' @examples
#' data(data_efficiency)
#' data_efficiency
#'
#' @export
data_efficiency <- read.table(text = "
dilutions    C2H2.26    GAPDH
1.00 25.57823 22.60794
1.00 25.53636 22.68348
1.00 25.50280 22.62602
0.50 26.70615 23.67162
0.50 26.72720 23.64855
0.50 26.86921 23.70494
0.20 28.16874 25.11064
0.20 28.06759 25.11985
0.20 28.10531 25.10976
0.10 29.19743 26.16919
0.10 29.49406 26.15119
0.10 29.07117 26.15019
0.05 30.16878 27.11533
0.05 30.14193 27.13934
0.05 30.11671 27.16338
0.02 31.34969 28.52016
0.02 31.35254 28.57228
0.02 31.34804 28.53100
0.01 32.55013 29.49048
0.01 32.45329 29.48433
0.01 32.27515 29.26234", header = T, check.names = F)


