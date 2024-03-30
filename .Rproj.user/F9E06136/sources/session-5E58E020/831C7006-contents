#' A sample two factor qPCR dataset
#'
#' A sample dataset for demonstration purposes.
#'
#' @format A data frame with 18 observations and 7 variables:
#' \describe{
#'   \item{Genotype}{First experimental factor}
#'   \item{Drought}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Where the data comes from (if applicable)
#'
#' @usage data(data_2factor)
#'
#' @examples
#' data(data_2factor)
#' data_2factor
#'
#' @export

data_2factor <- read.table(text = "
Genotype    Drought    Rep EPO POCt EGAPDH GAPDHCt
R   0   1   2   33.3    2   31.53
R   0   2   2   33.39   2   31.57
R   0   3   2   33.34   2   31.5
R   0.25    1   2   32.73   2   31.3
R   0.25    2   2   32.46   2   32.55
R   0.25    3   2   32.6    2   31.92
R   0.5 1   2   33.48   2   33.3
R   0.5 2   2   33.27   2   33.37
R   0.5 3   2   33.32   2   33.35
S   0   1   2   26.85   2   26.94
S   0   2   2   28.17   2   27.69
S   0   3   2   27.99   2   27.39
S   0.25    1   2   30.41   2   28.7
S   0.25    2   2   29.49   2   28.66
S   0.25    3   2   29.98   2   28.71
S   0.5 1   2   29.03   2   30.61
S   0.5 2   2   28.73   2   30.2
S   0.5 3   2   28.83   2   30.49", header = TRUE, check.names = FALSE)
