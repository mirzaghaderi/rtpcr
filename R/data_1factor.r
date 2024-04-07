#' Sample data (one factor three levels)
#'
#' A sample dataset for demonstration purposes.
#'
#' @format A data frame with 9 observations and 6 variables:
#' \describe{
#'   \item{SA}{First experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Where the data comes from (if applicable)
#'
#' @usage data(data_1factor)
#'
#' @examples
#' data_1factor
#'
#' @export

data_1factor <- read.table(text = "
SA	Rep	EPO	POCt	EGAPDH	GAPDHCt
L1	1	2	33.3	2	31.53
L1	2	2	33.39	2	31.57
L1	3	2	33.34	2	31.5
L2	1	2	32.73	2	31.3
L2	2	2	32.46	2	32.55
L2	3	2	32.6	2	31.92
L3	1	2	33.48	2	33.3
L3	2	2	33.27	2	33.37
L3	3	2	33.32	2	33.35", header = T)
