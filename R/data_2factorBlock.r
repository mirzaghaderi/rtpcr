#' Sample data (two factor with blocking factor)
#'
#' A sample qPCR data set with blocking factor.
#'
#' @format A data frame with 18 observations and 8 variables:
#' \describe{
#'   \item{factor1}{First experimental factor}
#'   \item{factor2}{Second experimental factor}
#'   \item{block}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Where the data comes from (if applicable)
#'
#' @usage data(data_2factorBlock)
#'
#' @examples
#' data(data_2factorBlock)
#' data_2factorBlock
#'
#' @export

data_2factorBlock <- read.table(text = "
factor1	factor2	block	Rep	EPO	POCt	EGAPDH	GAPDHCt
R	0	1	1	2	33.3	2	31.53
R	0	1	2	2	33.39	2	31.57
R	0	2	3	2	33.34	2	31.5
R	0.25	1	1	2	32.73	2	31.3
R	0.25	1	2	2	32.46	2	32.55
R	0.25	2	3	2	32.6	2	31.92
R	0.5	1	1	2	33.48	2	33.3
R	0.5	1	2	2	33.27	2	33.37
R	0.5	2	3	2	33.32	2	33.35
S	0	1	1	2	26.85	2	26.94
S	0	1	2	2	28.17	2	27.69
S	0	2	3	2	27.99	2	27.39
S	0.25	1	1	2	30.41	2	28.7
S	0.25	1	2	2	29.49	2	28.66
S	0.25	2	3	2	29.98	2	28.71
S	0.5	1	1	2	29.03	2	30.61
S	0.5	1	2	2	28.73	2	30.2
S	0.5	2	3	2	28.83	2	30.49
", header = TRUE, check.names = FALSE)
