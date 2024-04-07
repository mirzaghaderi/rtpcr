#' Sample data (with technical replicates)
#'
#' A sample data for calculating biological replicated.
#'
#' @format A data frame with 18 observations and 9 variables:
#' \describe{
#'   \item{factor1}{experimental factor}
#'   \item{factor2}{experimental factor}
#'   \item{factor3}{experimental factor}
#'   \item{biolrep}{biological replicate}
#'   \item{techrep}{technical replicates}
#'   \item{Eref}{Amplification efficiency of reference gene}
#'   \item{refCt}{Ct of refence gene}
#'   \item{Etarget}{Amplification efficiency of target gene}
#'   \item{targetCt}{Ct of target gene}
#' }
#'
#' @source Where the data comes from (if applicable)
#'
#' @usage data(data_withTechRep)
#'
#' @examples
#' data(data_withTechRep)
#' data_withTechRep
#'
#' @export

data_withTechRep <- read.table(text = "
factor1	factor2	factor3	biolrep	techrep	Etarget	targetCt	Eref	refCt
Line1	Heat	Ctrl	1	1	2	33.346	2	31.52
Line1	Heat	Ctrl	1	2	2	28.895	2	29.905
Line1	Heat	Ctrl	1	3	2	28.893	2	29.454
Line1	Heat	Ctrl	2	1	2	30.411	2	28.798
Line1	Heat	Ctrl	2	2	2	33.39	2	31.574
Line1	Heat	Ctrl	2	3	2	33.211	2	31.326
Line1	Heat	Ctrl	3	1	2	33.845	2	31.759
Line1	Heat	Ctrl	3	2	2	33.345	2	31.548
Line1	Heat	Ctrl	3	3	2	32.5	2	31.477
Line1	Heat	Treat	1	1	2	33.006	2	31.483
Line1	Heat	Treat	1	2	2	32.588	2	31.902
Line1	Heat	Treat	1	3	2	33.37	2	31.196
Line1	Heat	Treat	2	1	2	36.82	2	31.44
Line1	Heat	Treat	2	2	2	32.75	2	31.3
Line1	Heat	Treat	2	3	2	32.45	2	32.597
Line1	Heat	Treat	3	1	2	35.238	2	31.461
Line1	Heat	Treat	3	2	2	28.532	2	30.651
Line1	Heat	Treat	3	3	2	28.285	2	30.745
", header = TRUE)
