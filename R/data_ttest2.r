#' Sample data (one target, two reference)
#'
#' One target and two reference gens for demonstrating qPCR data analysis.
#'
#' @format A data frame with 18 observations and 4 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{Gene}{Genes}
#'   \item{E}{Amplification efficiency}
#'   \item{Ct}{Ct values}
#' }
#'
#' @source Na
#'
#' @usage data(data_ttest2)
#'
#' @examples
#' data(data_ttest2)
#' data_ttest2
#'
#' @export

data_ttest2 <- read.table(text = "
Condition	Gene	E	Ct 
control	g	2	31.26
control	g	2	31.01
control	g	2	30.97
treatment	g	2	32.65
treatment	g	2	32.03
treatment	g	2	32.4
control	ref1	2	28.86
control	ref1	2	28.42
control	ref1	2	28.56
treatment	ref1	2	28.31
treatment	ref1	2	29.13
treatment	ref1	2	28.62
control	ref2	2	28.43
control	ref2	2	28.9
control	ref2	2	28.17
treatment	ref2	2	28.64
treatment	ref2	2	29.12
treatment	ref2	2	28.95
", header = TRUE, check.names = FALSE)
