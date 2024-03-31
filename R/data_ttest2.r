#' A sample data with one target and two reference genes
#'
#' One target and two reference gens for demonstrating qPCR data analysis.
#'
#' @format A data frame with 24 observations and 4 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{E}{Amplification efficiency}
#'   \item{Gene}{Genes}
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
Condition	E	Gene	Ct 
control	2	g	31.26
control	2	g	31.01
control	2	g	30.97
treatment	2	g	32.65
treatment	2	g	32.03
treatment	2	g	32.4
control	2	ref1	28.86
control	2	ref1	28.42
control	2	ref1	28.56
treatment	2	ref1	28.31
treatment	2	ref1	29.13
treatment	2	ref1	28.62
control	2	ref2	28.43
control	2	ref2	28.90
control	2	ref2	28.17
treatment	2	ref2	28.64
treatment	2	ref2	29.12
treatment	2	ref2	28.95
", header = TRUE, check.names = FALSE)
