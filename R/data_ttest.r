#' Sample data (one factor-two level qPCR)
#'
#' A sample data for demonstrating qPCR data analysis.
#'
#' @format A data frame with 24 observations and 4 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{Gene}{Genes}
#'   \item{E}{Amplification efficiency}
#'   \item{Ct}{Ct values}
#' }
#'
#' @source University of Kurdistan
#'
#' @usage data(data_ttest)
#'
#' @examples
#' data(data_ttest)
#' data_ttest
#'
#' @export
data_ttest <- read.table(text = "
Condition	Gene	E	Ct
control	C2H2-26	2	31.26
control	C2H2-26	2	31.01
control	C2H2-26	2	30.97
treatment	C2H2-26	2	32.65
treatment	C2H2-26	2	32.03
treatment	C2H2-26	2	32.4
control	C2H2-01	2	31.06
control	C2H2-01	2	30.41
control	C2H2-01	2	30.97
treatment	C2H2-01	2	28.85
treatment	C2H2-01	2	28.93
treatment	C2H2-01	2	28.9
control	C2H2-12	2	28.5
control	C2H2-12	2	28.4
control	C2H2-12	2	28.8
treatment	C2H2-12	2	27.9
treatment	C2H2-12	2	28
treatment	C2H2-12	2	27.9
control	ref	2	28.87
control	ref	2	28.42
control	ref	2	28.53
treatment	ref	2	28.31
treatment	ref	2	29.14
treatment	ref	2	28.63
", header = T)
