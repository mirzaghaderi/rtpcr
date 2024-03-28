#' A sample two conditional qPCR data
#'
#' A sample data for demonstrating qPCR data analysis.
#'
#' @format A data frame with 24 observations and 4 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{E}{Amplification efficiency}
#'   \item{Gene}{Genes}
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
Condition	E	Gene	Ct
control	2	C2H2-26	31.26
control	2	C2H2-26	31.01
control	2	C2H2-26	30.97
treatment	2	C2H2-26	32.65
treatment	2	C2H2-26	32.03
treatment	2	C2H2-26	32.4
control	2	C2H2-01	31.06
control	2	C2H2-01	30.41
control	2	C2H2-01	30.97
treatment	2	C2H2-01	28.85
treatment	2	C2H2-01	28.93
treatment	2	C2H2-01	28.9
control	2	C2H2-12	28.5
control	2	C2H2-12	28.4
control	2	C2H2-12	28.8
treatment	2	C2H2-12	27.9
treatment	2	C2H2-12	28
treatment	2	C2H2-12	27.9
control	2	ref	28.87
control	2	ref	28.42
control	2	ref	28.53
treatment	2	ref	28.31
treatment	2	ref	29.14
treatment	2	ref	28.63
", header = T)
