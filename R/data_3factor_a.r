#' A sample multi-factorial qPCR dataset
#'
#' A sample dataset for demonstration purposes.
#'
#' @format A data frame with 36 observations and 8 variables:
#' \describe{
#'   \item{Genotype}{First experimental factor}
#'   \item{Drought}{Second experimental factor}
#'   \item{SA}{Third experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Where the data comes from (if applicable)
#'
#' @usage data(data_3factor_a)
#'
#' @examples
#' data(data_3factor_a)
#' data_3factor_a
#'
#' @export
data_3factor_a <- read.table(text = "
Genotype Drought SA Rep EPO POCt EGAPDH GAPDHCt
R	0	A1	1	1.839	33.3	1.918	31.53
R	0	A1	2	1.839	33.39	1.918	31.57
R	0	A1	3	1.839	33.34	1.918	31.5
R	0	A2	1	1.839	34.01	1.918	31.48
R	0	A2	2	1.839	36.82	1.918	31.44
R	0	A2	3	1.839	35.44	1.918	31.46
R	0.25	A1	1	1.839	32.73	1.918	31.3
R	0.25	A1	2	1.839	32.46	1.918	32.55
R	0.25	A1	3	1.839	32.6	1.918	31.92
R	0.25	A2	1	1.839	33.37	1.918	31.19
R	0.25	A2	2	1.839	33.12	1.918	31.94
R	0.25	A2	3	1.839	33.21	1.918	31.57
R	0.5	A1	1	1.839	33.48	1.918	33.3
R	0.5	A1	2	1.839	33.27	1.918	33.37
R	0.5	A1	3	1.839	33.32	1.918	33.35
R	0.5	A2	1	1.839	32.53	1.918	33.47
R	0.5	A2	2	1.839	32.61	1.918	33.26
R	0.5	A2	3	1.839	32.56	1.918	33.36
S	0	A1	1	1.839	26.85	1.918	26.94
S	0	A1	2	1.839	28.17	1.918	27.69
S	0	A1	3	1.839	27.99	1.918	27.39
S	0	A2	1	1.839	28.71	1.918	29.45
S	0	A2	2	1.839	29.01	1.918	29.46
S	0	A2	3	1.839	28.82	1.918	29.48
S	0.25	A1	1	1.839	30.41	1.918	28.7
S	0.25	A1	2	1.839	29.49	1.918	28.66
S	0.25	A1	3	1.839	29.98	1.918	28.71
S	0.25	A2	1	1.839	28.91	1.918	28.09
S	0.25	A2	2	1.839	28.6	1.918	28.65
S	0.25	A2	3	1.839	28.59	1.918	28.37
S	0.5	A1	1	1.839	29.03	1.918	30.61
S	0.5	A1	2	1.839	28.73	1.918	30.2
S	0.5	A1	3	1.839	28.83	1.918	30.49
S	0.5	A2	1	1.839	28.29	1.918	30.84
S	0.5	A2	2	1.839	28.53	1.918	30.65
S	0.5	A2	3	1.839	28.28	1.918	30.74
", header = T, check.names = F)


