#' Sample data (amplification efficiency)
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
#' @keywords internal
"data_efficiency"

#' Sample data (one factor three levels)
#'
#' A sample dataset for demonstration purposes. Each line belongs to a separate individual (non-repeated measure experiment).
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
#' @source Not applicable
#' @keywords internal
"data_1factor"

#' Sample data (two factor)
#'
#' A sample dataset for demonstration purposes. Each line belongs to a separate individual (non-repeated measure experiment).
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
#' @source Not applicable
#' @keywords internal
"data_2factor"

#' Sample data (two factor with blocking factor)
#'
#' A sample qPCR data set with blocking factor. Each line belongs to a separate individual (non-repeated measure experiment).
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
#' @source Not applicable
#' @keywords internal
"data_2factorBlock"

#' Sample data (three factor)
#'
#' A sample dataset for demonstration purposes. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 36 observations and 8 variables:
#' \describe{
#'   \item{Type}{First experimental factor}
#'   \item{Conc}{Second experimental factor}
#'   \item{SA}{Third experimental factor}
#'   \item{Replicate}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_3factor"


#' Sample qPCR data of an experiment conducted under two different conditions
#'
#' A sample data for demonstrating qPCR data analysis. In unpaired condition, each line belongs to a separate individual.
#' However t-test can be applied for paired samples where the data is acquired from same individuals in two conditions e.g. before and after a treatment. 
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
#' @keywords internal
"data_ttest"

#' Sample data (one target and two reference genes under two different conditions)
#'
#' One target and two reference gens for demonstrating qPCR data analysis. In unpaired condition, each line belongs to a separate individual.
#' However t-test can be applied for paired samples where the data is acquired from same individuals in two conditions e.g. before and after a treatment. 
#'
#' @format A data frame with 18 observations and 4 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{Gene}{Genes}
#'   \item{E}{Amplification efficiency}
#'   \item{Ct}{Ct values}
#' }
#'
#' @source Not applicable
#' @keywords internal
"Taylor_etal2019"

#' Sample data (with technical replicates)
#'
#' A sample data for calculating biological replicated. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 18 observations and 9 variables:
#' \describe{
#'   \item{factor1}{experimental factor}
#'   \item{factor2}{experimental factor}
#'   \item{factor3}{experimental factor}
#'   \item{biolrep}{biological replicate}
#'   \item{techrep}{technical replicates}
#'   \item{Etarget}{Amplification efficiency of target gene}
#'   \item{targetCt}{Ct of target gene}
#'   \item{Eref}{Amplification efficiency of reference gene}
#'   \item{refCt}{Ct of reference gene}
#' } 
#'
#' @source Not applicable
#' @keywords internal
"data_withTechRep"


#' Sample data (with technical replicates)
#'
#' A sample data for calculating biological replicated. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 72 observations and 8 variables:
#' \describe{
#'   \item{factor1}{experimental factor}
#'   \item{DS}{DS}
#'   \item{biolRep}{biological replicate}
#'   \item{techRep}{technical replicates}
#'   \item{APOE_efficiency}{Amplification efficiency of APOE gene}
#'   \item{APOE_Ct}{Ct of APOE gene}
#'   \item{GAPDH_efficiency}{Amplification efficiency of GAPDH gene}
#'   \item{GAPDH_Ct}{Ct of GAPDH gene}
#' }
#'
#' @source Lee et al, (2020) <doi:10.12688/f1000research.23580.2>
#' @keywords internal
"Lee_etal2020qPCR"



#' Repeated measure sample data
#'
#' A repeated measure sample data in which 6 individuals have been analysed.  In the "id" column, a unique number is assigned to each individual,  e.g. all the three number 1 indicate one individual.
#' samples are taken or measurements are scored over different time points (time column) from each individual.
#'
#' @format A data frame with 18 observations and 7 variables:
#' \describe{
#'   \item{id}{experimental factor}
#'   \item{treatment}{treatment}
#'   \item{time}{time course levels}
#'   \item{Eg}{Amplification efficiency of target gene}
#'   \item{Ctg}{Ct of target gene}
#'   \item{Eref}{Amplification efficiency of reference gene}
#'   \item{Ctref}{Ct of reference gene}
#' }
#' 
#' @source NA
#' @keywords internal

"data_repeated_measure_2"


#' Repeated measure sample data
#'
#' A repeated measure sample data in which 3 individuals have been analysed.  In the "id" column, a unique number is assigned to each individual, e.g. all the three number 1 indicate one individual.
#' samples are taken or measurements are scored over different time points (time column) from each individual.
#'
#' @format A data frame with 9 observations and 6 variables:
#' \describe{
#'   \item{id}{experimental factor}
#'   \item{time}{time course levels}
#'   \item{Eg}{Amplification efficiency of target gene}
#'   \item{Ctg}{Ct of target gene}
#'   \item{Eref}{Amplification efficiency of reference gene}
#'   \item{Ctref}{Ct of reference gene}
#' }
#' 
#' @source NA
#' @keywords internal

"data_repeated_measure_1"