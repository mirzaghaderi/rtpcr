#' Sample qPCR data: amplification efficiency
#'
#' A sample qPCR dataset for demonstrating efficiency calculation.
#'
#' @format A data frame with 21 observations and 4 variables:
#' \describe{
#'   \item{dilutions}{Dilution factor}
#'   \item{C2H2.26}{Target gene 1}
#'   \item{C2H2.01}{Target gene 2}
#'   \item{GAPDH}{Reference gene}
#' }
#'
#' @source Where the data comes from (if applicable)
#' @keywords internal
"data_efficiency"


#' Sample data (one factor three levels)
#'
#' A sample dataset. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 9 observations and 6 variables:
#' \describe{
#'   \item{SA}{An experimental factor here called SA}
#'   \item{Rep}{Biological replicates}
#'   \item{E_PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{E_GAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_1factor"

#' Sample data (two factor)
#'
#' A sample dataset. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 18 observations and 7 variables:
#' \describe{
#'   \item{Genotype}{First experimental factor}
#'   \item{Drought}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{E_PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{E_GAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' } 
#'
#' @source Not applicable
#' @keywords internal
"data_2factor"

#' Sample data in (two factor with blocking factor)
#'
#' A sample qPCR data set with blocking factor. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 18 observations and 8 variables:
#' \describe{
#'   \item{factor1}{First experimental factor}
#'   \item{factor2}{Second experimental factor}
#'   \item{block}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{E_PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{E_GAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_2factorBlock"

#' Sample data (three factor)
#'
#' A sample dataset. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 36 observations and 8 variables:
#' \describe{
#'   \item{Type}{First experimental factor}
#'   \item{Conc}{Second experimental factor}
#'   \item{SA}{Third experimental factor}
#'   \item{Replicate}{Biological replicates}
#'   \item{E_PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{E_GAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_3factor"


#' Sample qPCR data (two different conditions)
#'
#' Sample qPCR data (two different conditions)
#'
#' @format A data frame with 6 observations and 10 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{Rep}{Biological replicates}
#'   \item{C2H2_26_E}{amplification efficiency of C2H2_26 gene}
#'   \item{C2H2_26_Ct}{Ct values of C2H2_26 gene. Each is the mean of technical replicates}
#'   \item{C2H2_01_E}{Amplification efficiency of C2H2_01 gene}
#'   \item{C2H2_01_Ct}{Ct values of C2H2_01 gene. Each is the mean of technical replicates}
#'   \item{C2H2_12_E}{Amplification efficiency of C2H2_12 gene}
#'   \item{C2H2_12_Ct}{Ct values of C2H2_12 gene}
#'   \item{ref_E}{Amplification efficiency of ref gene}
#'   \item{ref_Ct}{Ct values of ref gene}
#' }
#'
#' @source University of Kurdistan
#' @keywords internal
"data_1factor_one_ref"


#' Sample qPCR data (two different conditions)
#'
#' Sample qPCR data (two different conditions)
#'
#' @format A data frame with 6 observations and 6 variables:
#' \describe{
#'   \item{Con}{Experimental conditions}
#'   \item{r}{Biological replicates}
#'   \item{target_E}{Amplification efficiency of target gene}
#'   \item{target_Ct}{Ct values of target gene}
#'   \item{Actin_E}{Amplification efficiency of reference gene}
#'   \item{Actin_Ct}{Ct values of reference gene}
#' }
#'
#' @source University of Kurdistan
#' @keywords internal
"data_1factor_one_ref_Eff"

#' Sample qPCR data (one target and two reference genes under two different conditions)
#'
#' Sample qPCR data (one target and two reference genes under two different conditions)
#'
#' @format A data frame with 6 observations and 8 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{Rep}{Biological replicates}
#'   \item{E_DER5}{amplification efficiency of DER5 gene}
#'   \item{Ct_DER5}{Ct values of DER5 gene. Each is the mean of technical replicates}
#'   \item{E_Actin}{Amplification efficiency of Actin gene}
#'   \item{Ct_Actin}{Ct values of Actin gene. Each is the mean of technical replicates}
#'   \item{E_HPRT}{Amplification efficiency of HPRT gene}
#'   \item{Ct_HPRT}{Ct values of HPRT gene}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_1factor_Two_ref"

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
#'   \item{E_target}{Amplification efficiency of target gene}
#'   \item{Ct_target}{Ct of target gene}
#'   \item{E_ref}{Amplification efficiency of reference gene}
#'   \item{Ct_ref}{Ct of reference gene}
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
#'   \item{E_APOE}{Amplification efficiency of APOE gene}
#'   \item{Ct_APOE}{Ct of APOE gene}
#'   \item{E_GAPDH}{Amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct of GAPDH gene}
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
#'   \item{E_target}{Amplification efficiency of target gene}
#'   \item{Ct_target}{Ct of target gene}
#'   \item{E_ref}{Amplification efficiency of reference gene}
#'   \item{Ct_ref}{Ct of reference gene}
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
#'   \item{E_target}{Amplification efficiency of target gene}
#'   \item{Ct_target}{Ct of target gene}
#'   \item{E_ref}{Amplification efficiency of reference gene}
#'   \item{Ct_ref}{Ct of reference gene}
#' }
#' 
#' @source NA
#' @keywords internal
"data_repeated_measure_1"


#' Sample data in (two factor with blocking factor)
#'
#' A sample qPCR data set with blocking factor. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 18 observations and 8 variables:
#' \describe{
#'   \item{factor1}{First experimental factor}
#'   \item{factor2}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{E_PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{E_GAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#'   \item{ref2E}{Mean amplification efficiency of ref2 gene}
#'   \item{ref2Ct}{Ct values of ref2 gene. Each is the mean of technical replicates}
#'   \item{ref3E}{Mean amplification efficiency of ref3 gene}
#'   \item{ref3Ct}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_2factor3ref"


#' Sample data in (two factor with blocking factor and 3 reference genes)
#'
#' A sample qPCR data set with blocking factor and 3 reference genes. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 18 observations and 8 variables:
#' \describe{
#'   \item{factor1}{First experimental factor}
#'   \item{factor2}{Second experimental factor}
#'   \item{block}{blocking factor}
#'   \item{Rep}{Biological replicates}
#'   \item{E_PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{E_GAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#'   \item{ref2E}{Mean amplification efficiency of ref2 gene}
#'   \item{ref2Ct}{Ct values of ref2 gene. Each is the mean of technical replicates}
#'   \item{ref3E}{Mean amplification efficiency of ref3 gene}
#'   \item{ref3Ct}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_2factorBlock3ref"
