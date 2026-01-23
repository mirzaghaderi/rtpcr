#' Sample qPCR data: amplification efficiency 1
#'
#' A sample qPCR dataset for demonstrating efficiency calculation.
#'
#' @format A data frame with 21 observations and 4 variables:
#' \describe{
#'   \item{Dilutions}{Dilution factor}
#'   \item{C2H2.26}{Target gene 1}
#'   \item{C2H2.01}{Target gene 2}
#'   \item{GAPDH}{Reference gene}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_efficiency1"



			

#' Sample qPCR data: amplification efficiency 2
#'
#' A sample qPCR dataset for demonstrating efficiency calculation.
#'
#' @format A data frame with 12 observations and 5 variables:
#' \describe{
#'   \item{dilutions}{dilution factor}
#'   \item{ref_control}{reference gene in control condition}
#'   \item{ref_treatment}{reference gene in treatment condition}
#'   \item{target_control}{target gene in control condition}
#'   \item{target_treatment}{Reference gene in treatment condition}
#' }
#'
#' @source Yuan2006PMCBioinf
#' @keywords internal
"data_efficiency_Yuan2006PMCBioinf"



#' Sample data (one factor three levels)
#'
#' A sample dataset. Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' @format A data frame with 9 observations and 6 variables:
#' \describe{
#'   \item{SA}{An experimental factor here called SA}
#'   \item{Rep}{Biological replicates}
#'   \item{PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{GAPDH}{Mean amplification efficiency of GAPDH gene}
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
#'   \item{PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{GAPDH}{Mean amplification efficiency of GAPDH gene}
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
#'   \item{Type}{First experimental factor}
#'   \item{Concentration}{Second experimental factor}
#'   \item{block}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{GAPDH}{Mean amplification efficiency of GAPDH gene}
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
#'   \item{PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{GAPDH}{Mean amplification efficiency of GAPDH gene}
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
#'   \item{C2H2_26}{amplification efficiency of C2H2_26 gene}
#'   \item{C2H2_26_Ct}{Ct values of C2H2_26 gene. Each is the mean of technical replicates}
#'   \item{C2H2_01}{Amplification efficiency of C2H2_01 gene}
#'   \item{C2H2_01_Ct}{Ct values of C2H2_01 gene. Each is the mean of technical replicates}
#'   \item{C2H2_12}{Amplification efficiency of C2H2_12 gene}
#'   \item{C2H2_12_Ct}{Ct values of C2H2_12 gene}
#'   \item{ref}{Amplification efficiency of ref gene}
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
#'   \item{target}{Amplification efficiency of target gene}
#'   \item{target_Ct}{Ct values of target gene}
#'   \item{Actin}{Amplification efficiency of reference gene}
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
#'   \item{DER5}{amplification efficiency of DER5 gene}
#'   \item{Ct_DER5}{Ct values of DER5 gene. Each is the mean of technical replicates}
#'   \item{Actin}{Amplification efficiency of Actin gene}
#'   \item{Ct_Actin}{Ct values of Actin gene. Each is the mean of technical replicates}
#'   \item{HPRT}{Amplification efficiency of HPRT gene}
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
#'   \item{Condition}{experimental factor}
#'   \item{biolrep}{biological replicate}
#'   \item{techrep}{technical replicates}
#'   \item{target}{Amplification efficiency of target gene}
#'   \item{Ct_target}{Ct of target gene}
#'   \item{ref}{Amplification efficiency of reference gene}
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
#'   \item{APOE}{Amplification efficiency of APOE gene}
#'   \item{Ct_APOE}{Ct of APOE gene}
#'   \item{GAPDH}{Amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct of GAPDH gene}
#' }
#'
#' @source Lee et al, (2020) <doi:10.12688/f1000research.23580.2>
#' @keywords internal
"data_Lee_etal2020qPCR"


#' Repeated measure sample data
#'
#' A repeated measure sample data in which 6 individuals have been analysed.  In the "id" column, a unique number is assigned to each individual,  e.g. all the three number 1 indicate one individual.
#' samples are taken or measurements are scored over different time points (time column) from each individual.
#'
#' @format A data frame with 18 observations and 7 variables:
#' \describe{
#'   \item{treatment}{treatment}
#'   \item{time}{time course levels}
#'   \item{id}{experimental factor}
#'   \item{Target}{Amplification efficiency of target gene}
#'   \item{Target_Ct}{Ct of target gene}
#'   \item{Ref}{Amplification efficiency of reference gene}
#'   \item{Ref_Ct}{Ct of reference gene}
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
#'   \item{time}{time course levels}
#'   \item{id}{experimental factor}
#'   \item{Target}{Amplification efficiency of target gene}
#'   \item{Ct_Target}{Ct of target gene}
#'   \item{Ref}{Amplification efficiency of reference gene}
#'   \item{Ct_Ref}{Ct of reference gene}
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
#'   \item{PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{GAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{Ct_GAPDH}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#'   \item{ref2}{Mean amplification efficiency of ref2 gene}
#'   \item{Ct_ref2}{Ct values of ref2 gene. Each is the mean of technical replicates}
#'   \item{ref3}{Mean amplification efficiency of ref3 gene}
#'   \item{Ct_ref3}{Ct values of GAPDH gene. Each is the mean of technical replicates}
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
#'   \item{Type}{First experimental factor}
#'   \item{Concentration}{Second experimental factor}
#'   \item{block}{blocking factor}
#'   \item{Rep}{Biological replicates}
#'   \item{PO}{Mean amplification efficiency of PO gene}
#'   \item{Ct_PO}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{NLM}{Mean amplification efficiency of NLM gene}
#'   \item{Ct_NLM}{Ct values of NLM gene. Each is the mean of technical replicates}
#'   \item{ref1}{Mean amplification efficiency of ref1 gene}
#'   \item{Ct_ref1}{Ct values of ref1 gene. Each is the mean of technical replicates}
#'   \item{ref2}{Mean amplification efficiency of ref2 gene}
#'   \item{Ct_ref2}{Ct values of ref2 gene. Each is the mean of technical replicates}
#'   \item{ref3}{Mean amplification efficiency of ref3 gene}
#'   \item{Ct_ref3}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
#' @keywords external
"data_2factorBlock3ref"


#' Sample data in (one factor with one reference gene)
#'
#' A sample qPCR data set with one experimental factor (condition)
#' and one reference gene. Each line belongs to a separate individual
#' (non-repeated measure experiment).
#'
#' @format A data frame with 24 observations and 6 variables:
#' \describe{
#'   \item{condition}{Experimental factor with two levels (control, treatment)}
#'   \item{rep}{Biological replicates}
#'   \item{target}{Mean amplification efficiency of target gene}
#'   \item{Ct_target}{Ct values of target gene. Each is the mean of technical replicates}
#'   \item{ref}{Mean amplification efficiency of reference gene}
#'   \item{Ct_ref}{Ct values of reference gene. Each is the mean of technical replicates}
#' }
#'
#' @source Yuan2006PMCBioinf
#' @keywords internal
"data_Yuan2006PMCBioinf"




#' Sample data in (two-group design with multiple target genes)
#'
#' A sample qPCR data set with two experimental groups (Control and Treatment)
#' and 18 target genes normalized to one reference gene (GAPDH).
#' Each line belongs to a separate individual (non-repeated measure experiment).
#' This dataset is suitable for two-sample t-test based expression analysis.
#'
#' @format A data frame with 8 observations and 38 variables:
#' \describe{
#'   \item{Condition}{Experimental group (Control, Treatment)}
#'   \item{Rep}{Biological replicates}
#'
#'   \item{ANGPT1, Ct.ANGPT1}{Efficiency and Ct values of ANGPT1 gene}
#'   \item{ANGPT2, Ct.ANGPT2}{Efficiency and Ct values of ANGPT2 gene}
#'   \item{CCL2, Ct.CCL2}{Efficiency and Ct values of CCL2 gene}
#'   \item{CCL5, Ct.CCL5}{Efficiency and Ct values of CCL5 gene}
#'   \item{CSF2, Ct.CSF2}{Efficiency and Ct values of CSF2 gene}
#'   \item{FGF2, Ct.FGF2}{Efficiency and Ct values of FGF2 gene}
#'   \item{IL1A, Ct.IL1A}{Efficiency and Ct values of IL1A gene}
#'   \item{IL1B, Ct.IL1B}{Efficiency and Ct values of IL1B gene}
#'   \item{IL6, Ct.IL6}{Efficiency and Ct values of IL6 gene}
#'   \item{IL8, Ct.IL8}{Efficiency and Ct values of IL8 gene}
#'   \item{PDGFA, Ct.PDGFA}{Efficiency and Ct values of PDGFA gene}
#'   \item{PDGFB, Ct.PDGFB}{Efficiency and Ct values of PDGFB gene}
#'   \item{TGFA, Ct.TGFA}{Efficiency and Ct values of TGFA gene}
#'   \item{TGFB, Ct.TGFB}{Efficiency and Ct values of TGFB gene}
#'   \item{TNF, Ct.TNF}{Efficiency and Ct values of TNF gene}
#'   \item{VEGFA, Ct.VEGFA}{Efficiency and Ct values of VEGFA gene}
#'   \item{VEGFB, Ct.VEGFB}{Efficiency and Ct values of VEGFB gene}
#'   \item{VEGFC, Ct.VEGFC}{Efficiency and Ct values of VEGFC gene}
#'   \item{GAPDH, Ct.GAPDH}{Efficiency and Ct values of GAPDH reference gene}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_ttest18genes"



#' Sample data in (repeated-measures design with two reference genes)
#'
#' A sample qPCR data set with repeated measurements over time,
#' a blocking factor, and two reference genes. The experiment includes
#' two treatments (untreated and treated) with repeated measures
#' on the same individuals (id).
#'
#' This dataset is suitable for repeated-measures mixed-model analysis
#' and normalization using multiple reference genes.
#'
#' @format A data frame with 18 observations and 10 variables:
#' \describe{
#'   \item{treatment}{Experimental treatment (untreated, treated)}
#'   \item{time}{Time point of measurement}
#'   \item{blk}{Blocking factor}
#'   \item{id}{Subject/individual identifier for repeated measures}
#'
#'   \item{Target}{Mean amplification efficiency of target gene}
#'   \item{Ct_Target}{Ct values of target gene}
#'
#'   \item{Ref1}{Mean amplification efficiency of reference gene 1}
#'   \item{Ct_Ref1}{Ct values of reference gene 1}
#'
#'   \item{Ref2}{Mean amplification efficiency of reference gene 2}
#'   \item{Ct_Ref2}{Ct values of reference gene 2}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_repeated_measure_3bLock"




#' Sample data in (multi-group design with 6 target genes)
#'
#' A sample qPCR data set with four experimental groups
#' (Uninjected, SNConly, DMSORPE, SNCRPE) and six target genes
#' normalized to one reference gene (Gapdh).
#' Each line belongs to a separate individual (non-repeated measure experiment).
#'
#' This dataset is suitable for t-test or one-way ANOVA based
#' expression analysis of multiple genes.
#'
#' @format A data frame with 18 observations and 16 variables:
#' \describe{
#'   \item{Treatment}{Experimental treatment group}
#'   \item{rep}{Biological replicates}
#'
#'   \item{Fn1, Ct.Fn1}{Efficiency and Ct values of Fn1 gene}
#'   \item{Col1a1, Ct.Col1a1}{Efficiency and Ct values of Col1a1 gene}
#'   \item{Acta2, Ct.Acta2}{Efficiency and Ct values of Acta2 gene}
#'   \item{TgfB, Ct.TgfB}{Efficiency and Ct values of TgfB gene}
#'   \item{Tnfa, Ct.Tnfa}{Efficiency and Ct values of Tnfa gene}
#'   \item{Mcp1, Ct.Mcp1}{Efficiency and Ct values of Mcp1 gene}
#'   \item{Gapdh, Ct.Gapdh}{Efficiency and Ct values of Gapdh reference gene}
#' }
#'
#' @source Not applicable
#' @keywords internal
"data_Heffer2020PlosOne"