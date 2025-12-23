# Mathematical basis of qPCR data analysis


Real-time polymerase chain reaction (real-time PCR) is widely used in biological studies. Various analysis methods are employed on the real-time PCR data to measure the mRNA levels under different experimental conditions. 
‘rtpcr’ package was developed for amplification efficiency calculation, statistical analysis and bar plot representation of real-time PCR data in R. By accounting for up to two reference genes and amplification efficiency values, a general calculation methodology described by <a href="https://doi.org/10.1186/s12859-017-1949-5">Ganger et al. (2017)</a>
and <a href="https://doi.org/10.1016/j.tibtech.2018.12.002">Taylor et al. (2019)</a>, matching both <a href="https://doi.org/10.1006/meth.2001.1262">Livak and Schmittgen (2001)</a> and <a href="https://doi.org/10.1093/nar/30.9.e36">Pfaffl et al. (2002) </a> methods  was used. Based on the experimental conditions, the functions of the ‘rtpcr’ package use t-test (for experiments with a two-level factor), analysis of variance, analysis of covariance (ANCOVA) or analysis of repeated measure data to calculate the fold change (FC, ${\Delta\Delta C_t}$ method) or relative expression (RE, ${\Delta C_t}$ method). The functions further provide standard errors and confidence interval for means, apply statistical mean comparisons and present significance. To facilitate function application, different data sets were used as examples and the outputs were explained. An outstanding feature of ‘rtpcr’ package is providing publication-ready bar plots with various controlling arguments for experiments with up to three different factors which are further editable by ggplot2 functions.

# Calculation methods
The basic method for expression estimation of a gene between conditions relies on the calculation of fold differences by applying the PCR amplification efficiency (E) and the threshold cycle (syn. crossing point or Ct). Among the various approaches developed for data analysis in real-time PCR, the Livak approach, also known as the $2^{-\Delta\Delta C_t}$ method, stands out for its simplicity and widespread use where the fold change (FC) exoression $(2^{-\Delta\Delta C_t})$ in Treatment (Tr) compared to Control (Co) condition is calculated according to equation:


$$\begin{align*}
\text{Fold change} & = 2^{-\Delta\Delta C_t} \\
& = \frac{2^{-(C_{t_{\text{target}}}-C_{t_{\text{ref}}})_{Tr}}}
{2^{-(C_{t_{\text{target}}}-C_{t_{\text{ref}}})_{Co}}} \\ 
& =2^{-[(C_{t_{\text{target}}}-C_{t_{\text{ref}}})_{\text{Tr}}-
{(C_{t_{\text{target}}}-C_{t_{\text{ref}}})}_{\text{Co}}]} \\ 
& = 2^{-[{(\Delta C_t)_{Tr} - (\Delta C_t)_{Co}}]}
\end{align*}$$


Here, $\Delta C_t$ is the difference between target Ct and reference Ct values for a given sample. Livak method assumes that both the target and reference genes are amplified with efficiencies close to 100%, allowing for the relative quantification of gene expression levels.

On the other hand, the Pfaffl method offers a more flexible approach by accounting for differences in amplification efficiencies between the target and reference genes. This method adjusts the calculated expression ratio by incorporating the specific amplification efficiencies, thus providing a more accurate representation of the relative gene expression levels.

$$\text{Fold change} = \frac{E^{-(C_{t_{\text{Tr}}}-C_{t_{\text{Co}}})_{target}}}
{E^{-(C_{t_{\text{Tr}}}-C_{t_{\text{Co}}})_{ref}}}$$

# A generalized calculation method
The `rtpcr` package was developed for the R environment in the major operating systems. The package functions are mainly based on the calculation of efficiency-weighted $\Delta C_t$ $(w\Delta C_t)$ values from target and reference gene Ct (equation 3). $w\Delta C_t$  values are weighted for the amplification efficiencies as described by Ganger et al. (2017) except that log2 is used instead of log10:


$$w\Delta Ct =\log_{2}(E_{target}).Ct_{target}-\log_{2}(E_{ref}).Ct_{ref}$$


The relative expression of the target gene normalized to that of reference gene(s) within the same sample or condition is called relative expression (RE). From the mean $w\Delta C_t$ values over biological replicates, RE of a target gene can be calculated for each condition according to the equation 

$$\text{Relative Expression} = 2^{-\overline{w\Delta Ct}}$$
Relative expression is only calibrated for the reference gene(s) and not for a control condition. However, often one condition is considered as calibrator and the fold change (FC) expression in other conditions is calculated relative to the calibrator. Examples are Treatment versus Control where Control is served as the calibrator, or time 0 versus time 1 (e.g. after 1 hour) and time 2 (e.g. after 2 hours) where time 0 is served as the reference or calibrator level. So, calibrator is the reference level or sample that all others are compared to. The fold change (FC) expression of a target gene for the reference or calibrator level is 1 because it is not changed compared to itself. The fold change expression of a target gene due to the treatment can be calculated as follows: 

$$\text{Fold Change due to Treatment}=2^{-(\overline{w\Delta Ct}_{\text{Tr}}-{\overline{w\Delta Ct}_{\text{Co}}})}$$

Standard error of the FC and RE means is calculated according to <a href="https://doi.org/10.1016/j.tibtech.2018.12.002">Taylor et al. (2019)</a> in `rtpcr` package.
Here, a brief methodology is presented but details about the $w\Delta C_t$  calculations and statistical analysis are available in <a href="https://doi.org/10.1186/s12859-017-1949-5">Ganger et al. (2017)</a>. Importantly, because relative expression values follow a lognormal distribution, a normal distribution is expected for the $w \Delta C_t$ values making it possible to apply t-test or analysis of variance. Following analysis, $w\Delta C_t$ values are statistically compared and standard deviations and confidence interval are calculated, but the transformation $y = 2^{-x}$ is applied in the final step in order to report the results.


# References
Livak, Kenneth J, and Thomas D Schmittgen. 2001. Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR and the Double Delta CT Method. Methods 25 (4). <a href="https://doi.org/10.1006/meth.2001.1262">doi.org/10.1006/meth.2001.1262</a>.

Ganger, MT, Dietz GD, Ewing SJ. 2017. A common base method for analysis of qPCR data and the application of simple blocking in qPCR experiments. BMC bioinformatics 18, 1-11. <a href="https://doi.org/10.1186/s12859-017-1949-5">doi.org/10.1186/s12859-017-1949-5</a>.

Pfaffl MW, Horgan GW, Dempfle L. 2002. Relative expression software tool (REST©) for group-wise comparison and statistical analysis of relative expression results in real-time PCR. Nucleic acids research 30, e36-e36. <a href="https://doi.org/10.1093/nar/30.9.e36">doi.org/10.1093/nar/30.9.e36</a>.

Taylor SC, Nadeau K, Abbasi M, Lachance C, Nguyen M, Fenrich, J. 2019. The ultimate qPCR experiment: producing publication quality, reproducible data the first time. Trends in Biotechnology, 37(7), 761-774<a href="https://doi.org/10.1016/j.tibtech.2018.12.002">doi.org/10.1016/j.tibtech.2018.12.002</a>.

Yuan, JS, Ann Reed, Feng Chen, and Neal Stewart. 2006. Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85). <a href="https://doi.org/10.1186/1471-2105-7-85">doi.org/10.1186/1471-2105-7-85</a>.