---
title: "rtpcr: a package for statistical analysis of real-time PCR data in R"
author: Ghader Mirzaghaderi
output:
  html_document: 
    toc: yes
    keep_md: yes
    output: rmarkdown::html_vignette
    df_print: default
  pdf_document: 
    toc: yes
    latex_engine: lualatex
  word_document:
    toc: No
vignette: |
  %\VignetteIndexEntry{Sending Messages With Gmailr}
  %\usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::knitr}
editor_options:
  markdown:
    wrap: 90
  chunk_output_type: console
---




# Overview


Real-time polymerase chain reaction (real-time PCR), is widely used in biological research. Various analysis methods are employed on the real-time PCR data to measure the mRNA levels under different experimental conditions. 
‘rtpcr’ package was developed for amplification efficiency calculation and statistical analysis of real-time PCR data in R. By accounting for up to two reference genes and amplification efficiency values, a general calculation methodology described by <a href="https://doi.org/10.1186/s12859-017-1949-5">Ganger et al. (2017)</a>, matching both <a href="https://doi.org/10.1006/meth.2001.1262">Livak and Schmittgen (2001)</a> and <a href="https://doi.org/10.1093/nar/30.9.e36">Pfaffl et al. (2002) </a> methods  was used. Based on the experimental conditions, the functions of the ‘rtpcr’ package use a t-test (for experiments with a two-level factor) or analysis of variance (for cases where more than two levels or factors or a blocking factor exist) to calculate the fold change (FC) or relative expression (RE). The functions further provide standard deviations and confidence limits for means, apply statistical mean comparisons and present letter mean grouping. To facilitate function application, different data sets were used as examples and the outputs were explained. An outstanding feature of ‘rtpcr’ package is providing publication-ready bar plots with various controlling arguments for experiments with up to three different factors.

# Calculation methods
The basic method for expression estimation of a gene between conditions relies on the calculation of fold differences by means of the PCR amplification efficiency (E) and the threshold cycle (syn. crossing point or Ct). Among the various approaches developed for data analysis in real-time PCR, the Livak approach, also known as the $2^{-\Delta\Delta C_t}$ method, stands out for its simplicity and widespread use where the fold change (FC) exoression ($2^{-\Delta\Delta C_t}$) in Treatment (Tr) compared to Control (Co) condition is calculated according to equation:


$$\begin{align*}
\text{Fold change} & = 2^{-\Delta\Delta C_t} \\
& = \frac{2^{-(C_{t_{\text{target}}}-C_{t_{\text{ref}}})_{Tr}}}
{2^{-(C_{t_{\text{target}}}-C_{t_{\text{ref}}})_{Co}}} \\ 
& =2^{-[(C_{t_{\text{target}}}-C_{t_{\text{ref}}})_{\text{Tr}}-
{(C_{t_{\text{target}}}-C_{t_{\text{ref}}})}_{\text{Co}}]} \\ 
& = 2^{-[{(\Delta C_t)_{Tr} - (\Delta C_t)_{Co}}]}
\end{align*}$$


Here, $\Delta C_t$ is the difference between two Ct values (e.g. Cttarget−Ctref) and target and ref are target gene and reference genes, respectively. This method assumes that both the target and reference genes are amplified with efficiencies close to 100%, allowing for the relative quantification of gene expression levels 
This method assumes that both the target and reference genes are amplified with efficiencies close to 100%, allowing for the relative quantification of gene expression levels.

On the other hand, the Pfaffl method offers a more flexible approach by accounting for differences in amplification efficiencies between the target and reference genes. This method adjusts the calculated expression ratio by incorporating the specific amplification efficiencies, thus providing a more accurate representation of the relative gene expression levels.

$$\text{Fold change} = \frac{E^{-(C_{t_{\text{Tr}}}-C_{t_{\text{Co}}})_{target}}}
{E^{-(C_{t_{\text{Tr}}}-C_{t_{\text{Co}}})_{ref}}}$$

# A generalized calculation method
The rtpcr package was developed for the R environment in the major operating systems. The packager functions are mainly based on the calculation of efficiency-weighted $\Delta C_t$ ($w\Delta C_t$) values from target and reference gene Ct (equation 3). $w\Delta C_t$  values are weighted for the amplification efficiencies as described by Ganger et al. (2017):


$$w\Delta Ct =\log_{10}(E_{target}).Ct_{target}-\log_{10}(E_{ref}).Ct_{ref}$$

From the mean wΔCt values over biological replicates, relative expression (RE) of a target gene can be calculated for each condition according to the equation 

$$\text{Relative Expression} = 10^{-\overline{w\Delta Ct}}$$
Often, one condition is considered as calibrator condition. Examples are Treatment versus Control where Control is served as the calibrator, or time 0 versus time 1 (e.g. after 1 hour) and time 2 (e.g. after 2 hours) where time 0 is served as the reference or calibrator level. So, calibrator is the reference level or sample that all others are compared to. The fold change (FC) expression of a target gene for the reference or calibrator level is 1 because it is not changed compared to itself. The fold change expression of a target gene due to the treatment can be calculated as follows: 

$$\text{Fold Change due to Treatment}=10^{-(\overline{w\Delta Ct}_{\text{Tr}}-{\overline{w\Delta Ct}_{\text{Co}}})}$$

`qpcrTTEST` and `qpcrTTESTplot` functions calculate FC for multi-genes-two conditional cases and `qpcrANCOVA` represents FC for single-gene-single factor or single-gene-multi factorial experiments. If wDCt values is calculated from the E valuses, these calculations match the formula of Pfaffl while if 2 (complete efficiency) be used instead, the result match the $2^{-\Delta\Delta C_t}$ method. In any case we called these as Fold Change in the outputs of `rtpcr`. Under factorial experiments where the calculation of the expression of the target gene relative to the reference gene (called Relative Expression) in each condition is desired, `qpcrANOVA`, `oneFACTORplot`, `twoFACTORplot` and `threeFACTORplot` functions were developed for ANOVA analysis, and representing the plots from single, double or triple factor experiments, respectively. The last three functions generate ‘ggplot2’-derived graphs based on the output of the `qpcrANOVA` function. If available, the blocking factor can alaso be handled by `qpcrTTEST`, `qpcrANOVA` and `qpcrANCOVA` functions.
Here, a brief methodology is presented but details about the $w\Delta C_t$  calculations and statistical analysis are available in. Importantly, because the relative or FC gene expression follow a lognormal distribution, a normal distribution is expected for the $w\Delta C_t$ or $w\Delta\Delta C_t$ values making it possible to apply t-test or analysis of variance to them. Following analysis, $w\Delta C_t$ values are statistically compared and standard deviations and confidence intervals is calculated, but the transformation $y = 10^x$ is applied in the final step in order to report the results (i.e. RE ratios, errors and confidence limits).





# Installing and loading



The rtpcr package can be installed and loaded using:
```r
install.packages("rtpcr")
library(rtpcr)
```

Alternatively, the `rtpcr` with the latest changes can be installed by running the following code in your R software: 


```r
# install `rtpcr` from github (under development)

devtools::install_github("mirzaghaderi/rtpcr")

# I strongly recommend to install the package with the vignette as it contains information about how to use the 'rtpcr' package. Through the following code, Vignette is installed as well.

devtools::install_github("mirzaghaderi/rtpcr", build_vignettes = TRUE)
```



# Data structure and column arrangement

To use the functions, input data should be prepared in the right format with appropriate column arrangement. The correct column arrangement is shown in Table 1.

*Table 1. Data structure and column arrangement required for ‘rtpcr’ package.  rep: technical replicate; targetE and refE: amplification efficiency columns for target and reference genes respectively. targetCt and refCt: target and reference Ct columns, respectively. factors (factor1, factor2 and/or factor3): experimental factors.*

 | Experiment type     |  Column arrangement of the input data   | Example in the package                   |
 |:---------------------|:---------------------------------------|:------------------------------------------|
 |Amplification efficiency             |Dilutions - targetCt - refCt | data_efficiency |
 |t-test (accepts multiple genes)      |condition (put the control level first) - gene (put reference gene(s) last.)- efficiency - Ct  | data_ttest |
 |Factorial (Up to three factors)      |factor1 - rep - targetE - targetCt - refE - refCt | data_1factor |
 |                                     |factor1 - factor2 - rep - targetE - targetCt - refE - refCt | data_2factor |
 |                                     |factor1 - factor2 - factor3 - rep - targetE - targetCt - refE - refCt | data_3factor |
 |Factorial with blocking              |factor1 - block - rep - targetE - targetCt - refE - refCt | |
 |                                     |factor1 - factor2 - block - rep - targetE - targetCt - refE - refCt	 | data_2factorBlock |
 |                                     |factor1 - factor2 - factor3 - block - rep - targetE - targetCt - refE - refCt | |
 |Two reference genes                  |. . . . . .  rep - targetE - targetCt - ref1E - ref1Ct - ref2E - ref2Ct | |
 |calculating biological replicated    |. . . . . .  biologicalRep - techcicalRep - Etarget - targetCt - Eref - refCt  | data_withTechRep |
 |                |. . . . . .  biologicalRep - techcicalRep - Etarget - targetCt - ref1E - ref1Ct - ref2E - ref2Ct  | |
 
 
Example data sets can be presented by runing the name of each data set. A description of the columns names in  each data set is called by "?" followed by the names of the data set.

# functions usage

To simplify 'rtpcr' usage, examples for using the functions are presented below.


           
*Table 2. Functions and examples for using them.*

 | function             |                 Analysis | Example (see package help for the more arguments) |
 |:---------------------|:-----------------------------------|:----------------------------------|
 | efficiency             | Efficiency, standard curves and related statistics | efficiency(data_efficiency) |
 | meanTech      | Calculating the mean of technical replicates | meanTech(data_withTechRep, groups = 1:4) |
 | qpcrANCOVA      | FC and bar plot of the target gene (one or multi-factorial experiments) | qpcrANCOVA(data_1factor, numberOfrefGenes = 1, mainFactor.column = 1, mainFactor.level.order = c("L1", "L2", "L3") |
 |  oneFACTORplot    | Bar plot of the relative gene expression from a one-factor experiment | out <- qpcrANOVA(data_1factor, numberOfrefGenes = 1);   oneFACTORplot(out) |
 |  qpcrANOVA  | Analysis of Variance of the qpcr data  | qpcrANOVA(data_3factor, numberOfrefGenes = 1)|
 | qpcrTTEST     | Computing the average fold change and related statistics | qpcrTTEST(data_ttest,  numberOfrefGenes = 1, paired = FALSE, var.equal = TRUE) |
 | qpcrTTESTplot  | Bar plot of the average fold change of the target genes	 | qpcrTTESTplot(data_ttest,  numberOfrefGenes = 1, order = c("C2H2-01", "C2H2-12", "C2H2-26")) |
 |  threeFACTORplot  | Bar plot of the relative gene expression from a three-factor experiment | res <- qpcrANOVA(data_3factor,  numberOfrefGenes = 1); threeFACTORplot(res, arrangement = c(3, 1, 2)) |
 | twoFACTORplot   | Bar plot of the relative gene expression from a two-factor experiment | res <- qpcrANOVA(data_2factor, numberOfrefGenes = 1); twoFACTORplot(res, x.axis.factor = Genotype, group.factor = Drought) |
 

 
*see package help for more arguments including the number of reference genes, levels arrangement, blocking, and arguments for adjusting the bar plots.*
 

# Amplification efficiency data analysis
## Sample data of amplification efficiency

To calculate the amplification efficiencies of a target and a reference gene, a data frame should be prepared with 3 columns of dilutions, target gene Ct values, and reference gene Ct values, respectively,  as shown below.


```r
data_efficiency
```

```
##    dilutions  C2H2.26    GAPDH
## 1       1.00 25.57823 22.60794
## 2       1.00 25.53636 22.68348
## 3       1.00 25.50280 22.62602
## 4       0.50 26.70615 23.67162
## 5       0.50 26.72720 23.64855
## 6       0.50 26.86921 23.70494
## 7       0.20 28.16874 25.11064
## 8       0.20 28.06759 25.11985
## 9       0.20 28.10531 25.10976
## 10      0.10 29.19743 26.16919
## 11      0.10 29.49406 26.15119
## 12      0.10 29.07117 26.15019
## 13      0.05 30.16878 27.11533
## 14      0.05 30.14193 27.13934
## 15      0.05 30.11671 27.16338
## 16      0.02 31.34969 28.52016
## 17      0.02 31.35254 28.57228
## 18      0.02 31.34804 28.53100
## 19      0.01 32.55013 29.49048
## 20      0.01 32.45329 29.48433
## 21      0.01 32.27515 29.26234
```


	



## Calculating amplification efficiency
The following `efficiency` function calculates the amplification efficiency of a target and a reference gene and presents the related standard curves along with the Slope, Efficiency, and R2 statistics. The function also compares the slopes of the two standard curves. For this, a regression line is fitted using the $\Delta C_t$ values. If $2^{-\Delta\Delta C_t}$ method is intended, the slope should not exceed 0.2!


```r
efficiency(data_efficiency)
```

```
## $plot
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-4-1.png" alt="Standard curve and the amplification efficiency analysis of target and reference genes. A sample data arrangement that is required as input for the calculation of amplification efficiency by the efficiency function."  />
<p class="caption">Standard curve and the amplification efficiency analysis of target and reference genes. A sample data arrangement that is required as input for the calculation of amplification efficiency by the efficiency function.</p>
</div>

```
## 
## $Efficiency_Analysis_Results
##      Gene  Slope     E    R2
## 1 C2H2.26 -3.388 1.973 0.997
## 2   GAPDH -3.415 1.963 0.999
## 
## $Slope_of_differences
## [1] 0.0264574
```

# Expression data analysis

## Target genes in two conditions (t-test)

### Example data
When a target gene is assessed under two different conditions (for example Control and treatment), it is possible to calculate the average fold change expression i.e. $2^{-\Delta \Delta C_t}$ of the target gene in treatment relative to control conditions. For this, the data should be prepared according to the following data set consisting of 4 columns belonging to condition levels, E (efficiency), genes and Ct values, respectively. Each Ct value is the mean of technical replicates. Complete amplification efficiencies of 2 have been assumed here for all wells but the calculated efficiencies can be used instead. 


```r
data_ttest
```

```
##    Condition    Gene E    Ct
## 1    control C2H2-26 2 31.26
## 2    control C2H2-26 2 31.01
## 3    control C2H2-26 2 30.97
## 4  treatment C2H2-26 2 32.65
## 5  treatment C2H2-26 2 32.03
## 6  treatment C2H2-26 2 32.40
## 7    control C2H2-01 2 31.06
## 8    control C2H2-01 2 30.41
## 9    control C2H2-01 2 30.97
## 10 treatment C2H2-01 2 28.85
## 11 treatment C2H2-01 2 28.93
## 12 treatment C2H2-01 2 28.90
## 13   control C2H2-12 2 28.50
## 14   control C2H2-12 2 28.40
## 15   control C2H2-12 2 28.80
## 16 treatment C2H2-12 2 27.90
## 17 treatment C2H2-12 2 28.00
## 18 treatment C2H2-12 2 27.90
## 19   control     ref 2 28.87
## 20   control     ref 2 28.42
## 21   control     ref 2 28.53
## 22 treatment     ref 2 28.31
## 23 treatment     ref 2 29.14
## 24 treatment     ref 2 28.63
```

### Data analysis under two conditions

Here, the above data set was used for the Fold Change expression analysis of the target genes using the `qpcrTTEST` function. This function performs a t-test-based analysis of any number of genes that 
have been evaluated under control and treatment conditions. The analysis can be done for unpaired or paired conditions. The output is a table of target gene names, fold changes confidence limits, and the t.test derived p-values. The `qpcrTTEST` function includes the `var.equal` argument. When set to `FALSE`, `t.test` is performed under the unequal variances hypothesis.


```r
qpcrTTEST(data_ttest, 
          numberOfrefGenes = 1,
          paired = F, 
          var.equal = T)
```

```
## $Raw_data
##    Var2       wDCt
## 1     1  0.7194617
## 2     1  0.7796677
## 3     1  0.7345132
## 4     1  1.3064702
## 5     1  0.8699767
## 6     1  1.1348831
## 7     2  0.6592557
## 8     2  0.5990497
## 9     2  0.7345132
## 10    2  0.1625562
## 11    2 -0.0632163
## 12    2  0.0812781
## 13    3 -0.1113811
## 14    3 -0.0060206
## 15    3  0.0812781
## 16    3 -0.1234223
## 17    3 -0.3431742
## 18    3 -0.2197519
## 
## $Result
##      Gene     dif     FC    LCL    UCL pvalue
## 1 C2H2-26  0.3592 0.4373 0.1926 0.9927 0.0488
## 2 C2H2-01 -0.6041 4.0185 2.4598 6.5649 0.0014
## 3 C2H2-12 -0.2167 1.6472 0.9595 2.8279 0.0624
```


### Generating plot
The `qpcrTTESTplot` function generates a bar plot of Fold Changes and confidence intervals for the target genes. the `qpcrTTESTplot` function accepts any gene name and any replicates. The `qpcrTTESTplot` function automatically puts appropriate signs of **, * on top of the plot columns based on the output p-values.



```r
# Producing the plot
t1 <- qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1,
              fontsizePvalue = 4)

# Producing the plot: specifying gene order
t2 <- qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1,
              order = c("C2H2-01", "C2H2-12", "C2H2-26"),
              paired = FALSE,
              var.equal = TRUE,
              width = 0.5,
              fill = "palegreen",
              y.axis.adjust = 0,
              y.axis.by = 2,
              ylab = "Average Fold Change (FC)",
              xlab = "none",
              fontsizePvalue = 4)

multiplot(t1, t2, cols = 2)
```

```
## $plot
## 
## $plot
```

```r
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-7-1.png" alt="Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests via ‘qpcrTTESTplot’ function."  />
<p class="caption">Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests via ‘qpcrTTESTplot’ function.</p>
</div>

## Analysis of covariance (ANCOVA)

analysis of covariance (ANCOVA) is a method based on both ANOVA and linear regression. It is basically suitable when the levels of a factor are also affected by an uncontrolled quantitative covariate. 
For example, suppose that wDCt of a target gene in a plant is affected by temperature. The gene may also be affected by drought. since we already know that temperature affects the target gene, we are interesting now if the gene expression is also altered by the drought levels. We can design an experiment to understand the gene behavior at both temperature and drought levels at the same time. The drought is another factor (the covariate) that may affect the expression of our gene under the levels of the first factor i.e. temperature. The data of such an experiment can be analyzed by ANCOVA or even ANOVA based on a factorial experiment using `qpcrANCOVA` function, if more than a factor exist. Bar plot of fold changes (FCs) along with the 95\% confidence interval is also returned by the `qpcrANCOVA` function. There is also a function called `oneFACTORplot` which returns FC values and related plot for a one-factor-experiment with more than two levels. 


```r
# See sample data
data_2factor
```

```
##    Genotype Drought Rep EPO  POCt EGAPDH GAPDHCt
## 1         R      D0   1   2 33.30      2   31.53
## 2         R      D0   2   2 33.39      2   31.57
## 3         R      D0   3   2 33.34      2   31.50
## 4         R      D1   1   2 32.73      2   31.30
## 5         R      D1   2   2 32.46      2   32.55
## 6         R      D1   3   2 32.60      2   31.92
## 7         R      D2   1   2 33.48      2   33.30
## 8         R      D2   2   2 33.27      2   33.37
## 9         R      D2   3   2 33.32      2   33.35
## 10        S      D0   1   2 26.85      2   26.94
## 11        S      D0   2   2 28.17      2   27.69
## 12        S      D0   3   2 27.99      2   27.39
## 13        S      D1   1   2 30.41      2   28.70
## 14        S      D1   2   2 29.49      2   28.66
## 15        S      D1   3   2 29.98      2   28.71
## 16        S      D2   1   2 29.03      2   30.61
## 17        S      D2   2   2 28.73      2   30.20
## 18        S      D2   3   2 28.83      2   30.49
```

```r
order <- unique(data_2factor$Drought)
qpcrANCOVA(data_2factor, 
           numberOfrefGenes = 1, 
           analysisType = "ancova",
           mainFactor.column = 2,
           mainFactor.level.order = order,
           fontsizePvalue = 4,
           x.axis.labels.rename = "none")
```

```
## $Final_data
##    Drought Genotype rep Etarget Cttarget Eref Ctref       wDCt
## 1       D0        R   1       2    33.30    2 31.53  0.5328231
## 2       D0        R   2       2    33.39    2 31.57  0.5478746
## 3       D0        R   3       2    33.34    2 31.50  0.5538952
## 10      D0        S   1       2    26.85    2 26.94 -0.0270927
## 11      D0        S   2       2    28.17    2 27.69  0.1444944
## 12      D0        S   3       2    27.99    2 27.39  0.1806180
## 4       D1        R   1       2    32.73    2 31.30  0.4304729
## 5       D1        R   2       2    32.46    2 32.55 -0.0270927
## 6       D1        R   3       2    32.60    2 31.92  0.2047004
## 13      D1        S   1       2    30.41    2 28.70  0.5147613
## 14      D1        S   2       2    29.49    2 28.66  0.2498549
## 15      D1        S   3       2    29.98    2 28.71  0.3823081
## 7       D2        R   1       2    33.48    2 33.30  0.0541854
## 8       D2        R   2       2    33.27    2 33.37 -0.0301030
## 9       D2        R   3       2    33.32    2 33.35 -0.0090309
## 16      D2        S   1       2    29.03    2 30.61 -0.4756274
## 17      D2        S   2       2    28.73    2 30.20 -0.4425141
## 18      D2        S   3       2    28.83    2 30.49 -0.4997098
## 
## $lm_ANOVA
## 
## Call:
## lm(formula = formula_ANOVA, data = x)
## 
## Coefficients:
## [1]  0.544864292151806   -0.342170761738059  -0.539847125557406
## [4]  -0.445524393582693  0.625138957662201   -0.032109866204158
## 
## 
## $lm_ANCOVA
## 
## Call:
## lm(formula = formula_ANCOVA, data = x)
## 
## Coefficients:
## [1]  0.446026110242133    -0.247848029763345   -0.0296012829069583
## [4]  -0.555902058659485 
## 
## 
## $ANOVA_table
## Analysis of Variance Table
## 
## Response: wDCt
##                  Df  Sum Sq Mean Sq F value    Pr(>F)    
## Drought           2 1.17379 0.58690  41.394 4.117e-06 ***
## Genotype          1 0.27643 0.27643  19.497 0.0008422 ***
## Drought:Genotype  2 0.41190 0.20595  14.526 0.0006239 ***
## Residuals        12 0.17014 0.01418                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $ANCOVA_table
## Analysis of Variance Table
## 
## Response: wDCt
##           Df  Sum Sq Mean Sq F value    Pr(>F)    
## Genotype   1 0.27643 0.27643   6.649 0.0218703 *  
## Drought    2 1.17379 0.58690  14.117 0.0004398 ***
## Residuals 14 0.58204 0.04157                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $FC_statistics_of_the_main_factor
##   contrast     FC pvalue sig       LCL      UCL    sddiff
## 1       D0      1 1.0000     0.0000000 0.000000 0.0000000
## 2 D1 vs D0 1.0705 0.8051  ns 0.5446459 2.104206 0.3466862
## 3 D2 vs D0 3.5967 0.0003 *** 1.6702550 7.744999 1.3629288
## 
## $FC_Plot_of_the_main_factor_levels
```

<img src="vignette_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />


## A target gene under more than two conditions (ANOVA)

Analysis of variance (ANOVA) of factorial experiments in the frame of a completely randomized design (CRD) can be done by the `qpcrANOVA` function. ANOVA of qPCR data is suitable when there is a factor with more than two levels, or when there is more than an experimental factor. The input data set should be prepared as shown below. Factor columns should be presented first followed by biological replicates and efficiency and Ct values of target and reference genes. The example data set below (`data_3factor`) represents amplification efficiency and Ct values for target and reference genes under three grouping factors (two different cultivars, three drought levels, and the presence or absence of bacteria). The table can contain any number of factor columns. The factor columns should be followed by five other columns assigned to biological replicates (r), the efficiency of the target gene, Ct values of the target gene, the efficiency of the reference gene, and Ct values of the reference gene, respectively. Here, the efficiency of 2 has been used for all wells, but the calculated efficiencies can be used instead.



```r
# See a sample dataset
data_3factor
```

```
##    Type Conc SA Replicate EPO  POCt EGAPDH GAPDHCt
## 1     R    L A1         1   2 33.30      2   31.53
## 2     R    L A1         2   2 33.39      2   31.57
## 3     R    L A1         3   2 33.34      2   31.50
## 4     R    L A2         1   2 34.01      2   31.48
## 5     R    L A2         2   2 36.82      2   31.44
## 6     R    L A2         3   2 35.44      2   31.46
## 7     R    M A1         1   2 32.73      2   31.30
## 8     R    M A1         2   2 32.46      2   32.55
## 9     R    M A1         3   2 32.60      2   31.92
## 10    R    M A2         1   2 33.37      2   31.19
## 11    R    M A2         2   2 33.12      2   31.94
## 12    R    M A2         3   2 33.21      2   31.57
## 13    R    H A1         1   2 33.48      2   33.30
## 14    R    H A1         2   2 33.27      2   33.37
## 15    R    H A1         3   2 33.32      2   33.35
## 16    R    H A2         1   2 32.53      2   33.47
## 17    R    H A2         2   2 32.61      2   33.26
## 18    R    H A2         3   2 32.56      2   33.36
## 19    S    L A1         1   2 26.85      2   26.94
## 20    S    L A1         2   2 28.17      2   27.69
## 21    S    L A1         3   2 27.99      2   27.39
## 22    S    L A2         1   2 28.71      2   29.45
## 23    S    L A2         2   2 29.01      2   29.46
## 24    S    L A2         3   2 28.82      2   29.48
## 25    S    M A1         1   2 30.41      2   28.70
## 26    S    M A1         2   2 29.49      2   28.66
## 27    S    M A1         3   2 29.98      2   28.71
## 28    S    M A2         1   2 28.91      2   28.09
## 29    S    M A2         2   2 28.60      2   28.65
## 30    S    M A2         3   2 28.59      2   28.37
## 31    S    H A1         1   2 29.03      2   30.61
## 32    S    H A1         2   2 28.73      2   30.20
## 33    S    H A1         3   2 28.83      2   30.49
## 34    S    H A2         1   2 28.29      2   30.84
## 35    S    H A2         2   2 28.53      2   30.65
## 36    S    H A2         3   2 28.28      2   30.74
```

The `qpcrANOVA` function performs ANOVA based on both factorial arrangement and completely randomized design (CRD). For the latter, a column of treatment combinations is made first as a grouping factor followed by ANOVA. You can call the input data set along with the added wCt and treatment combinations by `qpcrANOVA`. CRD-based analysis is especially useful when post-hoc tests and mean comparisons/grouping are desired for all treatment combinations. The final results along with the ANOVA tables can be called by `qpcrANOVA`.

### Reverse ordering of the grouping letters

One may be interested in presenting the statistical mean comparison result in the frame of grouping letters. This is rather challenging because in the grouping output of mean comparisons (via the `LSD.test` function of agricolae package), means are sorted into descending order so that the largest mean, is the first in the table and "a" letter is assigned to it. If `LSD.test` is applied to the wCt means, the biggest wCt mean receives "a" letter as expected, but this value turns into the smallest mean after its reverse log transformation by  $10^{-(\Delta Ct)}$. to solve this issue, I used a function that assigns the grouping letters appropriately.

### Output table of the analysis

The `qpcrANOVA` function produces the main analysis output including mean wDCt, LCL, UCL, grouping letters, and standard deviations. The standard deviation for each mean is derived from the back-transformed raw wDCt values from biological replicates for that mean.


```r
# If the data include technical replicates, means of technical replicates
# should be calculated first using meanTech function.

# Applying ANOVA analysis
res <- qpcrANOVA(data_2factor,
          numberOfrefGenes = 1,
          p.adj = "none")
res$Result
```

```
##      Genotype Drought      RE     LCL     UCL letters     std
## R:D0        R      D0 0.28519 0.40264 0.20200       d 0.00716
## R:D1        R      D1 0.62706 0.88530 0.44414      bc 0.35081
## R:D2        R      D2 0.98851 1.39563 0.70016       b 0.09786
## S:D0        S      D0 0.79554 1.12317 0.56347       b 0.21896
## S:D1        S      D1 0.41466 0.58543 0.29370      cd 0.12892
## S:D2        S      D2 2.96905 4.19183 2.10296       a 0.19549
```

```r
res$Post_hoc_Test
```

```
## NULL
```




```r
# Before plotting, the result needs to be extracted as below:
out2 <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)

f1 <- oneFACTORplot(out2,
              width = 0.2,
              fill = "skyblue",
              y.axis.adjust = 0.5,
              y.axis.by = 1,
              errorbar = "ci",
              show.letters = TRUE,
              letter.position.adjust = 0.1,
              ylab = "Relative Expression",
              xlab = "Factor Levels",
              fontsize = 12,
              fontsizePvalue = 4)


addline_format <- function(x,...){
    gsub('\\s','\n',x)
}
order <- unique(data_2factor$Drought)
f2 <- qpcrANCOVA(data_1factor,
                 numberOfrefGenes = 1,
                 mainFactor.column = 1,
                 mainFactor.level.order = c("L1","L2","L3"),
                 width = 0.5,
                 fill = c("skyblue","#79CDCD"),
                 y.axis.by = 1,
                 letter.position.adjust = 0,
                 y.axis.adjust = 1,
                 ylab = "Fold Change",
                 fontsize = 12,
                 x.axis.labels.rename = addline_format(c("Control", 
                                                       "Treatment_1 vs Control", 
                                                       "Treatment_2 vs Control")))


multiplot(f1, f2, cols = 2)
```

```
## $plot
## 
## $Final_data
##   SA rep Etarget Cttarget Eref Ctref       wDCt
## 1 L1   1       2    33.30    2 31.53  0.5328231
## 2 L1   2       2    33.39    2 31.57  0.5478746
## 3 L1   3       2    33.34    2 31.50  0.5538952
## 4 L2   1       2    32.73    2 31.30  0.4304729
## 5 L2   2       2    32.46    2 32.55 -0.0270927
## 6 L2   3       2    32.60    2 31.92  0.2047004
## 7 L3   1       2    33.48    2 33.30  0.0541854
## 8 L3   2       2    33.27    2 33.37 -0.0301030
## 9 L3   3       2    33.32    2 33.35 -0.0090309
## 
## $lm_ANOVA
## 
## Call:
## lm(formula = formula_ANOVA, data = x)
## 
## Coefficients:
## [1]  0.544864292151806   -0.342170761738059  -0.539847125557406
## 
## 
## $lm_ANCOVA
## 
## Call:
## lm(formula = formula_ANCOVA, data = x)
## 
## Coefficients:
## [1]  0.544864292151806   -0.342170761738059  -0.539847125557406
## 
## 
## $ANOVA_table
## Analysis of Variance Table
## 
## Response: wDCt
##           Df  Sum Sq  Mean Sq F value   Pr(>F)   
## SA         2 0.44759 0.223796  12.345 0.007473 **
## Residuals  6 0.10877 0.018129                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $ANCOVA_table
## Analysis of Variance Table
## 
## Response: wDCt
##           Df  Sum Sq  Mean Sq F value   Pr(>F)   
## SA         2 0.44759 0.223796  12.345 0.007473 **
## Residuals  6 0.10877 0.018129                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $FC_statistics_of_the_main_factor
##                   contrast     FC pvalue sig       LCL      UCL    sddiff
## 1                  Control      1 1.0000     0.0000000 0.000000 0.0000000
## 2 Treatment_1\nvs\nControl 2.1987 0.0208   * 0.5965179 8.104347 0.8762170
## 3 Treatment_2\nvs\nControl 3.4661 0.0027  ** 2.7448323 4.377019 0.2101155
## 
## $FC_Plot_of_the_main_factor_levels
```

```r
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-11-1.png" alt="A) A bar plot representing Relative expression of a gene under three levels of a factor generated using ‘oneFACTORplot’ function, B) Plot of average Fold changes produced by the ‘qpcrANCOVA’ function from the same data as ‘C’. Check level can be changed by user. Error bars represent 95% confidence interval."  />
<p class="caption">A) A bar plot representing Relative expression of a gene under three levels of a factor generated using ‘oneFACTORplot’ function, B) Plot of average Fold changes produced by the ‘qpcrANCOVA’ function from the same data as ‘C’. Check level can be changed by user. Error bars represent 95% confidence interval.</p>
</div>




### Barplot with the (1-alpha)% confidence interval as error bars


```r
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_2factor, numberOfrefGenes = 1)

# Plot of the 'res' data with 'Genotype' as grouping factor
q1 <- twoFACTORplot(res,
   x.axis.factor = Drought,
   group.factor = Genotype,
   width = 0.5,
   fill = "Greens",
   y.axis.adjust = 0.5,
   y.axis.by = 2,
   ylab = "Relative Expression",
   xlab = "Drought Levels",
   legend.position = c(0.15, 0.8),
   show.letters = TRUE,
   fontsizePvalue = 4)

# Plotting the same data with 'Drought' as grouping factor
q2 <- twoFACTORplot(res,
   x.axis.factor = Genotype,
   group.factor = Drought,
   xlab = "Genotype",
   fill = "Blues",
   legend.position = c(0.15, 0.8),
   show.letters = FALSE,
   show.errorbars = F,
   show.points = T,
   fontsizePvalue = 4)

multiplot(q1, q2, cols = 2)
```

```
## $plot
## 
## $plot
```

```r
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-12-1.png" alt="Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package."  />
<p class="caption">Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package.</p>
</div>





### A three-factorial experiment example

```r
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_3factor, numberOfrefGenes = 1)
res
```

```
## $Final_data
##    Type Conc SA rep Etarget Cttarget Eref Ctref       wDCt
## 1     R    L A1   1       2    33.30    2 31.53  0.5328231
## 2     R    L A1   2       2    33.39    2 31.57  0.5478746
## 3     R    L A1   3       2    33.34    2 31.50  0.5538952
## 4     R    L A2   1       2    34.01    2 31.48  0.7616059
## 5     R    L A2   2       2    36.82    2 31.44  1.6195414
## 6     R    L A2   3       2    35.44    2 31.46  1.1980994
## 7     R    M A1   1       2    32.73    2 31.30  0.4304729
## 8     R    M A1   2       2    32.46    2 32.55 -0.0270927
## 9     R    M A1   3       2    32.60    2 31.92  0.2047004
## 10    R    M A2   1       2    33.37    2 31.19  0.6562454
## 11    R    M A2   2       2    33.12    2 31.94  0.3552154
## 12    R    M A2   3       2    33.21    2 31.57  0.4936892
## 13    R    H A1   1       2    33.48    2 33.30  0.0541854
## 14    R    H A1   2       2    33.27    2 33.37 -0.0301030
## 15    R    H A1   3       2    33.32    2 33.35 -0.0090309
## 16    R    H A2   1       2    32.53    2 33.47 -0.2829682
## 17    R    H A2   2       2    32.61    2 33.26 -0.1956695
## 18    R    H A2   3       2    32.56    2 33.36 -0.2408240
## 19    S    L A1   1       2    26.85    2 26.94 -0.0270927
## 20    S    L A1   2       2    28.17    2 27.69  0.1444944
## 21    S    L A1   3       2    27.99    2 27.39  0.1806180
## 22    S    L A2   1       2    28.71    2 29.45 -0.2227622
## 23    S    L A2   2       2    29.01    2 29.46 -0.1354635
## 24    S    L A2   3       2    28.82    2 29.48 -0.1986798
## 25    S    M A1   1       2    30.41    2 28.70  0.5147613
## 26    S    M A1   2       2    29.49    2 28.66  0.2498549
## 27    S    M A1   3       2    29.98    2 28.71  0.3823081
## 28    S    M A2   1       2    28.91    2 28.09  0.2468446
## 29    S    M A2   2       2    28.60    2 28.65 -0.0150515
## 30    S    M A2   3       2    28.59    2 28.37  0.0662266
## 31    S    H A1   1       2    29.03    2 30.61 -0.4756274
## 32    S    H A1   2       2    28.73    2 30.20 -0.4425141
## 33    S    H A1   3       2    28.83    2 30.49 -0.4997098
## 34    S    H A2   1       2    28.29    2 30.84 -0.7676265
## 35    S    H A2   2       2    28.53    2 30.65 -0.6381836
## 36    S    H A2   3       2    28.28    2 30.74 -0.7405338
## 
## $lmCRD
## 
## Call:
## stats::lm(formula = wDCt ~ T, data = x)
## 
## Coefficients:
## (Intercept)      TR:H:A2      TR:L:A1      TR:L:A2      TR:M:A1      TR:M:A2  
##    0.005017    -0.244838     0.539847     1.188065     0.197676     0.496699  
##     TS:H:A1      TS:H:A2      TS:L:A1      TS:L:A2      TS:M:A1      TS:M:A2  
##   -0.477634    -0.720465     0.094323    -0.190652     0.377291     0.094323  
## 
## 
## $ANOVA_CRD
## Analysis of Variance Table
## 
## Response: wDCt
##           Df Sum Sq Mean Sq F value    Pr(>F)    
## T         11 8.5183 0.77439  29.188 3.248e-11 ***
## Residuals 24 0.6367 0.02653                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $Result
##        Type Conc SA      RE     LCL     UCL letters     std
## R:H:A1    R    H A1 0.98851 1.54552 0.63225      cd 0.09786
## R:H:A2    R    H A2 1.73708 2.71589 1.11104      bc 0.17469
## R:L:A1    R    L A1 0.28519 0.44589 0.18241       f 0.00716
## R:L:A2    R    L A2 0.06411 0.10023 0.04100       g 0.07728
## R:M:A1    R    M A1 0.62706 0.98039 0.40107      de 0.35081
## R:M:A2    R    M A2 0.31498 0.49246 0.20146       f 0.11049
## S:H:A1    S    H A1 2.96905 4.64204 1.89900      ab 0.19549
## S:H:A2    S    H A2 5.19335 8.11969 3.32167       a 0.78933
## S:L:A1    S    L A1 0.79554 1.24380 0.50883       d 0.21896
## S:L:A2    S    L A2 1.53333 2.39732 0.98072       c 0.15622
## S:M:A1    S    M A1 0.41466 0.64831 0.26522      ef 0.12892
## S:M:A2    S    M A2 0.79554 1.24380 0.50883       d 0.23677
```

```r
# releveling a factor levels first
res$Conc <- factor(res$Conc, levels = c("L","M","H"))
res$Type <- factor(res$Type, levels = c("S","R"))

# Arrange the first three colunms of the result table.
# This determines the columns order and shapes the plot output.
p1 <- threeFACTORplot(res,
    arrangement = c(3, 1, 2),
    legend.position = c(0.2, 0.85),
    xlab = "condition",
    fontsizePvalue = 4)


# When using ci as error, increase y.axis.adjust to see the plot correctly!
p2 <- threeFACTORplot(res,
   arrangement = c(2, 3, 1),
   bar.width = 0.8,
   fill = "Greens",
   xlab = "Drought",
   ylab = "Relative Expression",
   errorbar = "ci",
   y.axis.adjust = 8,
   y.axis.by = 2,
   letter.position.adjust = 0.6,
   legend.title = "Genotype",
   fontsize = 12,
   legend.position = c(0.2, 0.8),
   show.letters = TRUE,
   fontsizePvalue = 4)

multiplot(p1, p2, cols = 2)
```

```
## $plot
## 
## $plot
```

```r
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-13-1.png" alt="A and B) Relative expression (RE) of a target gene under two or three factors produced by ‘twoFACTORplot’ and ‘threeFACTORplot’ functions, respectively. Error bars represent standard deviations (can be set to confidence interval). Means (columns) lacking letters in common have significant differences at alpha = 0.05 as resulted from an ‘LSD.test’."  />
<p class="caption">A and B) Relative expression (RE) of a target gene under two or three factors produced by ‘twoFACTORplot’ and ‘threeFACTORplot’ functions, respectively. Error bars represent standard deviations (can be set to confidence interval). Means (columns) lacking letters in common have significant differences at alpha = 0.05 as resulted from an ‘LSD.test’.</p>
</div>


# An example of Showing point on the plot

```r
library(ggplot2)
b <- qpcrANOVA(data_3factor, numberOfrefGenes = 1)$Result
a <- qpcrANOVA(data_3factor, numberOfrefGenes = 1)$Final_data

ggplot(b, aes(x = Type, y = RE, fill = factor(Conc))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ SA) +
  scale_fill_brewer(palette = "Reds") +
  xlab("Type") +
  ylab("Relative Expression") +
  geom_point(data = a, aes(x = Type, y = (10^(-wDCt)), fill = factor(Conc)), 
             position = position_dodge(width = 0.9), color = "black") +
  ylab("ylab") +
  xlab("xlab") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
        axis.title  = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  theme(legend.position  = c(0.2, 0.7)) +
  theme(legend.title = element_text(size = 12, color = "black")) +
  scale_y_continuous(breaks = seq(0, max(b$RE) + max(b$std) + 0.1, by = 5), 
                     limits = c(0, max(b$RE) + max(b$std) + 0.1), expand = c(0, 0)) 
```



# Checking normality of residuals

If the residuals from a `t.test` or a `lm` object are not normally distributed, the grouping letters (deduced from the `LSD.test`) might be violated. In such cases, one could apply another data transformation to the wDCt data for ANOVA and mean comparison purposes or use non-parametric tests such as the Mann-Whitney test (also known as the Wilcoxon rank-sum test), `wilcox.test()`, which is an alternative to `t.test`, or the `kruskal.test()` test which alternative to one-way analysis of variance, to test the difference between medians of the populations using independent samples. However, the `t.test` function (along with the `qpcrTTEST` function described above) includes the `var.equal` argument. When set to `FALSE`, perform `t.test` under the unequal variances hypothesis.


```r
residuals <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)$lmCRD$residuals
shapiro.test(residuals) 
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals
## W = 0.84232, p-value = 0.06129
```

```r
qqnorm(residuals)
qqline(residuals, col = "red")
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-15-1.png" alt="QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."  />
<p class="caption">QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions.</p>
</div>

# Mean of technical replicates
Calculating the mean of technical replicates and getting an output table appropriate for subsequent ANOVA analysis can be done using the `meanTech` function. For this, the input data set should follow the column arrangement of the following example data. Grouping columns must be specified under the `groups` argument of the `meanTech` function.


```r
# See example input data frame:
data_withTechRep
```

```
##    factor1 factor2 factor3 biolrep techrep Etarget targetCt Eref  refCt
## 1    Line1    Heat    Ctrl       1       1       2   33.346    2 31.520
## 2    Line1    Heat    Ctrl       1       2       2   28.895    2 29.905
## 3    Line1    Heat    Ctrl       1       3       2   28.893    2 29.454
## 4    Line1    Heat    Ctrl       2       1       2   30.411    2 28.798
## 5    Line1    Heat    Ctrl       2       2       2   33.390    2 31.574
## 6    Line1    Heat    Ctrl       2       3       2   33.211    2 31.326
## 7    Line1    Heat    Ctrl       3       1       2   33.845    2 31.759
## 8    Line1    Heat    Ctrl       3       2       2   33.345    2 31.548
## 9    Line1    Heat    Ctrl       3       3       2   32.500    2 31.477
## 10   Line1    Heat   Treat       1       1       2   33.006    2 31.483
## 11   Line1    Heat   Treat       1       2       2   32.588    2 31.902
## 12   Line1    Heat   Treat       1       3       2   33.370    2 31.196
## 13   Line1    Heat   Treat       2       1       2   36.820    2 31.440
## 14   Line1    Heat   Treat       2       2       2   32.750    2 31.300
## 15   Line1    Heat   Treat       2       3       2   32.450    2 32.597
## 16   Line1    Heat   Treat       3       1       2   35.238    2 31.461
## 17   Line1    Heat   Treat       3       2       2   28.532    2 30.651
## 18   Line1    Heat   Treat       3       3       2   28.285    2 30.745
```

```r
# Calculating mean of technical replicates
meanTech(data_withTechRep, groups = 1:4)
```

```
##   factor1 factor2 factor3 biolrep Etarget targetCt Eref    refCt
## 1   Line1    Heat    Ctrl       1       2 30.37800    2 30.29300
## 2   Line1    Heat    Ctrl       2       2 32.33733    2 30.56600
## 3   Line1    Heat    Ctrl       3       2 33.23000    2 31.59467
## 4   Line1    Heat   Treat       1       2 32.98800    2 31.52700
## 5   Line1    Heat   Treat       2       2 34.00667    2 31.77900
## 6   Line1    Heat   Treat       3       2 30.68500    2 30.95233
```

# Citation

```r
citation("rtpcr")
```


# Contact 
Email: gh.mirzaghaderi at uok.ac.ir


# References
Livak, Kenneth J, and Thomas D Schmittgen. 2001. Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR and the Double Delta CT Method. Methods 25 (4). <a href="https://doi.org/10.1006/meth.2001.1262">doi.org/10.1006/meth.2001.1262</a>.


Ganger, MT, Dietz GD, and Ewing SJ. 2017. A common base method for analysis of qPCR data and the application of simple blocking in qPCR experiments. BMC bioinformatics 18, 1-11. <a href="https://doi.org/10.1186/s12859-017-1949-5">doi.org/10.1186/s12859-017-1949-5</a>.


Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006. Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85). <a href="https://doi.org/10.1186/1471-2105-7-85">doi.org/10.1186/1471-2105-7-85</a>.


.
