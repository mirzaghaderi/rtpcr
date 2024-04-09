---
title: "rtpcr: a package for statistical analysis of real-time PCR data in R"
author: Ghader Mirzaghaderi
output:
  html_document: 
    toc: yes
    keep_md: yes
    output: rmarkdown::html_vignette
    df_print: default
  word_document:
    toc: No
  pdf_document: 
    toc: yes
    latex_engine: xelatex
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
When there are only a two conditional factor (e.g. with treatment and control levels), or there is one multi-level factor that one of them is concidered as control, average fold change (FC) expression of target gene can be calculated according to:

$$\text{Fold Change}=10^{-(\overline{w\Delta Ct}_{\text{Tr}}-{\overline{w\Delta Ct}_{\text{Co}}})}$$

in rtpcr package, ‘qpcrTTEST’, ‘qpcrTTESTplot’ finctions calcualtes FC for multi-genes-two conditional cases and ‘qpcrANCOVA’ represents FC for single-gene-single factor or single-gene-multi factorial experiments. If wDCt values is calculated from the E valuses, these calculations match the formula of Pfaffl while if 2 (complete efficiency) be used instead, the result match the $2^{-\Delta\Delta C_t}$ method. In any case we called these as Fold Change in the outputs of ‘rtpcr’. Under factorial experiments where the calculation of the expression of the target gene relative to the reference gene (called Relative Expression) in each condition is desired, ‘qpcrANOVA’, ‘oneFACTORplot’, ‘twoFACTORplot’ and ‘threeFACTORplot’ functions were developed for ANOVA analysis, and representing the plots from single, double or triple factor experiments, respectively. The last three functions generate ‘ggplot2’-derived graphs based on the output of the ‘qpcrANOVA’ function. If available, the blocking factor can laso be handled by ‘qpcrTTEST’ and ‘qpcrANOVA’ functions.
Here, a brief methodology is presented but detailes about the $w\Delta C_t$  calculations and statistical analysis are available in. Importantly, because the relative or FC gene expression follow a lognormal distribution, a normal distribution is expected for the $w\Delta C_t$ or $w\Delta\Delta C_t$ values making it possible to apply t-test or analysis of variance to them. Following analysis, $w\Delta C_t$ values are statistically compared and standard deviations and confidence intervals is calculated, but the transformation $y = 10^x$ is applied in the final step in order to report the results (i.e. RE ratios, errors and confidence limits).





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

| Experiment type   |                  Column arrangement of the input data | Example in the package |
 |:---------------------|:-----------------------------------|:----------------------------------|
 |Amplification efficiency             |Dilutions - targetCt - refCt | data_efficiency |
 |t-test (accepts multiple genes)      |condition (put the control level first) - gene (put reference gene(s) last.)- efficiency - Ct  | data_ttest |
 |Factorial (Up to three factors)      |factor1 - rep - targetE - targetCt - refE - refCt | data_1factor |
 |                                     |factor1 - factor2 - rep - targetE - targetCt - refE - refCt | data_2factor |
 |                                     |factor1 - factor2 - factor3 - rep - targetE - targetCt - refE - refCt | data_3factor_b |
 |Factorial with blocking              |factor1 - block - rep - targetE - targetCt - refE - refCt | |
 |                                     |factor1 - factor2 - block - rep - targetE - targetCt - refE - refCt	 | data_2factorBlock |
 |                                     |factor1 - factor2 - factor3 - block - rep - targetE - targetCt - refE - refCt | |
 |Two reference genes                  |. . . . . .  rep - targetE - targetCt - ref1E - ref1Ct - ref2E - ref2Ct | |
 |calculating biological replicated    |. . . . . .  biologicalRep - techcicalRep - Etarget - targetCt - Eref - refCt  | data_withTechRep |
 |                                       |. . . . . .  biologicalRep - techcicalRep - Etarget - targetCt - ref1E - ref1Ct - ref2E - ref2Ct  | |
 

# functions usage

To simplify 'rtpcr' usage, examples for using the functions are presented below.


           
*Table 2. Functions and examples for using them.*
| function   |                 Analysis | Example (see package help for the more arguments) |
 |:---------------------|:-----------------------------------|:----------------------------------|
 | efficiency             | Efficiency, standard curves and related statistics | efficiency(data_efficiency) |
 | meanTech      | Calculating the mean of technical replicates | meanTech(data_withTechRep, groups = 1:4) |
 | qpcrANCOVA      | FC and bar plot of the target gene (one or multi-factorial experiments) | qpcrANCOVA(data_1factor, numberOfrefGenes = 1, mainFactor.column = 1, mainFactor.level.order = c("L1", "L2", "L3") |
 |  oneFACTORplot    | Bar plot of the relative gene expression from a one-factor experiment | out <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)$Result;   oneFACTORplot(out) |
 |  qpcrANOVA  | Analysis of Variance of the qpcr data  | qpcrANOVA(data_3factor_a, numberOfrefGenes = 1, p.adj = "none")|
 | qpcrTTEST     | Computing the average fold change and related statistics | qpcrTTEST(data_ttest,  numberOfrefGenes = 1, paired = FALSE, var.equal = TRUE) |
 | qpcrTTESTplot  | Bar plot of the average fold change of the target genes	 | qpcrTTESTplot(data_ttest,  numberOfrefGenes = 1, order = c("C2H2-01", "C2H2-12", "C2H2-26")) |
 |  threeFACTORplot  | Bar plot of the relative gene expression from a three-factor experiment | res <- qpcrANOVA(data_3factor_b,  numberOfrefGenes = 1)$Result; threeFACTORplot(res, arrangement = c(3, 1, 2)) |
 | twoFACTORplot   | Bar plot of the relative gene expression from a two-factor experiment | res <- qpcrANOVA(data_2factor,  numberOfrefGenes = 1)$Result; twoFACTORplot(res, x.axis.factor = Genotype, group.factor = Drought) |

 
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
              xlab = "Gene",
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





















