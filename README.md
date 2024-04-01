# iqpcr: a package for statistical analysis of real-time PCR data in R 

Real-time polymerase chain reaction (real-time PCR), is widely used in biological research. Various analysis methods are employed on the real-time PCR data to measure the mRNA levels under different experimental conditions. 
‘iqpcr’ package was developed for amplification efficiency calculation and statistical analysis of real-time PCR data in R. By accounting for up to two reference genes and amplification efficiency values, a general calculation methodology covering both the Livak and Pfaffl methods was developed. Based on the experimental conditions, the functions of the ‘iqpcr’ package use a t-test (for experiments with a two-level factor) or analysis of variance (for cases where more than two levels or factors or a blocking factor exist) to calculate the fold change (FC) or relative expression (RE). The functions further provide standard deviations and confidence limits for means, apply statistical mean comparisons and present letter mean grouping. To facilitate function application, different data sets were used as examples and the outputs were explained. An outstanding feature of ‘iqpcr’ package is providing publication-ready bar plots with various controlling arguments for experiments with up to three different factors. 
The ‘iqpcr’ package was built based on a general method with various capabilities. It is user-friendly and easy to work with and provides an applicable resource for analyzing real-time PCR data in R. ‘iqpcr’ is a CRAN package although the development version of the package can be obtained through Github.


# Overview

real-time PCR, also known as quantitative PCR (qPCR), is a powerful analytical tool that has revolutionized nucleic acid quantification. Among the various approaches developed for data analysis in real-time PCR, the Livak method (Livak and Schmittgen, 2001), also known as the $2^{-ΔΔCt}$ method, stands out for its simplicity and widespread use. This method assumes that both the target and reference genes are amplified with efficiencies close to 100 %. On the other hand, the Pfaffl method (Pfaffl et al., 2002) offers a more flexible approach by accounting for differences in amplification efficiencies between the target and reference genes. This method adjusts the calculated expression ratio by incorporating the specific amplification efficiencies, thus providing a more accurate representation of the relative gene expression levels (Pfaffl et al., 2002). Both methods have their merits and limitations, and the choice between them depends on the experimental design and the precision required for the study. This paper aims to delve into the mathematical underpinnings of these methodologies, providing a comprehensive understanding of their applications and implications in real-time PCR analysis. Among the various approaches developed for data analysis in real-time PCR, the Livak method, also known as the $2^{-\Delta\Delta C_t}$ method, stands out for its simplicity and widespread use.


$$\text{Fold Change} = $$

$$= \frac{2^{-(Ct_{target}-Ct_{ref})Tr}}{2^{-(Ct_{target}-Ct_{ref})Co}}$$

$$= 2^{[-(Ct_{target}-Ct_{ref})Tr - (Ct_{target}-Ct_{ref})Co]}$$

$$= 2^{-(ΔCt_{Tr} - ΔCt_{Co})}$$

$$= 2^{-ΔΔCt}$$



where Tr is Treatment and Co is Control conditions, respectively. This method assumes that both the target and reference genes are amplified with efficiencies close to 100%, allowing for the relative quantification of gene expression levels (Livak and Schmittgen, 2001). On the other hand, the Pfaffl method offers a more flexible approach by accounting for differences in amplification efficiencies between the target and reference genes. This method adjusts the calculated expression ratio by incorporating the specific amplification efficiencies, thus providing a more accurate representation of the relative gene expression levels (Pfaffl et al., 2002).


$$\text{Fold Change} = \frac{E^{-(Ct_{Tr}-Ct_{Co})target}}{E^{-(Ct_{Tr}-Ct_{Co})ref}}$$

# A generalized calculation method

ΔC_t is the difference between two Ct values (e.g. Cttarget−Ctref). The iqpcr packager functions are mainly based on the calculation of efficiency-weighted ΔC_t (wΔC_t) values from target and reference gene Ct values which are weighted for their amplification efficiencies as below:


$$wΔCt = log_{10}(E_{target}).Ct_{target} - log_{10}(E_{ref}).Ct_{ref}$$

From the mean wΔCt values over biological replicates, the relative expression of a target gene can be calculated for each condition as:

$$\text{Relative Expression} = 10^{-\overline{w\Delta Ct}}$$

When there are only two conditions (e.g. Treatment and Control), the average fold change expression of a target gene can be calculated as follows: 

$$\text{Fold Change} = 10^{-(\overline{w\Delta Ct}{Tr}-\overline{w\Delta Ct}{Co})}$$

if wDCt values are calculated from the E values, these calculations match the formula of Pfaffl while if 2 (complete efficiency) is used instead the result matches the $2^{-ΔΔCt}$ method. Here, a brief methodology is presented but details about the wΔC_t  calculations and statistical analysis are available in (Ganger et al., 2017). What is important here is that although the relative or fold change gene expression follows a lognormal distribution, a normal distribution is expected for the wΔCt or wΔΔCt values making it possible to apply t-test or analysis of variance to them. Following analysis, wΔCt values are statistically compared and standard deviations and confidence intervals are calculated, but the transformation y = 10x is applied in the final step to report the relative expression ratios, errors, and confidence limits.

# Installing and loading



<!-- `iqpcr` is a CRAN package and can be installed through `install.packages()`.
install.packages("iqpcr")
library(iqpcr) -->

The development version of the `iqpcr` package can be obtained through: 

```r
# Install from github
devtools::install_github("mirzaghaderi/iqpcr")

# Loading the package
library(iqpcr)
```


# Data structure and column arrangement

To use the functions, input data should be prepared in the right format with appropriate column arrangement. The correct column arrangement is shown in Table 1.

*Table 1. Data structure and column arrangement required for ‘iqpcr’ package.  rep: technical replicate; targetE and refE: amplification efficiency columns for target and reference genes respectively. targetCt and refCt: target and reference Ct columns, respectively. factors (factor1, factor2 and/or factor3): experimental factors.*
| Experiment type   |                  Column arrangement of the input data | Example in the package |
 |:---------------------|:-----------------------------------|:----------------------------------|
 |Amplification efficiency             |Dilutions - targetCt - refCt | data_efficiency |
 |t-test (accepts multiple genes)      |condition (put the control level first) - efficiency - gene (put reference gene(s) last.) - Ct  | data_ttest |
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

To simplify 'iqpcr, usage, examples for using the functions are presented below.

*Table 2. Functions and examples for using them.*
| function   |                 Analysis | Example (see package help for the more arguments) |
 |:---------------------|:-----------------------------------|:----------------------------------|
 | efficiency             | Efficiency, standard curves and related statistics | efficiency(data_efficiency) |
 | meanTech      | Calculating the mean of technical replicates | meanTech(data_withTechRep, groups = 1:4) |
 | oneFACTORfcplot      | Bar plot of the average fold change of one target gene with two or more levels | oneFACTORfcplot(data_1factor, levels = c(3, 2, 1), numberOfrefGenes = 1, level.names = c("A1", "A2", "A3")) |
 |  oneFACTORplot    | Bar plot of the relative gene expression from a one-factor experiment | out <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)$Result;   oneFACTORplot(out) |
 |  qpcrANOVA  | Analysis of Variance of the qpcr data  | qpcrANOVA(data_3factor_a, numberOfrefGenes = 1, p.adj = "none")|
 | qpcrTTEST     | Computing the average fold change and related statistics | qpcrTTEST(data_ttest,  numberOfrefGenes = 1, paired = FALSE, var.equal = TRUE) |
 | qpcrTTESTplot  | Bar plot of the average fold change of the target genes	 | qpcrTTESTplot(data_ttest,  numberOfrefGenes = 1, order = c("C2H2-01", "C2H2-12", "C2H2-26")) |
 |  threeFACTORplot  | Bar plot of the relative gene expression from a three-factor experiment | res <- qpcrANOVA(data_3factor_b,  numberOfrefGenes = 1)$Result; threeFACTORplot(res, arrangement = c(3, 1, 2)) |
 | twoFACTORplot   | Bar plot of the relative gene expression from a two-factor experiment | res <- qpcrANOVA(data_2factor,  numberOfrefGenes = 1)$Result; twoFACTORplot(res, x.axis.factor = Genotype, group.factor = Drought) |
 
*see package help for more arguments including the number of reference genes, levels arrangement, blocking, and arguments for adjusting the bar plots.*
 
   

# Amplification efficiency data analysis
## Sample data of amplification efficiency

To calculate the amplification efficiencies
of a target and a reference gene, a data frame should be prepared with 3 columns 
of dilutions, target gene Ct values, and reference gene Ct values, respectively, 
as shown below.


```r
data_efficiency
```

```
dilutions  C2H2.26    GAPDH
     1.00 25.57823 22.60794
     1.00 25.53636 22.68348
     1.00 25.50280 22.62602
     0.50 26.70615 23.67162
     0.50 26.72720 23.64855
     0.50 26.86921 23.70494
     0.20 28.16874 25.11064
     0.20 28.06759 25.11985
     0.20 28.10531 25.10976
     0.10 29.19743 26.16919
     0.10 29.49406 26.15119
     0.10 29.07117 26.15019
     0.05 30.16878 27.11533
     0.05 30.14193 27.13934
     0.05 30.11671 27.16338
     0.02 31.34969 28.52016
     0.02 31.35254 28.57228
     0.02 31.34804 28.53100
     0.01 32.55013 29.49048
     0.01 32.45329 29.48433
     0.01 32.27515 29.26234
```

## Calculating amplification efficiency
The following `efficiency` function calculates the amplification efficiency of a target and a
reference gene, and presents the related standard curves along with the Slope, Efficiency, and R2 statistics. The function also compares the slopes of the two standard curves. For this, a regression line is fitted using the $\Delta C_t$ values. If $2^{-\Delta\Delta C_t}$ method is intended, the slope should not exceed 0.2!


```r
efficiency(data_efficiency)
```

<img src="https://raw.githubusercontent.com/mirzaghaderi/iqpcr/c02e476cf7c69783fa9c65bc625a7ea4a0e3abcf/Figure%201.jpg">

*Figure 1 - Standard curve and the amplification efficiency analysis of target and reference genes.*


```
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
when a target gene is assessed under two different conditions (for example Control and treatment), it is possible to calculate the average fold change expression i.e. $2^{-\Delta \Delta C_t}$ of the target gene in treatment relative to control conditions. For this, the data should be prepared according to the following dataset consisting of 4 columns belonging to condition levels, E (efficiency), genes, and Ct values, respectively. Each Ct value is the mean of technical replicates. Complete amplification efficiencies of 2 have been assumed here for all wells but the calculated efficiencies can be used instead. 


```r
data_ttest
```

```
Condition E    Gene    Ct
  control 2 C2H2-26 31.26
  control 2 C2H2-26 31.01
  control 2 C2H2-26 30.97
treatment 2 C2H2-26 32.65
treatment 2 C2H2-26 32.03
treatment 2 C2H2-26 32.40
  control 2 C2H2-01 31.06
  control 2 C2H2-01 30.41
  control 2 C2H2-01 30.97
treatment 2 C2H2-01 28.85
treatment 2 C2H2-01 28.93
treatment 2 C2H2-01 28.90
  control 2 C2H2-12 28.50
  control 2 C2H2-12 28.40
  control 2 C2H2-12 28.80
treatment 2 C2H2-12 27.90
treatment 2 C2H2-12 28.00
treatment 2 C2H2-12 27.90
  control 2     ref 28.87
  control 2     ref 28.42
  control 2     ref 28.53
treatment 2     ref 28.31
treatment 2     ref 29.14
treatment 2     ref 28.63
```


Here, the above data set was used for the Fold Change expression analysis of the target genes using the `qpcrTTEST` function. This function performs a t-test-based analysis of any number of genes that 
have been evaluated under control and treatment conditions. The analysis can be done for unpaired or paired conditions. The output is a table of target gene names, fold changes confidence limits and the t.test derived p-values. The `qpcrTTEST` function includes the `var.equal` argument. When set to `FALSE`,
`t.test` is performed under the unequal variances hypothesis.


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
##      Gene     dif Fold_Change Lower.Er Upper.Er pvalue
## 1 C2H2-26  0.3592      0.4373   0.1926   0.9927 0.0488
## 2 C2H2-01 -0.6041      4.0185   2.4598   6.5649 0.0014
## 3 C2H2-12 -0.2167      1.6472   0.9595   2.8279 0.0624
```


### Generating plot
`qpcrTTESTplot` function generates a bar plot of Fold Changes and confidence intervals for the target genes. the `qpcrTTESTplot` function uses raw Ct data and accepts any gene name and any replicates. `qpcrTTESTplot` function automatically puts appropriate signs of **, * on top of the plot columns based on the output p-values.


```r
# Producing the plot
qpcrTTESTplot(data_ttest, numberOfrefGenes = 1)
```


```r
# Producing the plot: specifying gene order
qpcrTTESTplot(data_ttest,
   numberOfrefGenes = 1,
   order = c("C2H2-01", "C2H2-12", "C2H2-26"),
   paired = FALSE,
   var.equal = TRUE,
   width = 0.5,
   fill = "skyblue",
   y.axis.adjust = 0,
   y.axis.by = 2,
   ylab = "Average Fold Change",
   xlab = "Gene")
```


## A target gene under more than two conditions (ANOVA)

Analysis of variance (ANOVA) of factorial experiments in the frame of a completely randomized design (CRD) can be done by the `qpcrANOVA` function. ANOVA of qPCR data is suitable when there is a factor with more than two levels, or when there is more than one experimental factor. The input data set should be prepared as shown below. Factor columns should be presented first followed by biological replicates and efficiency and Ct values of target and reference genes. The example data set below (`data_3factor_a`) represents amplification efficiency and Ct values for target and reference genes under three grouping factors (two different cultivars, three drought levels, and the presence or absence of bacteria). The table can contain any number of factor columns. The factor columns should be followed by five other columns assigned to biological replicates (r), the efficiency of the target gene, Ct values of the target gene, the efficiency of the reference gene, and Ct values of the reference gene, respectively. Here, the efficiency of 2 has been used for all wells, but the calculated efficiencies can be used instead.



```r
# See a sample dataset
data_3factor_a
```

```
Genotype Drought SA Rep   EPO  POCt EGAPDH GAPDHCt
       R    0.00 A1   1 1.839 33.30  1.918   31.53
       R    0.00 A1   2 1.839 33.39  1.918   31.57
       R    0.00 A1   3 1.839 33.34  1.918   31.50
       R    0.00 A2   1 1.839 34.01  1.918   31.48
       R    0.00 A2   2 1.839 36.82  1.918   31.44
       R    0.00 A2   3 1.839 35.44  1.918   31.46
       R    0.25 A1   1 1.839 32.73  1.918   31.30
       R    0.25 A1   2 1.839 32.46  1.918   32.55
       R    0.25 A1   3 1.839 32.60  1.918   31.92
       R    0.25 A2   1 1.839 33.37  1.918   31.19
       R    0.25 A2   2 1.839 33.12  1.918   31.94
       R    0.25 A2   3 1.839 33.21  1.918   31.57
       R    0.50 A1   1 1.839 33.48  1.918   33.30
       R    0.50 A1   2 1.839 33.27  1.918   33.37
       R    0.50 A1   3 1.839 33.32  1.918   33.35
       R    0.50 A2   1 1.839 32.53  1.918   33.47
       R    0.50 A2   2 1.839 32.61  1.918   33.26
       R    0.50 A2   3 1.839 32.56  1.918   33.36
       S    0.00 A1   1 1.839 26.85  1.918   26.94
       S    0.00 A1   2 1.839 28.17  1.918   27.69
       S    0.00 A1   3 1.839 27.99  1.918   27.39
       S    0.00 A2   1 1.839 28.71  1.918   29.45
       S    0.00 A2   2 1.839 29.01  1.918   29.46
       S    0.00 A2   3 1.839 28.82  1.918   29.48
       S    0.25 A1   1 1.839 30.41  1.918   28.70
       S    0.25 A1   2 1.839 29.49  1.918   28.66
       S    0.25 A1   3 1.839 29.98  1.918   28.71
       S    0.25 A2   1 1.839 28.91  1.918   28.09
       S    0.25 A2   2 1.839 28.60  1.918   28.65
       S    0.25 A2   3 1.839 28.59  1.918   28.37
       S    0.50 A1   1 1.839 29.03  1.918   30.61
       S    0.50 A1   2 1.839 28.73  1.918   30.20
       S    0.50 A1   3 1.839 28.83  1.918   30.49
       S    0.50 A2   1 1.839 28.29  1.918   30.84
       S    0.50 A2   2 1.839 28.53  1.918   30.65
       S    0.50 A2   3 1.839 28.28  1.918   30.74
```

`qpcrANOVA` function performs  ANOVA based on both factorial arrangement and completely randomized design (CRD). For the latter, a column of treatment combinations is made first as a grouping factor followed by ANOVA. You can call the input data set along with the added wCt and treatment combinations by `qpcrANOVA(data_3factor_a)$Final_data`. CRD-based analysis is especially useful when post-hoc tests and mean comparisons/grouping is desired for all treatment combinations. The final results along with the ANOVA tables can be called by `qpcrANOVA(data_3factor_a)`.

### Reverse ordering of the grouping letters

One may be interested in presenting the statistical mean comparison result in the frame of grouping letters. This is rather challenging because in the grouping output of mean comparisons (via `LSD.test` function of agricolae package), means are sorted into descending order so that the largest mean, is the first in the
table and "a" letter is assigned to it. If `LSD.test` is applied to the wCt means,
the biggest wCt mean receives "a" letter as expected, but this value turns into
the smallest mean after its reverse log transformation by  $10^{-(\Delta Ct)}$. to solve
this issue, I used a function that assigns the grouping letters appropriately.

### Output table of the analysis

The `qpcrANOVA` function produces the main analysis output including mean wDCt, LCL, UCL, grouping letters, and standard deviations. The standard deviation for each mean is derived from the back-transformed raw wDCt values from biological replicates for that mean. If the data includes technical replicates, the means of technical replicates should be calculated first using `ManTech` function.


```r
# Applying ANOVA analysis
qpcrANOVA(data_3factor_a, 
	  numberOfrefGenes = 1,
	  p.adj = "none")
```

```
## $ANOVA_factorial
## Analysis of Variance Table
## 
## Response: wDCt
##                     Df  Sum Sq Mean Sq F value    Pr(>F)    
## Genotype             1 1.32885 1.32885 23.9559 3.701e-05 ***
## Drought              1 3.04576 3.04576 54.9076 4.537e-08 ***
## SA                   1 0.00400 0.00400  0.0720  0.790365    
## Genotype:Drought     1 0.21286 0.21286  3.8373  0.060151 .  
## Genotype:SA          1 0.47334 0.47334  8.5332  0.006821 ** 
## Drought:SA           1 0.19253 0.19253  3.4708  0.072982 .  
## Genotype:Drought:SA  1 0.27533 0.27533  4.9635  0.034099 *  
## Residuals           28 1.55318 0.05547                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $ANOVA_CRD
## Analysis of Variance Table
## 
## Response: wDCt
##           Df Sum Sq Mean Sq F value    Pr(>F)    
## T         11 6.5790 0.59809  28.323 4.502e-11 ***
## Residuals 24 0.5068 0.02112                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $Result
##           Genotype Drought SA    wDCt     LCL     UCL letters    std
## R:0.25:A1        R    0.25 A1  2.5409  3.7857  1.7054      de 1.3099
## R:0.25:A2        R    0.25 A2  1.3666  2.0362  0.9173       f 0.4422
## R:0.5:A1         R     0.5 A1  4.0235  5.9947  2.7005      cd 0.3568
## R:0.5:A2         R     0.5 A2  6.6104  9.8488  4.4368      bc 0.6135
## R:0:A1           R       0 A1  1.2506  1.8633  0.8394       f 0.0280
## R:0:A2           R       0 A2  0.3358  0.5003  0.2254       g 0.3414
## S:0.25:A1        S    0.25 A1  1.5419  2.2973  1.0349      ef 0.4179
## S:0.25:A2        S    0.25 A2  2.6972  4.0186  1.8103      de 0.7382
## S:0.5:A1         S     0.5 A1  9.3608 13.9467  6.2829      ab 0.6034
## S:0.5:A2         S     0.5 A2 15.5027 23.0975 10.4052       a 2.1351
## S:0:A1           S       0 A1  2.5829  3.8482  1.7336      de 0.5779
## S:0:A2           S       0 A2  5.0276  7.4906  3.3745       c 0.4507
## 
## $Post_hoc_Test
##                            FC pvalue signif.     LCL     UCL
## R:0.25:A1 - R:0.25:A2  1.8592 0.0325       *  1.0579  3.2675
## R:0.25:A1 - R:0.5:A1   0.6315 0.1055          0.3593  1.1098
## R:0.25:A1 - R:0.5:A2   0.3844 0.0018      **  0.2187  0.6755
## R:0.25:A1 - R:0:A1     2.0317 0.0159       *  1.1561  3.5706
## R:0.25:A1 - R:0:A2     7.5672 0.0000     ***  4.3058 13.2990
## R:0.25:A1 - S:0.25:A1  1.6479 0.0800       .  0.9377  2.8961
## R:0.25:A1 - S:0.25:A2  0.9420 0.8288          0.5360  1.6556
## R:0.25:A1 - S:0.5:A1   0.2714 0.0001     ***  0.1545  0.4770
## R:0.25:A1 - S:0.5:A2   0.1639 0.0000     ***  0.0933  0.2880
## R:0.25:A1 - S:0:A1     0.9837 0.9527          0.5598  1.7289
## R:0.25:A1 - S:0:A2     0.5054 0.0197       *  0.2876  0.8882
## R:0.25:A2 - R:0.5:A1   0.3397 0.0006     ***  0.1933  0.5969
## R:0.25:A2 - R:0.5:A2   0.2067 0.0000     ***  0.1176  0.3633
## R:0.25:A2 - R:0:A1     1.0928 0.7482          0.6218  1.9205
## R:0.25:A2 - R:0:A2     4.0701 0.0000     ***  2.3159  7.1530
## R:0.25:A2 - S:0.25:A1  0.8863 0.6627          0.5043  1.5577
## . . . 
```

### Barplot with the (1-alpha)% confidence interval as error bars

Although the plot of fold changes from a single factor experiment with two levels is done directly from the raw Ct data, before plotting data from two or three factorial experiments, the result table needs to be extracted as below. 
```r
res <- qpcrANOVA(data_2factor, numberOfrefGenes = 1)$Result
res


##                Q75 Genotype Drought   wDCt    LCL    UCL letters    std
## R:0     0.55088489        R       0 0.2852 0.4026 0.2020       d 0.0072
## R:0.25  0.31758665        R    0.25 0.6271 0.8853 0.4441      bc 0.3508
## R:0.5   0.02257725        R     0.5 0.9885 1.3956 0.7002       b 0.0979
## S:0     0.16255620        S       0 0.7955 1.1232 0.5635       b 0.2190
## S:0.25  0.44853469        S    0.25 0.4147 0.5854 0.2937      cd 0.1289
## S:0.5  -0.45907074        S     0.5 2.9690 4.1918 2.1030       a 0.1955
```

The plot of the Result data with 'Genotype' as a grouping factor is produced as follows. Note that the input file in the result table of the `qpcrANOVA` function.

```r
# Plot of the Result data with 'Genotype' as grouping factor
twoFACTORplot(res,
   x.axis.factor = Drought,
   group.factor = Genotype,
   width = 0.5,
   fill = "Greens",
   y.axis.adjust = 0.5,
   y.axis.by = 2,
   ylab = "Relative Expression",
   xlab = "Drought Levels",
   legend.position = c(0.09, 0.8),
   show.letters = TRUE)
```

<img src="https://raw.githubusercontent.com/mirzaghaderi/iqpcr/7ff54cdbd095c731ba12eaa4a6aa474ef2b039eb/Figure%202.jpg">

*Figure 2 - Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests via ‘qpcrTTESTplot’ function (A). Plot of average Fold changes of one gene under a three-level conditions which level1 has been selected as check. Check level can be changed by the user. The plot produced by the ‘oneFACTORfcplot’ function (B). plot of the same data of ‘B’ represented as Relative expression using ‘oneFACTORplot’ function (C). Error bars represent 95% confidence interval.*


<img src="https://raw.githubusercontent.com/mirzaghaderi/iqpcr/7ff54cdbd095c731ba12eaa4a6aa474ef2b039eb/Figure%203.jpg">

*Figure 3 - A) Plot of Fold change of three genes (differentially arranged by an argument compared to Figure 2a). B-D) Average relative expression (RE) of a target gene under two or three factors produced by ‘twoFACTORplot’ (C) and ‘twoFACTORplot’ (B and D) functions. Error bars represent standard deviations, albeit, error type can be set to confidence interval. Means (columns) lacking letters in common have significant differences at alpha = 0.05 as resulted from an ‘LSD.test’.*



```r
# Plotting the same data with 'Drought' as grouping factor
twoFACTORplot(res,
   x.axis.factor = Genotype,
   group.factor = Drought,
   xlab = "Genotype",
   fill = "Blues",
   show.letters = FALSE)
```


### A three-factorial experiment example

```r
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_3factor_b, numberOfrefGenes = 1)$Result
res
```

```
##        Type Conc SA   wDCt    LCL    UCL letters    std
## R:H:A1    R    H A1 0.9885 1.5455 0.6323      cd 0.0979
## R:H:A2    R    H A2 1.7371 2.7159 1.1110      bc 0.1747
## R:L:A1    R    L A1 0.2852 0.4459 0.1824       f 0.0072
## R:L:A2    R    L A2 0.0641 0.1002 0.0410       g 0.0773
## R:M:A1    R    M A1 0.6271 0.9804 0.4011      de 0.3508
## R:M:A2    R    M A2 0.3150 0.4925 0.2015       f 0.1105
## S:H:A1    S    H A1 2.9690 4.6420 1.8990      ab 0.1955
## S:H:A2    S    H A2 5.1934 8.1197 3.3217       a 0.7893
## S:L:A1    S    L A1 0.7955 1.2438 0.5088       d 0.2190
## S:L:A2    S    L A2 1.5333 2.3973 0.9807       c 0.1562
## S:M:A1    S    M A1 0.4147 0.6483 0.2652      ef 0.1289
## S:M:A2    S    M A2 0.7955 1.2438 0.5088       d 0.2368
```

```r
# Arrange the first three colunms of the result table.
# This determines the columns order and shapes the plot output.
threeFACTORplot(res,
    arrangement = c(3, 1, 2),
    legend.position = c(0.9, 0.85),
    xlab = "condition")
```


```r
threeFACTORplot(res,
   arrangement = c(1, 2, 3),
   bar.width = 0.5,
   fill = "Greys",
   xlab = "Genotype",
   ylab = "Relative Expression")
```


```r
# releveling a factor levels first
res$Conc <- factor(res$Conc, levels = c("L","M","H"))
res$Type <- factor(res$Type, levels = c("S","R"))


# When using ci as error, increase y.axis.adjust to see the plot correctly!
threeFACTORplot(res,
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
   show.letters = TRUE)
```



# Checking normality of residuals

If the residuals from a `t.test` or an `lm` object are not normally distributed, the grouping letters (deduced from the `LSD.test`) might be violated. In such cases, one could apply another data transformation to the wDCt data for ANOVA and mean comparison purposes or use non-parametric tests such as the Mann-Whitney test (also known as the Wilcoxon rank-sum test), `wilcox.test()`, which is an alternative to `t.test`, or the `kruskal.test()` test which alternative to one-way analysis of variance, to test the difference between medians of the populations using independent samples. However, the `t.test` function (along with the `qpcrTTEST` function described above) includes the `var.equal` argument. When set to `FALSE`, perform `t.test` under the unequal variances hypothesis.


```r
residualsCRD <- qpcrANOVA(data_3factor_b, numberOfrefGenes = 1)$lmCRD$residuals
shapiro.test(residualsCRD) 
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residualsCRD
## W = 0.88893, p-value = 0.001724
```

```r
qqnorm(residualsCRD)
qqline(residualsCRD, col = "red")
```


# Mean of technical replicates
Calculating the mean of technical replicates and getting an output table appropriate for subsequent ANOVA analysis can be done using the `meanTech` function. For this, the input data set should follow the column arrangement of the following example data. Grouping columns need to be specified under the `groups` argument of the `meanTech` function.


```r
# See example input data frame:
data_withTechRep
```

```
factor1 factor2 factor3 biolrep techrep Etarget targetCt Eref  refCt
  Line1    Heat    Ctrl       1       1       2   33.346    2 31.520
  Line1    Heat    Ctrl       1       2       2   28.895    2 29.905
  Line1    Heat    Ctrl       1       3       2   28.893    2 29.454
  Line1    Heat    Ctrl       2       1       2   30.411    2 28.798
  Line1    Heat    Ctrl       2       2       2   33.390    2 31.574
  Line1    Heat    Ctrl       2       3       2   33.211    2 31.326
  Line1    Heat    Ctrl       3       1       2   33.845    2 31.759
  Line1    Heat    Ctrl       3       2       2   33.345    2 31.548
  Line1    Heat    Ctrl       3       3       2   32.500    2 31.477
  Line1    Heat   Treat       1       1       2   33.006    2 31.483
  Line1    Heat   Treat       1       2       2   32.588    2 31.902
  Line1    Heat   Treat       1       3       2   33.370    2 31.196
  Line1    Heat   Treat       2       1       2   36.820    2 31.440
  Line1    Heat   Treat       2       2       2   32.750    2 31.300
  Line1    Heat   Treat       2       3       2   32.450    2 32.597
  Line1    Heat   Treat       3       1       2   35.238    2 31.461
  Line1    Heat   Treat       3       2       2   28.532    2 30.651
  Line1    Heat   Treat       3       3       2   28.285    2 30.745
      .       .       .       .       .       .        .           .
```

Calculating the mean of technical replicates using the `meanTech` function:
```r
meanTech(data_withTechRep, groups = 1:4)
```

```
## # A tibble: 6 × 8
## # Groups:   factor1, factor2, factor3 [2]
##   factor1 factor2 factor3 biolrep Etarget targetCt  Eref refCt
##   <chr>   <chr>   <chr>     <int>   <dbl>    <dbl> <dbl> <dbl>
## 1 Line1   Heat    Ctrl          1       2     30.4     2  30.3
## 2 Line1   Heat    Ctrl          2       2     32.3     2  30.6
## 3 Line1   Heat    Ctrl          3       2     33.2     2  31.6
## 4 Line1   Heat    Treat         1       2     33.0     2  31.5
## 5 Line1   Heat    Treat         2       2     34.0     2  31.8
## 6 Line1   Heat    Treat         3       2     30.7     2  31.0
```

# Citation

```r
citation("iqpcr")
```


# Contact 
Email: gh.mirzaghaderi@uok.ac.ir


# References
Livak, Kenneth J, and Thomas D Schmittgen. 2001. Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR and the Double Delta CT Method. Methods 25 (4). <a href="https://doi.org/10.1006/meth.2001.1262">doi.org/10.1006/meth.2001.1262</a>.


Ganger, MT, Dietz GD, and Ewing SJ. 2017. A common base method for analysis of qPCR data and the application of simple blocking in qPCR experiments. BMC bioinformatics 18, 1-11. <a href="https://doi.org/10.1186/s12859-017-1949-5">doi.org/10.1186/s12859-017-1949-5</a>.


Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006. Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85). <a href="https://doi.org/10.1186/1471-2105-7-85">doi.org/10.1186/1471-2105-7-85</a>.


.
