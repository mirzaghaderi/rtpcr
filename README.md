# iqpcr 

iqpcr is a package for statistical analysis of real-time PCR data in R.

# Overview

iqpcr is a package for statistical analysis and graphical presentation of real-time PCR (quantitative PCR  or qPCR) data. Various analysis methods are employed for real-time PCR data to measure the mRNA levels under different experimental conditions. iqpcr handles amplification efficiency calculation and statistical analysis of real-time PCR data based on one reference gene. By accounting for amplification efficiency values, iqpcr was developed using a general calculation method described by <a href="https://doi.org/10.1186/s12859-017-1949-5">Ganger et al. 2017</a>., covering both the Livak and Pfaffl methods. Based on the experimental conditions, the functions of the iqpcr package use t-test (for experiments with a two-level factor) or analysis of variance (for cases where more than two levels or factors exist) to calculate the fold change or relative expression. The functions further provide standard deviations and confidence limits for means, apply statistical mean comparisons, and present letter mean grouping. To facilitate functions application, different dataset have been used in the examples and the outputs are explained. An outstanding feature of iqpcr package is providing publication-ready bar plots with various controlling arguments for experiments with up to three different factors. The iqpcr package is user-friendly and easy to work with and provides an applicable resource for analyzing real-time PCR data in R.

# Calculation methods
Among the various approaches developed for data analysis in real-time PCR, the Livak method, also known as the $2^{-\Delta\Delta C_t}$ method, stands out for its simplicity and widespread use. 


$$Fold change = 2^{-ΔΔCt}$$



Fold change = 2^(-ΔΔCt)

= 2^(-(Ct_target-Ct_ref)_Tr) / 2^(-(Ct_target-Ct_ref)_Co)

= 2^[-(Ct_target-Ct_ref)_Tr - (Ct_target-Ct_ref)_Co]

= 2^[-(ΔCt_Tr - ΔCt_Co)]



where Tr is Treatment and Co is Control conditions, respectively. This method assumes that both the target and reference genes are amplified with efficiencies close to 100%, allowing for the relative quantification of gene expression levels.

On the other hand, the Pfaffl method offers a more flexible approach by accounting for differences in amplification efficiencies between the target and reference genes. This method adjusts the calculated expression ratio by incorporating the specific amplification efficiencies, thus providing a more accurate representation of the relative gene expression levels.

Fold change = E^(-(Ct_Tr-Ct_Co)_target) / E^(-(Ct_Tr-Ct_Co)_ref)

# A generalized calculation method
the iqpcr packager  calculates wDCt values for each experimental condition based on the efficiency and Ct values as described by Ganger et al. (2017).


$$wΔCt = log_{10}(E_{target}).Ct_{target} - log_{10}(E_{ref}).Ct_{ref}$$

if wDCt values are calculated using the efficiencies (E) values, calculations match the formula of Pfaffl. if 2 (100% amplification efficiency) is used, calculations match the formula of Livak. Analysis of variance (ANOVA) and t-test (paired or unpaired) is done on the wDCt values. Subsequently, the mean of wCt over each condition along with the resulting statistics (e.g. confidence limits for the means) is converted by a power of ten conversion for presentation. For example, relative expression is calculated as:
$$\text{Relative Expression} = 10^{-\overline{w\Delta Ct}}$$
When there are only two conditions (Tr and Co), fold change can also be calculated:

$$Fold Change = 10^{-(\overline{w\Delta Ct}{Tr}-\overline{w\Delta Ct}{Co})}$$


# Installing and loading



<!-- `iqpcr` is a CRAN package and can be installed through `install.packages()`.
install.packages("iqpcr")
library(iqpcr) -->

The development version of the `iqpcr` package can be obtained through: 

```r
# Install `iqpcr` from github (under development)
devtools::install_github("mirzaghaderi/iqpcr", build_vignettes = TRUE)

# Loading the package
library(iqpcr)
```

# Amplification efficiency data analysis
## Sample data of amplification efficiency

In order to calculate the amplification efficiencies
of a target and a reference gene, a dataframe should be prepared with 3 columns 
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
The following `efficiency` function calculates amplification efficiency of a target and a
reference gene, and present the related standard curves along with the Slope, Efficiency, and R2 statistics. The function also compares the slopes of the two standard curves. For this, a regression line is fitted using the $\Delta C_t$ values. If $2^{-\Delta\Delta C_t}$ method is intended, the slope should not exceed 0.2!


```r
efficiency(data_efficiency)
```

![alt text](https://raw.githubusercontent.com/mirzaghaderi/iqpcr/721b5879bd0ac967b0666e1b40ff4c5ae23b4a3f/Figure1.png)
*Standard curve and the amplification efficiency analysis of target and reference genes.*


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
when a target gene is assessed under two different conditions (for example Control and treatment), it is possible to calculate the average fold change expression i.e. $2^{-\Delta \Delta C_t}$ of the target gene in treatment relative to control conditions. For this, the data should be prepared according to the following dataset consisting of 4 columns belonging to condition levels, E (efficiency), genes and Ct values, respectively. Each Ct value is the mean of technical replicates. Complete amplification efficiencies of 2 have been assumed here for all wells but the calculated efficiencies can be used instead. 


```r
data_ttest
```

```
##    Condition E    Gene    Ct
## 1    control 2 C2H2-26 31.26
## 2    control 2 C2H2-26 31.01
## 3    control 2 C2H2-26 30.97
## 4  treatment 2 C2H2-26 32.65
## 5  treatment 2 C2H2-26 32.03
## 6  treatment 2 C2H2-26 32.40
## 7    control 2 C2H2-01 31.06
## 8    control 2 C2H2-01 30.41
## 9    control 2 C2H2-01 30.97
## 10 treatment 2 C2H2-01 28.85
## 11 treatment 2 C2H2-01 28.93
## 12 treatment 2 C2H2-01 28.90
## 13   control 2 C2H2-12 28.50
## 14   control 2 C2H2-12 28.40
## 15   control 2 C2H2-12 28.80
## 16 treatment 2 C2H2-12 27.90
## 17 treatment 2 C2H2-12 28.00
## 18 treatment 2 C2H2-12 27.90
## 19   control 2     ref 28.87
## 20   control 2     ref 28.42
## 21   control 2     ref 28.53
## 22 treatment 2     ref 28.31
## 23 treatment 2     ref 29.14
## 24 treatment 2     ref 28.63
```

### Data analysis under two conditions

Here, the above data set was used for the Fold Change expression analysis of the target genes using the `qpcrTTEST` function. This function performs a t.test based analysis to any number of genes that 
have been evaluated under control and treatment conditions. The analysis can be done for unpaired or paired condition. The output is a table of target gene names, fold changes confidence limits and the t.test derived p-values. The `qpcrTTEST` function include the `var.equal` argument. When set to `FALSE`,
`t.test` is performed under unequal variances hypothesis.


```r
qpcrTTEST(data_ttest, paired = F, var.equal = T)
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
`qpcrTTESTplot` function generates a bar plot of Fold Changes and confidence intervals for the target genes. the `qpcrTTESTplot` function accepts any gene name and any replicates. `qpcrTTESTplot` function automatically put appropriate signs of **, * on top of the plot columns based on the output p-values.



```r
# Producing the plot
qpcrTTESTplot(data_ttest)
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-8-1.png" alt="Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests. Error bars represent 95% confidence interval."  />
<p class="caption">Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests. Error bars represent 95% confidence interval.</p>
</div>

```r
# Producing the plot: specifying gene order
qpcrTTESTplot(data_ttest,
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

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-8-2.png" alt="Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests. Error bars represent 95% confidence interval."  />
<p class="caption">Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests. Error bars represent 95% confidence interval.</p>
</div>

## A target gene under more than two conditions (ANOVA)

Analysis of variance (ANOVA) of factorial experiments in the frame of a completely randomized design (CRD) can be done by the `qpcrANOVA` function. ANOVA of qPCR data is suitable when there is a factor with more than two levels, or when there are more than a experimental factor. The input data set should be prepared as shown below. Factor columns should be presented first followed by biological replicates and efficiency and Ct values of target and reference genes. The example data set below (`data_3factor_a`) represents amplification efficiency and Ct values for target and reference genes under three grouping factors (two different cultivars, three drought levels and the presence or absence of bacteria). The table can contain any number of factor columns. The factor columns, should be followed by five other columns assigned to biological replicates (r), efficiency of the target gene, Ct values of target gene, efficiency of reference gene, and Ct values of the reference gene, respectively. Here, the efficiency of 2 has been used for all wells, but the calculated efficiencies can be used instead.



```r
# See a sample dataset
data_3factor_a
```

```
##    Genotype Drought SA Rep   EPO  POCt EGAPDH GAPDHCt
## 1         R    0.00 A1   1 1.839 33.30  1.918   31.53
## 2         R    0.00 A1   2 1.839 33.39  1.918   31.57
## 3         R    0.00 A1   3 1.839 33.34  1.918   31.50
## 4         R    0.00 A2   1 1.839 34.01  1.918   31.48
## 5         R    0.00 A2   2 1.839 36.82  1.918   31.44
## 6         R    0.00 A2   3 1.839 35.44  1.918   31.46
## 7         R    0.25 A1   1 1.839 32.73  1.918   31.30
## 8         R    0.25 A1   2 1.839 32.46  1.918   32.55
## 9         R    0.25 A1   3 1.839 32.60  1.918   31.92
## 10        R    0.25 A2   1 1.839 33.37  1.918   31.19
## 11        R    0.25 A2   2 1.839 33.12  1.918   31.94
## 12        R    0.25 A2   3 1.839 33.21  1.918   31.57
## 13        R    0.50 A1   1 1.839 33.48  1.918   33.30
## 14        R    0.50 A1   2 1.839 33.27  1.918   33.37
## 15        R    0.50 A1   3 1.839 33.32  1.918   33.35
## 16        R    0.50 A2   1 1.839 32.53  1.918   33.47
## 17        R    0.50 A2   2 1.839 32.61  1.918   33.26
## 18        R    0.50 A2   3 1.839 32.56  1.918   33.36
## 19        S    0.00 A1   1 1.839 26.85  1.918   26.94
## 20        S    0.00 A1   2 1.839 28.17  1.918   27.69
## 21        S    0.00 A1   3 1.839 27.99  1.918   27.39
## 22        S    0.00 A2   1 1.839 28.71  1.918   29.45
## 23        S    0.00 A2   2 1.839 29.01  1.918   29.46
## 24        S    0.00 A2   3 1.839 28.82  1.918   29.48
## 25        S    0.25 A1   1 1.839 30.41  1.918   28.70
## 26        S    0.25 A1   2 1.839 29.49  1.918   28.66
## 27        S    0.25 A1   3 1.839 29.98  1.918   28.71
## 28        S    0.25 A2   1 1.839 28.91  1.918   28.09
## 29        S    0.25 A2   2 1.839 28.60  1.918   28.65
## 30        S    0.25 A2   3 1.839 28.59  1.918   28.37
## 31        S    0.50 A1   1 1.839 29.03  1.918   30.61
## 32        S    0.50 A1   2 1.839 28.73  1.918   30.20
## 33        S    0.50 A1   3 1.839 28.83  1.918   30.49
## 34        S    0.50 A2   1 1.839 28.29  1.918   30.84
## 35        S    0.50 A2   2 1.839 28.53  1.918   30.65
## 36        S    0.50 A2   3 1.839 28.28  1.918   30.74
```

`qpcrANOVA` function performs  anova based on both factorial arrangement and completely randomized design (CRD). For the latter, a column of treatment combinations is made first as grouping factor followed by ANOVA. You can call the input data set along with the added wCt and treatment combinations by `qpcrANOVA(data_3factor_a)$Final_data`. CRD based analysis is especially useful when post-hoc test and mean comparisons/grouping is desired for all treatment combinations. The final results along with the ANOVA tables can be called by `qpcrANOVA(data_3factor_a)`.

### Reverse ordering of the grouping letters

One may be interested to present the statistical mean comparison result in the frame of grouping letters. This is rather challenging because in the grouping output of mean comparisons (via `LSD.test` function of agricolae package), means are sorted into descending order so that the largest mean, is the first in the
table and "a" letter is assigned to it. If `LSD.test` is applied to the wCt means,
the biggest wCt mean receive "a" letter as expected, but this value turns into
the smallest mean after its reverse log transformation by  $10^{-(\Delta Ct)}$. to solve
this issue, I used a function that assigns the grouping letters appropriately.

### Output table of the analysis

The `qpcrANOVA` function produces the main analysis output 
including mean wDCt, LCL, UCL, grouping letters and standard deviations. The
standard deviation for each mean is derived from the back transformed raw wDCt values from biological
replicates for that mean.


```r
# If the data include technical replicates, means of technical replicates
# should be calculated first using meanTech function.

# Applying ANOVA analysis
qpcrANOVA(data_3factor_a, p.adj = "none")
```

```
## $Final_data
##    Genotype Drought SA rep Etarget Cttarget  Eref Ctref         wDCt
## 1         R    0.00 A1   1   1.839    33.30 1.918 31.53 -0.107644864
## 2         R    0.00 A1   2   1.839    33.39 1.918 31.57 -0.095146452
## 3         R    0.00 A1   3   1.839    33.34 1.918 31.50 -0.088576136
## 4         R    0.00 A2   1   1.839    34.01 1.918 31.48  0.094350594
## 5         R    0.00 A2   2   1.839    36.82 1.918 31.44  0.849139197
## 6         R    0.00 A2   3   1.839    35.44 1.918 31.46  0.478359439
## 7         R    0.25 A1   1   1.839    32.73 1.918 31.30 -0.193401271
## 8         R    0.25 A1   2   1.839    32.46 1.918 32.55 -0.618399091
## 9         R    0.25 A1   3   1.839    32.60 1.918 31.92 -0.403163029
## 10        R    0.25 A2   1   1.839    33.37 1.918 31.19  0.007044382
## 11        R    0.25 A2   2   1.839    33.12 1.918 31.94 -0.271237502
## 12        R    0.25 A2   3   1.839    33.21 1.918 31.57 -0.142771163
## 13        R    0.50 A1   1   1.839    33.48 1.918 33.30 -0.560662180
## 14        R    0.50 A1   2   1.839    33.27 1.918 33.37 -0.636023745
## 15        R    0.50 A1   3   1.839    33.32 1.918 33.35 -0.617137686
## 16        R    0.50 A2   1   1.839    32.53 1.918 33.47 -0.860099085
## 17        R    0.50 A2   2   1.839    32.61 1.918 33.26 -0.779534340
## 18        R    0.50 A2   3   1.839    32.56 1.918 33.36 -0.821048287
## 19        S    0.00 A1   1   1.839    26.85 1.918 26.94 -0.515921930
## 20        S    0.00 A1   2   1.839    28.17 1.918 27.69 -0.378810500
## 21        S    0.00 A1   3   1.839    27.99 1.918 27.39 -0.341580630
## 22        S    0.00 A2   1   1.839    28.71 1.918 29.45 -0.733749907
## 23        S    0.00 A2   2   1.839    29.01 1.918 29.46 -0.657203874
## 24        S    0.00 A2   3   1.839    28.82 1.918 29.48 -0.713131375
## 25        S    0.25 A1   1   1.839    30.41 1.918 28.70 -0.071824515
## 26        S    0.25 A1   2   1.839    29.49 1.918 28.66 -0.303925762
## 27        S    0.25 A1   3   1.839    29.98 1.918 28.71 -0.188423145
## 28        S    0.25 A2   1   1.839    28.91 1.918 28.09 -0.296159461
## 29        S    0.25 A2   2   1.839    28.60 1.918 28.65 -0.536575015
## 30        S    0.25 A2   3   1.839    28.59 1.918 28.37 -0.460023224
## 31        S    0.50 A1   1   1.839    29.03 1.918 30.61 -0.977188133
## 32        S    0.50 A1   2   1.839    28.73 1.918 30.20 -0.940594725
## 33        S    0.50 A1   3   1.839    28.83 1.918 30.49 -0.996162646
## 34        S    0.50 A2   1   1.839    28.29 1.918 30.84 -1.238033791
## 35        S    0.50 A2   2   1.839    28.53 1.918 30.65 -1.120792942
## 36        S    0.50 A2   3   1.839    28.28 1.918 30.74 -1.212394748
## 
## $lmCRD
## 
## Call:
## stats::lm(formula = wDCt ~ T, data = x)
## 
## Coefficients:
## (Intercept)   TR:0.25:A2    TR:0.5:A1    TR:0.5:A2      TR:0:A1      TR:0:A2   TS:0.25:A1   TS:0.25:A2    TS:0.5:A1  
##   -0.404988     0.269333    -0.199620    -0.415239     0.307865     0.878938     0.216930    -0.025931    -0.566327  
##   TS:0.5:A2      TS:0:A1      TS:0:A2  
##   -0.785419    -0.007117    -0.296374  
## 
## 
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
## R:0.25:A2 - S:0.25:A2  0.5067 0.0202       *  0.2883  0.8905
## R:0.25:A2 - S:0.5:A1   0.1460 0.0000     ***  0.0831  0.2566
## R:0.25:A2 - S:0.5:A2   0.0882 0.0000     ***  0.0502  0.1549
## R:0.25:A2 - S:0:A1     0.5291 0.0285       *  0.3011  0.9299
## R:0.25:A2 - S:0:A2     0.2718 0.0001     ***  0.1547  0.4777
## R:0.5:A1 - R:0.5:A2    0.6087 0.0817       .  0.3463  1.0697
## R:0.5:A1 - R:0:A1      3.2173 0.0003     ***  1.8306  5.6541
## R:0.5:A1 - R:0:A2     11.9828 0.0000     ***  6.8183 21.0591
## R:0.5:A1 - S:0.25:A1   2.6095 0.0018      **  1.4848  4.5860
## R:0.5:A1 - S:0.25:A2   1.4917 0.1562          0.8488  2.6216
## R:0.5:A1 - S:0.5:A1    0.4298 0.0050      **  0.2446  0.7554
## R:0.5:A1 - S:0.5:A2    0.2595 0.0000     ***  0.1477  0.4561
## R:0.5:A1 - S:0:A1      1.5578 0.1178          0.8864  2.7377
## R:0.5:A1 - S:0:A2      0.8003 0.4228          0.4554  1.4065
## R:0.5:A2 - R:0:A1      5.2857 0.0000     ***  3.0076  9.2894
## R:0.5:A2 - R:0:A2     19.6869 0.0000     *** 11.2020 34.5986
## R:0.5:A2 - S:0.25:A1   4.2872 0.0000     ***  2.4394  7.5344
## R:0.5:A2 - S:0.25:A2   2.4508 0.0032      **  1.3945  4.3071
## R:0.5:A2 - S:0.5:A1    0.7062 0.2151          0.4018  1.2411
## R:0.5:A2 - S:0.5:A2    0.4264 0.0047      **  0.2426  0.7494
## R:0.5:A2 - S:0:A1      2.5593 0.0021      **  1.4563  4.4978
## R:0.5:A2 - S:0:A2      1.3148 0.3264          0.7481  2.3107
## R:0:A1 - R:0:A2        3.7245 0.0001     ***  2.1193  6.5457
## R:0:A1 - S:0.25:A1     0.8111 0.4509          0.4615  1.4254
## R:0:A1 - S:0.25:A2     0.4637 0.0096      **  0.2638  0.8149
## R:0:A1 - S:0.5:A1      0.1336 0.0000     ***  0.0760  0.2348
## R:0:A1 - S:0.5:A2      0.0807 0.0000     ***  0.0459  0.1418
## R:0:A1 - S:0:A1        0.4842 0.0139       *  0.2755  0.8509
## R:0:A1 - S:0:A2        0.2487 0.0000     ***  0.1415  0.4372
## R:0:A2 - S:0.25:A1     0.2178 0.0000     ***  0.1239  0.3827
## R:0:A2 - S:0.25:A2     0.1245 0.0000     ***  0.0708  0.2188
## R:0:A2 - S:0.5:A1      0.0359 0.0000     ***  0.0204  0.0630
## R:0:A2 - S:0.5:A2      0.0217 0.0000     ***  0.0123  0.0381
## R:0:A2 - S:0:A1        0.1300 0.0000     ***  0.0740  0.2285
## R:0:A2 - S:0:A2        0.0668 0.0000     ***  0.0380  0.1174
## S:0.25:A1 - S:0.25:A2  0.5717 0.0518       .  0.3253  1.0047
## S:0.25:A1 - S:0.5:A1   0.1647 0.0000     ***  0.0937  0.2895
## S:0.25:A1 - S:0.5:A2   0.0995 0.0000     ***  0.0566  0.1748
## S:0.25:A1 - S:0:A1     0.5970 0.0711       .  0.3397  1.0491
## S:0.25:A1 - S:0:A2     0.3067 0.0002     ***  0.1745  0.5390
## S:0.25:A2 - S:0.5:A1   0.2881 0.0001     ***  0.1640  0.5064
## S:0.25:A2 - S:0.5:A2   0.1740 0.0000     ***  0.0990  0.3058
## S:0.25:A2 - S:0:A1     1.0443 0.8753          0.5942  1.8353
## S:0.25:A2 - S:0:A2     0.5365 0.0318       *  0.3053  0.9428
## S:0.5:A1 - S:0.5:A2    0.6038 0.0772       .  0.3436  1.0612
## S:0.5:A1 - S:0:A1      3.6242 0.0001     ***  2.0622  6.3693
## S:0.5:A1 - S:0:A2      1.8619 0.0321       *  1.0594  3.2722
## S:0.5:A2 - S:0:A1      6.0021 0.0000     ***  3.4152 10.5483
## S:0.5:A2 - S:0:A2      3.0835 0.0004     ***  1.7545  5.4191
## S:0:A1 - S:0:A2        0.5137 0.0226       *  0.2923  0.9029
```

### Barplot with the (1-alpha)% confidence interval as error bars


```r
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_2factor)$Result
res
```

```
##                Q75 Genotype Drought   wDCt    LCL    UCL letters    std
## R:0     0.55088489        R       0 0.2852 0.4026 0.2020       d 0.0072
## R:0.25  0.31758665        R    0.25 0.6271 0.8853 0.4441      bc 0.3508
## R:0.5   0.02257725        R     0.5 0.9885 1.3956 0.7002       b 0.0979
## S:0     0.16255620        S       0 0.7955 1.1232 0.5635       b 0.2190
## S:0.25  0.44853469        S    0.25 0.4147 0.5854 0.2937      cd 0.1289
## S:0.5  -0.45907074        S     0.5 2.9690 4.1918 2.1030       a 0.1955
```

```r
# Plot of the 'res' data with 'Genotype' as grouping factor
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

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-11-1.png" alt="Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package."  />
<p class="caption">Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package.</p>
</div>

```r
# Plotting the same data with 'Drought' as grouping factor
twoFACTORplot(res,
   x.axis.factor = Genotype,
   group.factor = Drought,
   xlab = "Genotype",
   fill = "Blues",
   show.letters = FALSE)
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-11-2.png" alt="Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package."  />
<p class="caption">Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package.</p>
</div>


### A three-factorial experiment example

```r
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_3factor_b)$Result
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

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-12-1.png" alt="Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package."  />
<p class="caption">Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package.</p>
</div>

```r
threeFACTORplot(res,
   arrangement = c(1, 2, 3),
   bar.width = 0.5,
   fill = "Greys",
   xlab = "Genotype",
   ylab = "Relative Expression")
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-12-2.png" alt="Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package."  />
<p class="caption">Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package.</p>
</div>

```r
# releveling a factor levels first
res$Conc <- factor(res$Conc, levels = c("L","M","H"))
res$Type <- factor(res$Type, levels = c("S","R"))

# Producing the plot
threeFACTORplot(res,
   arrangement = c(2, 3, 1),
   bar.width = 0.5,
   fill = "Reds",
   xlab = "Drought",
   ylab = "Relative Expression",
   errorbar = "std",
   legend.title = "Genotype",
   legend.position = c(0.2, 0.8))
```

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-12-3.png" alt="Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package."  />
<p class="caption">Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package.</p>
</div>

```r
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

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-12-4.png" alt="Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package."  />
<p class="caption">Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package.</p>
</div>



# Checking normality of residuals

If the residuals from a `t.test` or a `lm` object are not normaly distributed, the grouping letters (deduced from the `LSD.test`) might be violated. In such cases, one
could apply another data transformation to the wDCt data for ANOVA and mean comparison purposes or use non-parametric tests such as the Mann-Whitney test (also known as the Wilcoxon rank-sum test), `wilcox.test()`, which is alternative to `t.test`, or the
`kruskal.test()` test which alternative to one-way analysis of variance, to test
the difference between medians of the populations using independent samples.
However, the `t.test` function (along with the `qpcrTTEST` function described
above) include the `var.equal` argument. When set to `FALSE`, performs
`t.test` under unequal variances hypothesis.


```r
residualsCRD <- qpcrANOVA(data_3factor_b)$lmCRD$residuals
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

<div class="figure" style="text-align: center">
<img src="vignette_files/figure-html/unnamed-chunk-13-1.png" alt="QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."  />
<p class="caption">QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions.</p>
</div>

# Mean of technical replicates
Calculating the mean of technical replicates and getting an output table appropriate for subsequent ANOVA analysis can be done using `meanTech` function. For this, the input data set should follow column arrangement of the following example data. Grouping columns needs to be specified under the `groups` argument of the `meanTech` function.


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
