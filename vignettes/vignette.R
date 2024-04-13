## ----setup, include = FALSE, fig.align='center', warning = F, message=F-------
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(echo = TRUE)

## ----eval= T, include= F, message=FALSE, warning = FALSE----------------------
library(rtpcr)
library(agricolae)
library(dplyr)
library(reshape2)
library(tidyr)
library(lme4)
library(ggplot2)
library(grid)

## ----eval= F, include= T, message=FALSE, warning = FALSE----------------------
#  
#  # install `rtpcr` from github (under development)
#  
#  devtools::install_github("mirzaghaderi/rtpcr")
#  
#  # I strongly recommend to install the package with the vignette as it contains information about how to use the 'rtpcr' package. Through the following code, Vignette is installed as well.
#  
#  devtools::install_github("mirzaghaderi/rtpcr", build_vignettes = TRUE)

## ----eval= T------------------------------------------------------------------
data_efficiency

## ----eval = T, fig.height = 3, fig.align = 'center', fig.cap = 'Standard curve and the amplification efficiency analysis of target and reference genes. A sample data arrangement that is required as input for the calculation of amplification efficiency by the efficiency function.', warning = FALSE, message = FALSE----
efficiency(data_efficiency)

## ----eval= T, fig.height = 3, fig.width = 5, fig.align = 'center'-------------
data_ttest

## ----eval= T------------------------------------------------------------------
qpcrTTEST(data_ttest, 
          numberOfrefGenes = 1,
          paired = F, 
          var.equal = T)

## ----eval= T, fig.height=3, fig.width=8, fig.align='center', fig.cap = "Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests via ‘qpcrTTESTplot’ function.", warning = F, message = F----

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
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

