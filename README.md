# ConditionalScoreDecomp

This repository contains R code to reproduce the results presented in the paper  

> Allen, S., Ferro, C.A.T. and Kwasniok, F. (2023). 
> A conditional decomposition of proper scores: quantifying the sources of information in a forecast.
> Quarterly Journal of the Royal Meteorological Society.

## Conditional score decompositions

This paper concerns decompositions of proper scores into terms that quantify the forecast uncertainty, resolution, and reliability. 

Proper scoring rules provide a well-established framework with which to assess probabilistic forecasts. While proper scoring rules return a single measure of forecast accuracy, it is well-known that the expected score can be decomposed into terms that quantify the uncertainty, resolution, and reliability of the forecasts. This work extends this by introducing decompositions of proper scores into uncertainty, resolution, and reliability components conditional on certain events, or states, having occurred. This allows practitioners to identify the sources of conditional miscalibration in their forecasts.

This repository provides the functionality to implement these score decompositions in practice. In particular, functionality is currently available to decompose the Brier score into uncertainty, resolution, and reliability components, conditional on certain events having occurred. Functions are also available that allow the user to visualise these decomposition terms graphically. The vignette thoroughly documents the usage of these functions, and reproduces the results in the above paper.

## Data

The data used in this study is not publicly available, so the repository instead contains a 'noisy' data set that seeks to artificially mimic the features of the actual data. Hence, the results presented in the vignette do not correspond exactly to those presented in the paper. Further details regarding the data can be found in the paper and the vignette.

## Installation and development

This package has not been submitted to CRAN, and can therefore be installed in R using devtools
```r
# install.packages("devtools")
library(devtools)
install_github("sallen12/ConditionalScoreDecomp")
```
The package is still in active development, and the vignette lists several possible extensions that could be implemented. Comments, suggestions, and input are more than welcome.
