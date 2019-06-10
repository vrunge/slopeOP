<a id="top"></a>

<!-- 
%\VignetteEngine{knitr::rmarkdown} 
%\VignetteIndexEntry{An Introduction to gfpop}
--> 

# slopeOP Vignette
### Vincent Runge
#### LaMME, Evry University
### June 10, 2019

> [Introduction](#intro)

> [Quick Start](#mf)

> [Some examples](#se)

<a id="intro"></a>

## Introduction

<a id="mf"></a>

## The main function

We install the package from Github:

```r
#devtools::install_github("vrunge/slopeOP")
library(slopeOP)
```

We simulate data with the function `slopeData`

```r
data <- slopeData(c(1,100,200,300,500),c(0,1,0,3,2),1)
```
And the changepoint detection is done by using the function `slopeOP`


```r
slopeOP(c(1,2,3,2,1,2,3,2,1,2,3),c(1,2,3),0.5)
```

```
## $changepoints
## [1]  1  3  5  7  9 11
## 
## $parameters
## [1] 1 3 1 3 1 3
## 
## $globalCost
## [1] 2.5
```


<a id="se"></a>

## Some examples

To be done ;) 

[Back to Top](#top)
