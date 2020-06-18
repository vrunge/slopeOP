<!-- 
%\VignetteEngine{knitr::rmarkdown} 
%\VignetteIndexEntry{An Introduction to slopeOP}
--> 


[![Build Status](https://travis-ci.com/vrunge/slopeOP.svg?branch=master)](https://travis-ci.com/vrunge/slopeOP)
[![](https://img.shields.io/badge/docs-vignettes-blue.svg)](https://github.com/vrunge/slopeOP)


# slopeOP Vignette
### Vincent Runge
#### LaMME, Evry University
### June 13, 2019

> [Introduction](#intro)

> [The slopeOP function](#sf)

> [Options for constraining inference](#options)

> [The slopeSN function](#sn)

> [plot function](#plot)

> [Python bindings](#python)

<a id="intro"></a>

## Introduction

The package `slopeOP` is designed to segment univariate data <img src="/tex/b704a970e46e8418f3ef56718438b122.svg?invert_in_darkmode&sanitize=true" align=middle width=126.38913869999998pt height=24.65753399999998pt/> by a continuous piecewise linear signal with restrictions on starting/ending values for the inferred segments. The finite set of states <img src="/tex/cef39aeb23a61b09d838693a0897fe03.svg?invert_in_darkmode&sanitize=true" align=middle width=11.187179849999989pt height=22.465723500000017pt/> contains these values. 

When we write <img src="/tex/dca84180777f523a6d1cb11ceff47536.svg?invert_in_darkmode&sanitize=true" align=middle width=124.85777534999997pt height=14.15524440000002pt/>, the variable <img src="/tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode&sanitize=true" align=middle width=7.7054801999999905pt height=14.15524440000002pt/> goes through all the values of <img src="/tex/cef39aeb23a61b09d838693a0897fe03.svg?invert_in_darkmode&sanitize=true" align=middle width=11.187179849999989pt height=22.465723500000017pt/> from the smallest one to the biggest one. For computational efficiency we recommend to have <img src="/tex/13f705a91487b0bb06669ee7c6c3552f.svg?invert_in_darkmode&sanitize=true" align=middle width=105.80650409999998pt height=22.831056599999986pt/> but this is not mandatory. 

The cost for data <img src="/tex/7d54d6946c15cd803b31e7da5077c813.svg?invert_in_darkmode&sanitize=true" align=middle width=40.78280084999999pt height=14.15524440000002pt/>, <img src="/tex/a7e59809c70c5654d0732669ce4d3cf6.svg?invert_in_darkmode&sanitize=true" align=middle width=36.90056204999999pt height=20.221802699999984pt/>, with linear interpolation from value <img src="/tex/286f7d4815c0996530bda7973b1ec5ea.svg?invert_in_darkmode&sanitize=true" align=middle width=14.25802619999999pt height=14.15524440000002pt/> to value <img src="/tex/97c7f491f7ac1623c0a86b1fb656029b.svg?invert_in_darkmode&sanitize=true" align=middle width=14.25802619999999pt height=14.15524440000002pt/> is given by

<p align="center"><img src="/tex/32da28b15640aa1c521680c305ded2fc.svg?invert_in_darkmode&sanitize=true" align=middle width=401.57801639999997pt height=48.39056475pt/></p>

The value <img src="/tex/286f7d4815c0996530bda7973b1ec5ea.svg?invert_in_darkmode&sanitize=true" align=middle width=14.25802619999999pt height=14.15524440000002pt/> is "unseen" as the cost <img src="/tex/6598c6a71c20e965dba8e739fd65d865.svg?invert_in_darkmode&sanitize=true" align=middle width=70.78262729999999pt height=26.76175259999998pt/> obtained at index <img src="/tex/0fe1677705e987cac4f589ed600aa6b3.svg?invert_in_darkmode&sanitize=true" align=middle width=9.046852649999991pt height=14.15524440000002pt/> is not present in the summation.

Data are generated by the model 

<p align="center"><img src="/tex/8525275768014f9e64d0244f61ab600f.svg?invert_in_darkmode&sanitize=true" align=middle width=481.47380489999995pt height=35.82121785pt/></p>

with <img src="/tex/9e2fef61c286c438282ce2bf9513dd8d.svg?invert_in_darkmode&sanitize=true" align=middle width=239.6036841pt height=21.18721440000001pt/>, <img src="/tex/e5e607c35cb2b5fa02be9a25cbb7733b.svg?invert_in_darkmode&sanitize=true" align=middle width=107.10602924999999pt height=22.465723500000017pt/> and <img src="/tex/5118258da41c2cdc7d0fd96f2955fe02.svg?invert_in_darkmode&sanitize=true" align=middle width=95.95543319999999pt height=26.76175259999998pt/> identically and independently distributed (iid). The vector <img src="/tex/aaf8528aaf165a68b1715ed210c179b4.svg?invert_in_darkmode&sanitize=true" align=middle width=118.53886274999999pt height=24.65753399999998pt/> is called a changepoint vector. The optimization problem is then the following:

<p align="center"><img src="/tex/c3f1419da1436d3dcaf4bb12df979464.svg?invert_in_darkmode&sanitize=true" align=middle width=453.2107701pt height=72.48949455pt/></p>

where the states defined inside the cost function yield the continuity constraint between successive segments.

<img src="/tex/99751e94989c68f9be0f6aa442bc80d5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.302373649999986pt height=22.831056599999986pt/> is a penalty parameter, understood as an additional cost when introducing a new segment. 

Notice that the cost can be computed in constant time with the formula 

<p align="center"><img src="/tex/65dc4200c058ac688623e4791c2e8939.svg?invert_in_darkmode&sanitize=true" align=middle width=536.11095285pt height=34.3600389pt/></p>

<p align="center"><img src="/tex/a1d9b5c990cdc52fb8027b2688f7cd9c.svg?invert_in_darkmode&sanitize=true" align=middle width=312.55418204999995pt height=39.887022449999996pt/></p>

where

<p align="center"><img src="/tex/2b92786806bccea954f56077a6b80ab1.svg?invert_in_darkmode&sanitize=true" align=middle width=512.84137905pt height=47.02068525pt/></p>

To address the continuity constraint by a dynamic programming algorithm, we introduce the function <img src="/tex/2f0bbd76007598852957bb89252ea7cb.svg?invert_in_darkmode&sanitize=true" align=middle width=72.55010564999998pt height=24.65753399999998pt/> which is the optimal penalized cost up to position <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> with a last infered value equal to <img src="/tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode&sanitize=true" align=middle width=7.7054801999999905pt height=14.15524440000002pt/> (at position t). The idea is then to update a set

<p align="center"><img src="/tex/6f0b14e5caf963602d6ad84a46571ed2.svg?invert_in_darkmode&sanitize=true" align=middle width=237.13390965pt height=16.438356pt/></p>

at any time step <img src="/tex/985dbad85d61b07e704840368824ee09.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/>. <img src="/tex/95e4b0db7bae72fed7d46a11d7f120e1.svg?invert_in_darkmode&sanitize=true" align=middle width=32.14725194999999pt height=14.15524440000002pt/> and <img src="/tex/17a6319f14ff4191293702fe0b373c37.svg?invert_in_darkmode&sanitize=true" align=middle width=33.95508434999999pt height=14.15524440000002pt/> are the bounds of the interval of possible ending values for the considered data to segment. They can be determined in a preprocessing step.

The new update with continuity constraint takes the form

<p align="center"><img src="/tex/a699538b2ec9cb612fe187e38f41a204.svg?invert_in_darkmode&sanitize=true" align=middle width=432.32085104999993pt height=39.452455349999994pt/></p>

where the presence of the same value <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/> in <img src="/tex/e05578c58a4f9e4f1301d4cf24e3234c.svg?invert_in_darkmode&sanitize=true" align=middle width=20.387619449999992pt height=22.465723500000017pt/> and the cost realizes the continuity constraint. At initial step we simply have <img src="/tex/cf8c7572e6ac52e95e69475a2b4281e3.svg?invert_in_darkmode&sanitize=true" align=middle width=86.58177495pt height=24.65753399999998pt/>. 

The slopeOP function computes <img src="/tex/8dfc04a9fbcadb584faf331368654540.svg?invert_in_darkmode&sanitize=true" align=middle width=39.274025999999985pt height=24.65753399999998pt/> for all <img src="/tex/d4982b4fcd5203e3c32cdf5c3ba3211f.svg?invert_in_darkmode&sanitize=true" align=middle width=38.98379759999999pt height=22.465723500000017pt/> and <img src="/tex/28166798932f09dc74160bd665b9663e.svg?invert_in_darkmode&sanitize=true" align=middle width=74.25025244999999pt height=21.18721440000001pt/>. The argminimum state into the set <img src="/tex/44314730b45a8f895f72c002bb5251d5.svg?invert_in_darkmode&sanitize=true" align=middle width=21.550742399999987pt height=22.465723500000017pt/> gives the last value of the last inferred segment. A backtracking procedure eventually returns the optimal changepoint vector with all its associated state values.

<a id="sf"></a>

## The slopeOP function

We install the package from Github:


```r
#devtools::install_github("vrunge/slopeOP")
library(slopeOP)
```

We simulate data with the function `slopeData` with arguments `index` (a changepoint vector), `states` its associated state values and the `noise` level which is the standard deviation of a normal standard noise (iid).

```r
data <- slopeData(index = c(1,100,200,300,500), states = c(0,1,0,3,2), noise = 1)
```

The changepoint detection is achieved by using the function `slopeOP`


```r
slopeOP(data = data, states = c(0,1,2,3), penalty = 10)
```

```
## $changepoints
## [1]   1 101 201 307 500
## 
## parameters
## [1] 0 1 0 3 2
## 
## globalCost
## [1] 496.4254
## 
## attr(,"class")
## [1] "slopeOP"
```

In `slopeOP` function, the parameter `type` is `channel` by default. With type equal to `channel` we use the monotonicity property in optimal cost matrix to reduce time complexity. If it is equal to `pruning` we prune some positions using a theorem taking into account unseen data. The pruning option is similar to PELT pruning but less effective than channel in this case.

<a id="options"></a>

## Options for constraining inference


Parameter `constraint` can be set to `isotonic` which corresponds to a restriction to nondecreasing state vectors.


```r
myData <- slopeData(index = c(1,150,200,350,500,750,1000), states = c(71,73,70,75,77,73,80), noise = 1)
slopeOP(data = myData, states = 71:80, penalty = 5, constraint = "isotonic")
```

```
## changepoints
## [1]    1   63  254  356  828 1000
## 
## parameters
## [1] 71 72 72 75 75 80
## 
## globalCost
## [1] 1768.407
## 
## attr(,"class")
## [1] "slopeOP"
```

With `constraint` equal to `unimodal` the infered signal is increasing and then decreasing.


```r
myData <- slopeData(index = c(1,150,200,350,500,750,1000), states = c(71,73,70,75,78,73,75), noise = 1)
slopeOP(data = myData, states = 71:80, penalty = 5, constraint = "unimodal")
```

```
## changepoints
## [1]    1  316  317  502  697 1000
## 
## parameters
## [1] 71 73 74 78 74 74
## 
## globalCost
## [1] 1407.936
## 
## attr(,"class")
## [1] "slopeOP"
```

We also can limit the angles between successive segments with `constraint` equal to `smoothing` and the parameter `minAngle` in degree.


```r
myData <- slopeData(c(1,30,40,70,100,150,200),c(70,80,70,80,70,80,70), noise = 0.5)
slopeOP(data = myData, states = 70:80, penalty = 5, constraint = "smoothing", minAngle = 170)
```

```
## changepoints
##  [1]   1   7  19  25  28  34  40  46  47  53  62  68  73  79  94 100 103
## [18] 149 154 160 200
## 
## parameters
##  [1] 71 72 76 77 77 76 74 73 73 74 77 78 78 77 72 71 71 79 79 78 70
## 
## globalCost
## [1] 319.9376
## 
## attr(,"class")
## [1] "slopeOP"
```

<a id="sn"></a>

## The slopeSN function

With `slopeSN`, we are able to constrain to number of segments in the inference


```r
myData <- slopeData(index = c(1,10,20,30), states = c(0,5,3,6), noise = 1)
res <- slopeSN(data = myData, states = 0:6, nbSegments = 2)
res
```

```
## changepoints
## [1]  1  6 30
## 
## parameters
## [1] 0 3 5
## 
## globalCost
## [1] 35.42239
## 
## attr(,"class")
## [1] "slopeOP"
```

<a id="plot"></a>

## Plot function

A simple plot function can be used to show raw data with the inferred segments on the same graph. Option `data =` should be always present in the call of the plot function.


```r
data <- slopeData(c(1,11,21),c(70,80,70), noise = 2)
slope <- slopeOP(data, states = 70:80, penalty = 10, constraint = "null", type = "channel")
plot(slope, data = data)
```

![](Rplot.png)

<a id="python"></a>

## Python Bindings

Instruction to generate slopeOP python module. The following commands have to be run in the base directory of the slopeOP repo.

- Update git submodules

```{sh}
git submodule init
git submodule update
```

- Build and install the module. You have 2 options:
    1. `pip install -e .` to generate a python importable `".so"` module in this folder
    2. `pip install .` to install the module in the python environment.

- run python and import the module with `import slopeOP`.

The `SlopeOP` python module has 2 functions:

- `op2D(x, y, penalty)` for pice-wise linear OP
  
- `slopeOP(data, states, penality, constraint="null", minAngle=0, type="channel")` for SlopeOP

### troubleshooting
In case of error during build, try to delete the `build` folder (if any) and try to build again.

[Back to Top](#top)

