<!-- 
%\VignetteEngine{knitr::rmarkdown} 
%\VignetteIndexEntry{An Introduction to slopeOP}
--> 

# slopeOP Vignette
### Vincent Runge
#### LaMME, Evry University
### June 13, 2019

> [Introduction](#intro)

> [The slopeOP function](#sf)

> [Options for constraining inference](#options)

<a id="intro"></a>

## Introduction

The package `slopeOP` is designed to segment univariate data <img src="/tex/b704a970e46e8418f3ef56718438b122.svg?invert_in_darkmode&sanitize=true" align=middle width=126.38913869999998pt height=24.65753399999998pt/> by a continuous piecewise linear signal with restriction on starting/ending values for the segments. The changepoint vector <img src="/tex/ed535ca37308856c5f0c6d85b4d4c676.svg?invert_in_darkmode&sanitize=true" align=middle width=209.11492305pt height=27.91243950000002pt/> defines the <img src="/tex/33359de825e43daa97171e27f6558ae9.svg?invert_in_darkmode&sanitize=true" align=middle width=37.38576269999999pt height=22.831056599999986pt/> segments <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/>, <img src="/tex/59680c09e0e2d723d0bcf2005047b028.svg?invert_in_darkmode&sanitize=true" align=middle width=73.18587374999998pt height=22.831056599999986pt/> with fixed bounds <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/> and  <img src="/tex/afdb85da7c3e7b7c3d226050994dbf5f.svg?invert_in_darkmode&sanitize=true" align=middle width=63.70246739999999pt height=14.15524440000002pt/>. We use the set <img src="/tex/36dd5900c84c0ddd4a48f4858bbd6e8f.svg?invert_in_darkmode&sanitize=true" align=middle width=264.26770754999995pt height=27.91243950000002pt/> to define the minimal global cost given by

<p align="center"><img src="/tex/a1dbd44a6b2fbeb770fcf69f85a4a3df.svg?invert_in_darkmode&sanitize=true" align=middle width=277.1488236pt height=49.315569599999996pt/></p>

where <img src="/tex/99751e94989c68f9be0f6aa442bc80d5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.302373649999986pt height=22.831056599999986pt/> is a penalty parameter and <img src="/tex/a44ff4154fc3bc708e9e752a14051324.svg?invert_in_darkmode&sanitize=true" align=middle width=49.762892849999986pt height=24.65753399999998pt/> is the minimal cost over the segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/>. The penalty <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> is understood as an additional cost when introducing a new segment. 

The cost for data <img src="/tex/7d54d6946c15cd803b31e7da5077c813.svg?invert_in_darkmode&sanitize=true" align=middle width=40.78280084999999pt height=14.15524440000002pt/> with linear interpolation from value <img src="/tex/44bc9d542a92714cac84e01cbbb7fd61.svg?invert_in_darkmode&sanitize=true" align=middle width=8.68915409999999pt height=14.15524440000002pt/> to value <img src="/tex/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?invert_in_darkmode&sanitize=true" align=middle width=7.054796099999991pt height=22.831056599999986pt/> is given by

<p align="center"><img src="/tex/77e662a6ebeab23830009ce77fbd9557.svg?invert_in_darkmode&sanitize=true" align=middle width=358.5926619pt height=48.39056475pt/></p>

which can be computed in constant time with the formula 

<p align="center"><img src="/tex/da574314e6b90740e74ddfb986bd58cb.svg?invert_in_darkmode&sanitize=true" align=middle width=484.25056514999994pt height=39.887022449999996pt/></p>

<p align="center"><img src="/tex/3625658a3a08ace419cc04fe7315f856.svg?invert_in_darkmode&sanitize=true" align=middle width=351.10582155pt height=34.3600389pt/></p>

where

<p align="center"><img src="/tex/1aab79a04b5886855193911c1cff1343.svg?invert_in_darkmode&sanitize=true" align=middle width=512.84137905pt height=47.02068525pt/></p>

To address the continuity constraint, we introduce the function <img src="/tex/ae890e198f4a771aa52230b69322c12f.svg?invert_in_darkmode&sanitize=true" align=middle width=74.25482954999998pt height=24.65753399999998pt/> which is the optimal penalized cost up to position <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> with a last infered value equal to <img src="/tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode&sanitize=true" align=middle width=8.55786029999999pt height=14.15524440000002pt/> (at position t). The idea is then to update a set

<p align="center"><img src="/tex/64b17ac0eaa9a23ba8f6e8beeef8b6b3.svg?invert_in_darkmode&sanitize=true" align=middle width=239.36380874999998pt height=16.438356pt/></p>

at any time step <img src="/tex/985dbad85d61b07e704840368824ee09.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/>. <img src="/tex/6005746712b75a74769ef4baa4d8fae3.svg?invert_in_darkmode&sanitize=true" align=middle width=32.40983954999999pt height=14.15524440000002pt/> and <img src="/tex/c1a3766b77f3cfa5e53205c84f3a7b1b.svg?invert_in_darkmode&sanitize=true" align=middle width=34.21767194999999pt height=14.15524440000002pt/> are the bounds of the interval of possible ending values for the considered data to segment. They can be determined in a preprocessing step.

The new update with continuity constraint takes the form

<p align="center"><img src="/tex/ab4ef0926386dff2055b015b92895c58.svg?invert_in_darkmode&sanitize=true" align=middle width=368.09542769999996pt height=31.6657044pt/></p>


where the presence of the same value <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/> in <img src="/tex/e05578c58a4f9e4f1301d4cf24e3234c.svg?invert_in_darkmode&sanitize=true" align=middle width=20.387619449999992pt height=22.465723500000017pt/> and the cost realizes the continuity constraint. At initial step we simply have <img src="/tex/cf8c7572e6ac52e95e69475a2b4281e3.svg?invert_in_darkmode&sanitize=true" align=middle width=86.58177495pt height=24.65753399999998pt/>. 


<a id="sf"></a>

## The slopeOP function

We install the package from Github:


```r
#devtools::install_github("vrunge/slopeOP")
library(slopeOP)
```

We simulate data with the function `slopeData` given the indices for extremal values and in second vector these values to reach for corresponding indices. The last parameter is the standard deviation of a normal standard noise.


```r
data <- slopeData(c(1,100,200,300,500), c(0,1,0,3,2), 1)
```

The changepoint detection is done by using the function `slopeOP`


```r
slopeOP(data, c(0,1,2,3), 10)
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

In `slopeOP` function, the parameter `type` is `channel` by default. With type equal to `channel` we use the monotonicity property in optimal cost matrix to reduce time complexity. If it is equal to `pruning` we prune some positions using a theorem taking into account unseen data.

<a id="options"></a>

## Options for constraining inference


Parameter `constraint` can be set to `up` which corresponds to a restriction to nondecreasing vector parameter.


```r
data <- slopeData(c(1,150,200,350,500,750,1000), c(71,73,70,75,77,73,80), 1)
slopeOP(data, 71:80, 5, constraint = "up")
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

With `constraint` equal to `updown` the infered vector is unimodal (with a maximum).


```r
data <- slopeData(c(1,150,200,350,500,750,1000), c(71,73,70,75,78,73,75), 1)
slopeOP(data, 71:80, 5, constraint = "updown")
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
data <- slopeData(c(1,30,40,70,100,150,200),c(70,80,70,80,70,80,70),0.5)
slopeOP(data,70:80,5, constraint = "smoothing", minAngle = 170)
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



[Back to Top](#top)

