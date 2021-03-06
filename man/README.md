---
output:
  html_document: default
  pdf_document: default
  word_document: default
---
<a id="top"></a>

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

The package `slopeOP` is designed to segment univariate data $y_{1:n} = \{y_1,...,y_n\}$ by a continuous piecewise linear signal with restriction on starting/ending values for the segments. The changepoint vector $\overline{\tau} = (\tau_0 < \cdots < \tau_{k+1}) \in \mathbb{N}^{k+2}$ defines the $k+1$ segments $\{\tau_i+1,...,\tau_{i+1}\}$, $i = 0,...,k$ with fixed bounds $\tau_0 = 0$ and  $\tau_{k+1} = n$. We use the set $S_n = \{\hbox{changepoint vector } \overline{\tau} \in \mathbb{N}^{k+2}\}$ to define the minimal global cost given by

$$Q_n = \min_{\overline{\tau} \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) + \beta \rbrace \right]\,,$$

where $\beta > 0$ is a penalty parameter and $\mathcal{C}(y_{u:v})$ is the minimal cost over the segment $\{u,...,v\}$. The penalty $\beta$ is understood as an additional cost when introducing a new segment. 

The cost for data $y_{\tau+1:t}$ with linear interpolation from value $a$ to value $b$ is given by

$$
\mathcal{C}(y_{\tau+1:t},a,b) = \sum_{i=\tau+1}^{t}\left(y_i - [a + (b-a)\frac{i-\tau}{t-\tau}]\right)^2\,,
$$

which can be computed in constant time with the formula 

$$\mathcal{C}(y_{\tau+1:t},a,b) = S_{t}^2 - S_{\tau}^2 +\frac{b^2-a^2}{2} + \frac{a^2 + ab +b^2}{3}(t-\tau) + \frac{(b-a)^2}{6(t-\tau)}$$

$$\quad\quad\quad\quad\quad\quad - \frac{2}{t-\tau} \left( (ta-b\tau)[S_{t}^1 - S_{\tau}^1] + (b-a)[S_{t}^+ - S_{\tau}^+]\right)\,,$$

where

$$
S^1_t = \sum_{i=1}^t y_i\quad , \quad S^2_t = \sum_{i=1}^t y_i^2\quad \hbox{and} \quad S^+_t = \sum_{i=1}^t iy_i\quad \hbox{for all} \,\, t \in \{1,...,n\}\,.
$$

To address the continuity constraint, we introduce the function $v \mapsto Q_t(v)$ which is the optimal penalized cost up to position $t$ with a last infered value equal to $v$ (at position t). The idea is then to update a set

$$
\mathcal{Q}_t = \{Q_t(v), v= v_{min},...,v_{max}\}\,,
$$

at any time step $t \in \{1,...,n\}$. $v_{min}$ and $v_{max}$ are the bounds of the interval of possible ending values for the considered data to segment. They can be determined in a preprocessing step.

The new update with continuity constraint takes the form

$$
Q_t(v) = \min_{0 \le \tau < t}\left( \min_{u}\{Q_{\tau}(u) + \mathcal{C}(y_{\tau+1:t},u,v) + \beta\}\right)\,,
$$


where the presence of the same value $u$ in $Q_{\tau}$ and the cost realizes the continuity constraint. At initial step we simply have $Q_0(v) = -\beta$. 


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
## $parameters
## [1] 0 1 0 3 2
## 
## $globalCost
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
## $changepoints
## [1]    1   63  254  356  828 1000
## 
## $parameters
## [1] 71 72 72 75 75 80
## 
## $globalCost
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
## $changepoints
## [1]    1  316  317  502  697 1000
## 
## $parameters
## [1] 71 73 74 78 74 74
## 
## $globalCost
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
## $changepoints
##  [1]   1   7  19  25  28  34  40  46  47  53  62  68  73  79  94 100 103
## [18] 149 154 160 200
## 
## $parameters
##  [1] 71 72 76 77 77 76 74 73 73 74 77 78 78 77 72 71 71 79 79 78 70
## 
## $globalCost
## [1] 319.9376
## 
## attr(,"class")
## [1] "slopeOP"
```



[Back to Top](#top)
