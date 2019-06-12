<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{An Introduction to gfpop}
-->


<a id="top"></a>

# slopeOP Vignette
### Vincent Runge
#### LaMME, Evry University
### June 11, 2019

> [Introduction](#intro)

> [The R functions](#rf)

<a id="intro"></a>

## Introduction

The package `slopeOP` is designed to segment univariate data $y_{1:n} = \{y_1,...,y_n\}$ by a continuous piecewise linear signal. The changepoint vector $\overline{\tau} = (\tau_0 < \cdots < \tau_{k+1}) \in \mathbb{N}^{k+2}$ defines the $k+1$ segments $\{\tau_i+1,...,\tau_{i+1}\}$, $i = 0,...,k$ with fixed bounds $\tau_0 = 0$ and  $\tau_{k+1} = n$. We use the set $S_n = \{\hbox{changepoint vector } \overline{\tau} \in \mathbb{N}^{k+2}\}$ to define the nonconstrained minimal global cost given by

$$Q_n = \min_{\overline{\tau} \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) + \beta \rbrace \right]\,,$$

where $\beta > 0$ is a penalty parameter and $\mathcal{C}(y_{u:v})$ is the minimal cost over the segment $\{u,...,v\}$. The penalty $\beta$ is understood as an additional cost when introducing a new segment.

The cost for data $y_{\tau+1:t}$ with linear interpolation from value $a$ to value $b$ is given by

$$
\mathcal{C}(y_{\tau+1:t},a,b) = \sum_{i=\tau+1}^{t}\left(y_i - [a + (b-a)\frac{i-\tau}{t-\tau}]\right)^2\,.
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


<a id="rf"></a>

## The R functions

We install the package from Github:


```r
#devtools::install_github("vrunge/slopeOP")
library(slopeOP)
```

We simulate data with the function `slopeData` given the indices for extremal values and in second vector these values to reach corresponding indices. The last parameter is the standard deviation of a normal standard noise.


```r
data <- slopeData(c(1,100,200,300,500), c(0,1,0,3,2), 1)
```

The changepoint detection is done by using the function `slopeOP`


```r
slopeOP(data, c(0,1,2,3), 10)
```

```
## changepoints
## [1]   1  78 197 293 500
##
## parameters
## [1] 0 1 0 3 2
##
## globalCost
## [1] 573.4121
##
## attr(,"class")
## [1] "slopeOP"
```

In `slopeOP` function, the parameter `type` is `null` by default. If it is equal to `channel` we use the monotonicity property in optimal cost matrix to reduce time complexity. If it is equal to `pruning` we prune some positions using a theorem taking into account unseen data.


[Back to Top](#top)