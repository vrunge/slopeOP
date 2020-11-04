
# the sdDiff function gives an estimate of the standard deviation with two possible methods (HALL, MAD)

sdDiff <- function(x, method = 'HALL')
{
  if(is.numeric(x) == FALSE || length(x) < 2){stop('x is not a numeric vector of length > 1')}
  if(method == "HALL")
  {
    n = length(x)
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(x)
    mat[2, -n] = mat[2, -1]
    mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]
    return(sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3)))
  }
  if(method == "MAD")
  {
    return(mad(diff(x)/sqrt(2)))
  }
  return(NULL)
}

# the estim_sigma function gives an estimate of the standard deviation 
# with three possible methods (HALL, MAD, sdHallDiff) based on the previous sdDiff function

estim_sigma <- function(signal, sigma, nbSimus)
{
  # the new normalization coefficient
  d0 <- 0.1942
  d1 <- 0.2809
  d2 <- 0.3832
  d3 <- -0.8582
  corrector <- sqrt(d3^2 + (d2-d3)^2 + (d1-d2)^2 + (d0-d1)^2 + d0^2)
  
  #data generation
  n <- length(signal)
  y <- signal + rnorm(n, 0, sigma)
  
  sdMad <- sdDiff(y, method = "MAD")
  sdHall <- sdDiff(y, method = "HALL")
  sdHallDiff <- sdDiff(diff(y), method = "HALL")/corrector
  
  response <- c(sdMad = sdMad, sdHall = sdHall, sdHallDiff = sdHallDiff)
  return(response)
}
 
########             ########
######## SIMULATIONS ########
########             ########

nbSimus <- 100 #number of simulations

sigma <- 1 # We choose a noise level (standard deviation)
n <- 333 # data length is 3*n

signal1 <- c(seq(from = 0, to = n, length.out = n), seq(from = n, to = 0, length.out = n),  seq(from = 0, to = n, length.out = n))
signal2 <- n * cos(3*pi/(3*n) * 0:(3*n-1))

##### estimation nbSimus times of the three sd estimates #####
res_signal1 <- replicate(nbSimus, estim_sigma(signal1, sigma))
res_signal2 <- replicate(nbSimus, estim_sigma(signal2, sigma))

df_s1 <- as.data.frame(t(res_signal1))
df_s2 <- as.data.frame(t(res_signal2))

########         ########
######## results ########
########         ########

colMeans(df_s1)
colMeans(df_s2)

apply(df_s1, 2, sd)
apply(df_s2, 2, sd)


########          ########
######## boxplots ######## 
########          ########

## Simple boxplots
boxplot(df_s1)
boxplot(df_s2)

## Fancy boxplots
library(ggplot2)
library(reshape2)

df <- cbind(df_s1, df_s2)
colnames(df) <- c("sdMad_s1", "sdHall_s1", "sdHallDiff_s1", "sdMad_s2", "sdHall_s2", "sdHallDiff_s2")

ggplot(data = melt(df), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable)) + theme(legend.position="none")

