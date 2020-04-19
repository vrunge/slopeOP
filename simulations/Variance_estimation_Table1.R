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

####

estim_sigma <- function(signal, sigma, p)
{
  d0 <- 0.1942
  d1 <- 0.2809
  d2 <- 0.3832
  d3 <- -0.8582
  corrector <- sqrt(d3^2 + (d2-d3)^2 + (d1-d2)^2 + (d0-d1)^2 + d0^2)
  
  sdMad <- NULL
  sdHall <- NULL
  sdHallDiff <- NULL

  n <- length(signal)
  for(i in 1:p)
  {
    y <- signal + rnorm(n,0,sigma)
    sdMad <- c(sdMad,sdDiff(y,method = "MAD"))
    sdHall <- c(sdHall,sdDiff(y,method = "HALL"))
    sdHallDiff <- c(sdHallDiff,sdDiff(diff(y),method = "HALL")/corrector)
  }

  response <- list(sdMad = sdMad, sdHall = sdHall, sdHallDiff = sdHallDiff)
  return(response)
}



###############

sigma <- 1

n <- 333
signal1 <- c(seq(from = 0, to = n, length.out = n), seq(from = n, to = 0, length.out = n),  seq(from = 0, to = n, length.out = n))
signal2 <- n*cos(3*pi/(3*n)*0:(3*n-1))
p <- 100

#estimation
res1 <- estim_sigma(signal1, sigma, p)
res2 <- estim_sigma(signal2, sigma, p)

#results
mean(res1$sdMad)
mean(res1$sdHall)
mean(res1$sdHallDiff)
sd(res1$sdMad)
sd(res1$sdHall)
sd(res1$sdHallDiff)

mean(res2$sdMad)
mean(res2$sdHall)
mean(res2$sdHallDiff)
sd(res2$sdMad)
sd(res2$sdHall)
sd(res2$sdHallDiff)

boxplot(res1$sdMad, res1$sdHall, res1$sdHallDiff)
boxplot(res2$sdMad, res2$sdHall, res2$sdHallDiff)
