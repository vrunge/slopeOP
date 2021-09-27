

##### ##### ##### ##### ##### ##### #####
##### standard deviation estimation #####
##### ##### ##### ##### ##### ##### #####

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
  if(method == "SD")
  {
    return(sd(diff(x)/sqrt(2)))
  }
  return(NULL)
}



##### ##### ##### ##### ##### ##### #####
##### estim_sigma change in slope PB ####
##### ##### ##### ##### ##### ##### #####

# the estim_sigma function gives an estimate of the standard deviation
# with three possible methods (HALL, MAD, sdHallDiff) based on the previous sdDiff function

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
    y <- signal + rnorm(n, 0, sigma)
    sdMad <- c(sdMad, sdDiff(y, method = "MAD"))
    sdHall <- c(sdHall, sdDiff(y, method = "HALL"))
    sdHallDiff <- c(sdHallDiff, sdDiff(diff(y), method = "HALL")/corrector)
  }

  response <- list(sdMad = sdMad, sdHall = sdHall, sdHallDiff = sdHallDiff)
  return(response)
}





########                ########
######## define signals ########
########                ########

sigma <- 1

n <- 33
signal1 <- c(seq(from = 0, to = n, length.out = n), seq(from = n, to = 0, length.out = n),  seq(from = 0, to = n, length.out = n))
signal2 <- -n*cos(3*pi/(3*n)*0:(3*n-1))  ## with the minus

########       ########
######## PLOTS ########
########       ########

par(mfrow = c(2, 2))
par(mar=c(2,5,1,1))

plot(signal1 + rnorm(n, 0, 1),  type ="p", pch = 3, cex = 1.5, cex.axis = 1.3, col = 1, xlab = "", ylab = "Signal 1", main = "", cex.lab = 1.8)
plot(signal1 + rnorm(n, 0, 5),  type ="p", pch = 3, cex = 1.5, cex.axis = 1.3, col = 1, xlab = "", ylab = "", main = "")
plot(-signal2 + rnorm(n, 0, 1),  type ="p", pch = 3, cex = 1.5, cex.axis = 1.3, col = 1, xlab = "", ylab = "Signal 2", main = "", cex.lab = 1.8)
plot(-signal2 + rnorm(n, 0, 5),  type ="p", pch = 3, cex = 1.5, cex.axis = 1.3, col = 1, xlab = "", ylab = "", main = "")


########             ########
######## SIMULATIONS ########
########             ########


p <- 10000 ###NB SIMULATIONS

#results
M1 <- NULL; M2 <- NULL; M3 <- NULL; M4 <- NULL;
H1 <- NULL; H2 <- NULL; H3 <- NULL; H4 <- NULL;
HD1 <- NULL; HD2 <- NULL; HD3 <- NULL; HD4 <- NULL;

for(sigma in 1:5)
{
  print(sigma)
  #estimation
  res1 <- estim_sigma(signal1, sigma, p)
  res2 <- estim_sigma(signal2, sigma, p)

  #results
  M1 <- c(M1, mean(res1$sdMad))
  M2 <- c(M2, sd(res1$sdMad))
  M3 <- c(M3, mean(res2$sdMad))
  M4 <- c(M4, sd(res2$sdMad))

  H1 <- c(H1, mean(res1$sdHall))
  H2 <- c(H2, sd(res1$sdHall))
  H3 <- c(H3, mean(res2$sdHall))
  H4 <- c(H4, sd(res2$sdHall))

  HD1 <- c(HD1, mean(res1$sdHallDiff))
  HD2 <- c(HD2, sd(res1$sdHallDiff))
  HD3 <- c(HD3, mean(res2$sdHallDiff))
  HD4 <- c(HD4, sd(res2$sdHallDiff))

  boxplot(res1$sdMad, res1$sdHall, res1$sdHallDiff)
  boxplot(res2$sdMad, res2$sdHall, res2$sdHallDiff)
}

M1; M2; M3; M4
H1; H2; H3; H4
HD1; HD2; HD3; HD4

paste(signif(M1,3),collapse = " & ");paste(signif(M2,2),collapse = " & ");paste(signif(M3,3),collapse = " & ");paste(signif(M4,2),collapse = " & ")
paste(signif(H1,3),collapse = " & ");paste(signif(H2,2),collapse = " & ");paste(signif(H3,3),collapse = " & ");paste(signif(H4,2),collapse = " & ")
paste(signif(HD1,3),collapse = " & ");paste(signif(HD2,2),collapse = " & ");paste(signif(HD3,3),collapse = " & ");paste(signif(HD4,2),collapse = " & ")



####################################################################################
##############
############## more changes (review)
##############


########                ########
######## define signals ########
########                ########

signal1 <- NULL
s <- 5*sample(0:6, size = 11, replace = TRUE)
for(i in 1:10) signal1 <-  c(signal1[-length(signal1)], seq(from = s[i], to = s[i+1], length = 11))
signal1 <-  signal1[-length(signal1)]

n <- 100
nbChange <- 10
signal2 <- -n/3*cos(nbChange*pi/(n)*0:(n-1))


########       ########
######## PLOTS ########
########       ########

par(mfrow = c(1, 2))
par(mar=c(2,5,1,1))

plot(signal1,  type ="b", cex = 1, cex.axis = 1, col = 1, xlab = "", ylab = "Signal 1", main = "", cex.lab = 2)
plot(signal2,  type ="b", cex = 1, cex.axis = 1, col = 1, xlab = "", ylab = "Signal 2", main = "", cex.lab = 2)


########             ########
######## SIMULATIONS ########
########             ########

p <- 10000 ###NB SIMULATIONS

#results
M1 <- NULL; M2 <- NULL; M3 <- NULL; M4 <- NULL;
H1 <- NULL; H2 <- NULL; H3 <- NULL; H4 <- NULL;
HD1 <- NULL; HD2 <- NULL; HD3 <- NULL; HD4 <- NULL;

for(sigma in 1:5)
{
  print(sigma)
  #estimation
  res1 <- estim_sigma(signal1, sigma, p)
  res2 <- estim_sigma(signal2, sigma, p)

  #results
  M1 <- c(M1, mean(res1$sdMad))
  M2 <- c(M2, sd(res1$sdMad))
  M3 <- c(M3, mean(res2$sdMad))
  M4 <- c(M4, sd(res2$sdMad))

  H1 <- c(H1, mean(res1$sdHall))
  H2 <- c(H2, sd(res1$sdHall))
  H3 <- c(H3, mean(res2$sdHall))
  H4 <- c(H4, sd(res2$sdHall))

  HD1 <- c(HD1, mean(res1$sdHallDiff))
  HD2 <- c(HD2, sd(res1$sdHallDiff))
  HD3 <- c(HD3, mean(res2$sdHallDiff))
  HD4 <- c(HD4, sd(res2$sdHallDiff))

  boxplot(res1$sdMad, res1$sdHall, res1$sdHallDiff)
  boxplot(res2$sdMad, res2$sdHall, res2$sdHallDiff)
}

M1; M2; M3; M4
H1; H2; H3; H4
HD1; HD2; HD3; HD4

paste(signif(M1,3),collapse = " & ");paste(signif(M2,2),collapse = " & ");paste(signif(M3,3),collapse = " & ");paste(signif(M4,2),collapse = " & ")
paste(signif(H1,3),collapse = " & ");paste(signif(H2,2),collapse = " & ");paste(signif(H3,3),collapse = " & ");paste(signif(H4,2),collapse = " & ")
paste(signif(HD1,3),collapse = " & ");paste(signif(HD2,2),collapse = " & ");paste(signif(HD3,3),collapse = " & ");paste(signif(HD4,2),collapse = " & ")




