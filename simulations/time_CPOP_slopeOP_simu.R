

##############################
###### CPOP and slopeOP ######
##############################

install.packages("not")
library(devtools)
install_github("hadley/l1tf")
devtools::install_github("vrunge/slopeOP", force = TRUE)

dyn.load("coeff.updateR.so")
dyn.load("prune2R.so")

library("not")
library("l1tf")
source("CPOPexample_functions.R")
library(slopeOP)
library(parallel)


################################
###### one.simu functions ######
################################

one.simuCPOP <- function(i, n, sigma, signal = 1)
{
  penalty <- 2 * sigma^2 * log(n)

  ### CHOOSE SIGNAL SHAPE
  if(signal == 1)
  {y <- c(seq(from = 10, to = 50, length.out = n/2), seq(from = 50, to = 10, length.out = n/2)) + rnorm(n, 0, sigma)}
  if(signal == 2)
  {
    signal <- NULL
    s <- c(rep(c(10,50), 10), 10)
    for(i in 1:19) signal <-  c(signal[-length(signal)], seq(from = s[i], to = s[i+1], length = (n %/% 19) + 1))
    signal <-  c(signal[-length(signal)], seq(from = s[20], to = s[20] + (n %% 19)*( s[21]- s[20])/((n %/% 19)), length = n %% 19))
    y <- signal + rnorm(n, 0, sigma)
  }

  start_time1 <- Sys.time()
  result <- CPOP.run(y / sqrt(penalty))
  end_time1 <- Sys.time()
  return(unclass(end_time1 - start_time1)[1])
}

one.simuSLOPE_OP <- function(i, n, sigma, signal = 1)
{
  penalty <- 2 * sigma^2 * log(n)

  ### CHOOSE SIGNAL SHAPE
  if(signal == 1)
    {y <- c(seq(from = 10, to = 50, length.out = n/2), seq(from = 50, to = 10, length.out = n/2)) + rnorm(n, 0, sigma)}
  if(signal == 2)
  {
    signal <- NULL
    s <- c(rep(c(10,50), 10), 10)
    for(i in 1:19) signal <-  c(signal[-length(signal)], seq(from = s[i], to = s[i+1], length = (n %/% 19) + 1))
    signal <-  c(signal[-length(signal)], seq(from = s[20], to = s[20] + (n %% 19)*( s[21]- s[20])/((n %/% 19)), length = n %% 19))
    y <- signal + rnorm(n, 0, sigma)
  }

  start_time2 <- Sys.time()
  result <- slopeOP(y, 0:60, penalty, type = "channel")
  end_time2 <- Sys.time()
  return(unclass(end_time2 - start_time2)[1])
}




################################
###### data length vector ######
################################

#mytest <- seq(from = 100, to = 1500, by = 10)
my_n_vector_LOG <- seq(from = log(100), to = log(1500), by = log(15)/50)
my_n_vector <- round(exp(my_n_vector_LOG))
my_n_vector <- 2*round(my_n_vector/2)
my_n_vector



##################################
###### Initialize dataframe ######
##################################
p <- 100
df <- data.frame(matrix( nrow = 4 * length(mytest), ncol = 3 + p))
colnames(df) <- c("type", "n", "sigma", 1:p)
dim(df)


####################################
###### simulations multicores ######
####################################

signal_type <- 2 #swith to signal_type = 1
nbCores <- 50
j <- 1

for(n in my_n_vector)
{
  print(c("n =", n))
  liste1 <- mclapply(1:p, FUN = one.simuCPOP,
                     n = n,
                     sigma = 3,
                     signal = signal_type,
                     mc.cores = nbCores)
  liste2 <- mclapply(1:p, FUN = one.simuSLOPE_OP,
                     n = n,
                     sigma = 3,
                     signal = signal_type,
                     mc.cores = nbCores)

  df[j ,] <- c("CPOP", n, 3, do.call(cbind, liste1))
  df[j+1, ] <- c("slopeOP", n, 3, do.call(cbind, liste2))
  j <- j + 2

  print(c("n =", n))
  liste3 <- mclapply(1:p, FUN = one.simuCPOP,
                     n = n,
                     sigma = 24,
                     signal = signal_type,
                     mc.cores = nbCores)
  liste4 <- mclapply(1:p, FUN = one.simuSLOPE_OP,
                     n = n,
                     sigma = 24,
                     signal = signal_type,
                     mc.cores = nbCores)

  df[j ,] <- c("CPOP", n, 24, do.call(cbind, liste3))
  df[j+1, ] <- c("slopeOP", n, 24, do.call(cbind, liste4))
  j <- j + 2
}


##########################################
###### saving dataframe in csv file ######
##########################################

write.csv(df, "timeSIGNAL.csv", row.names = FALSE)

