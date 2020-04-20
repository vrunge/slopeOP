
## CPOP and slopeOP
library("not")
library("l1tf")
source("CPOP/CPOPexample_functions.R")
dyn.load("CPOP/coeff.updateR.so")
dyn.load("CPOP/prune2R.so")

#devtools::install_github("vrunge/slopeOP", force = TRUE)
library(slopeOP)

################################################

T1 <- NULL
T2 <- NULL
sigma <- 3  #### THE LEVEL OF NOISE
nbSimus <- 10

################################################

mytest <- seq(from = 100, to = 1500, by = 10)
for(n in mytest)
{
  #print(n)
  penalty <- 2 * sigma^2 * log(n)
  res1 <- NULL
  res2 <- NULL
  for(i in 1:nbSimus)
  {
    y <- c(seq(from = 10, to = 50, length.out = n/2), seq(from = 50, to = 10, length.out = n/2)) + rnorm(n,0,sigma)

    #TIME CPOP
    start_time1 <- Sys.time()
    result <- CPOP.run(y / sqrt(penalty))
    end_time1 <- Sys.time()
    res1 <- c(res1, unclass(end_time1 - start_time1)[1])

    #TIME slopeOP
    start_time2 <- Sys.time()
    s <- slopeOP(data = y, states = 0:60, penalty = penalty, type = "channel")
    end_time2 <- Sys.time()
    res2 <- c(res2, unclass(end_time2 - start_time2)[1])
  }
  T1 <- c(T1, mean(res1))
  T2 <- c(T2, mean(res2))
}

df <- data.frame(mytest,T1,T2)
write.csv(x = df, file = "time3.csv")
