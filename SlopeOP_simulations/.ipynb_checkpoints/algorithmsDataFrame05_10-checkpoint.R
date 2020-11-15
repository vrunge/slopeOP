args = commandArgs(trailingOnly=TRUE)
columns = NULL

# test if there is at least one argument: if not, return an error
if (length(args) < 4) {
  stop("Four arguments needed: algorithm, data file, bêta, working directory", call.=FALSE)
} else {
  # default output file
  	algo = args[1];
  	fileCSV = args[2];
  	beta = as.double(args[3]);
  	wd = args[4];

	  if (length(args) > 4) {
	  	column_from = as.numeric(args[5]);
	  	column_to = as.numeric(args[6]); 
	  	columns = column_from:column_to;
	  }
}

### ### ### IMPORT packages ### ### ### 
#devtools::install_github("vrunge/slopeOP", force = TRUE)
library(slopeOP)
#devtools::install_github("vrunge/gfpop", force = TRUE)
library(gfpop)
library(splines)

### ### ### ### ### ### ### ### ### ### ### ### 
setwd(dir = wd)
### ### ### ### ### ### ### ### ### ### ### ### 

### CPOP ### 
dyn.load("coeff.updateR.so")
dyn.load("prune2R.so")
library("not")
library("l1tf") #Trouver sur github
source("CPOPexample_functions.R")

#### build the signal from changepoints and states
buildSignal <- function(cp, states)
{
  steps <- diff(states)/diff(cp)
  res <- rep(steps, diff(cp))
  return(c(states[1], cumsum(res) + states[1]))
}

#### build the signal from changepoints by gfpop function
buildSignal_fpop <- function(cp, means, signal)
{
  res <- cumsum(rep(means, diff(c(0,cp + 1))))
  const <- mean(res - signal) #minimum of a quadratics
  return(res - const)
}



algoChpt <- function(fileCSV, algo = "slopeOP", states = 0:80, beta = 2, columns = NULL)
{
  df <- read.csv2(fileCSV, header = TRUE, sep = ",", dec=".")
  dfresponse <- df
  n <- dim(df)[1]
  
  if(is.null(columns)){columns <- 2:dim(df)[2]}
  for(i in columns)
  {
    strin <- strsplit(colnames(df), "_")[[i]][1] #first element = sigma
    sigma <- as.numeric(substr(strin,2,nchar(strin))) #delete X + numeric
    penalty <- beta * sigma * sigma * log(n) #get penalty value
    
    if(algo == "SlopeOP")
    {
      res <- slopeOP(df[,i], states, penalty)
      chgpts <- res$changepoints
      stateValues <- res$parameters
      dfresponse[,i] <- buildSignal(chgpts, stateValues)
    }
    
    if(algo == "OP2D")
    {
      res <- linearOP(0:(n-1), df[,i], penalty)
      chgpts <- res$x + 1
      stateValues <- res$y
      dfresponse[,i] <- buildSignal(chgpts, stateValues)
    }
    
    if(algo == "CPOP")
    {
      cpop_penalty <- beta * log(n)
      res <- CPOP.run(df[,i] / sigma, beta = cpop_penalty)
      dfresponse[,i] <- res$fit * sigma
    }
    
    if(algo == "FPOP")
    {
      newData <- df[2:n,i] - df[1:(n-1),i] #lag1 data
      res <- gfpop(newData, mygraph = graph(penalty = penalty, type = "std"))
      dfresponse[,i] <- buildSignal_fpop(res$changepoints, res$parameters, df[,1])
    }
    
    if(algo == "RFPOP")
    {
      newData <- df[2:n,i] - df[1:(n-1),i] #lag1 data
      res <- gfpop(newData, mygraph = graph(penalty = penalty, type = "std", K = 3*sigma))
      dfresponse[,i] <- buildSignal_fpop(res$changepoints, res$parameters, df[,1])
    }
    
  }
  write.csv(x = dfresponse, file = paste0("simu_",substr(fileCSV,1,nchar(fileCSV)-4),"_",algo,"_",sprintf("%.3f", beta),".csv"))
}


##################### Execution #########################

algoChpt(fileCSV=fileCSV, algo=algo, beta=beta, columns=columns)




