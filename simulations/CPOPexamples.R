###
# Example use of code
#
###

###example of analysis for simulation study

source("CPOPexample_functions.R") ##load functions

x=rwave1() ##simulate data under basic wave1 scenario
z=x[[1]] ##data
mu=x[[2]] ##mean
taus=x[[3]] ##changepoints
ns=max(diff(taus)) ##largest segment
n=length(z) ##data length

sig.hat=mad(diff(diff(z)))/sqrt(6) ##crude estimate of residual sd

##NOT OUTPUT
w=not.run(z,M=1e5) ##run NOT

cat("NOT MSE=",sum((mu-w$f)^2)/n,"\n")
cat("NOT number changepoints=",length(w$c),"\n")

DIFF=abs(matrix(w$c,nrow=length(taus),ncol=length(w$c),byr=T)-matrix(taus,nrow=length(taus),ncol=length(w$c),byr=F))
cat("NOT d_H =",max(apply(DIFF,1,min),apply(DIFF,2,min))/ns,"\n")

##Trend-filtering
 w=tf.run(z/sig.hat) ##standardised data to have residual sd approx 1  
w$f=w$f*sig.hat ##correct estimate of mean for standardisation

cat("TF MSE=",sum((mu-w$f)^2)/n,"\n")
cat("TF number changepoints=",length(w$c),"\n")

DIFF=abs(matrix(w$c,nrow=length(taus),ncol=length(w$c),byr=T)-matrix(taus,nrow=length(taus),ncol=length(w$c),byr=F))
cat("TF d_H =",max(apply(DIFF,1,min),apply(DIFF,2,min))/ns,"\n")

##CPOP
w=CPOP.run(z/sig.hat)##standardised data to have residual sd approx 1  
w$f=w$f*sig.hat ##correct estimate of mean for standardisation

cat("CPOP MSE=",sum((mu-w$f)^2)/n,"\n")
cat("CPOP number changepoints=",length(w$c),"\n")

DIFF=abs(matrix(w$c,nrow=length(taus),ncol=length(w$c),byr=T)-matrix(taus,nrow=length(taus),ncol=length(w$c),byr=F))
cat("CPOP d_H =",max(apply(DIFF,1,min),apply(DIFF,2,min))/ns,"\n")

###############Example use of CROPS with CPOP
source("CROPS_CPOP.R")

n=length(z)
out=CROPS.CPOP(z,sig.hat^2,min_pen=log(n),max_pen=2*log(n),PRINT=T,useC=T,useCprune=T) 

out[[1]] ##summary
out[[2]][[1]] ##one optimal segmentation


