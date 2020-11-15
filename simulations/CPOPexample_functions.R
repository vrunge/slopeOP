###
#Functions used for the examples in the CPOP paper:
# Maidstone, Fearnhead and Letchford
# "Detecting changes in slope with an $L_0$ penalty"
#
###

###
# Simulate Data
###

#simulate wave1
#K determines length -- overall length is K*1408
#sig is residual sd
##output is data and mean function and changepoints
rwave1=function(K=1,sig=1)
{
K=as.integer(K)
if(K<1) K=1
tau=c(0,256,512,768,1024,1152,1280,1344,1408)*K
beta=(-cumsum(c(0,1,-2,3,-4,5,-6,7)*2^-6)+2^-8)/K ##slopes for each segment
alpha=1
##calculate mu
m=length(tau)-1 ##number of segments
n=tau[m+1]
mu=rep(0,n);
for(i in 1:m){
    mu[(tau[i]+1):tau[i+1]]=alpha+beta[i]*(1:(tau[i+1]-tau[i]))
    alpha=mu[tau[i+1]]
    }
return(list(y=rnorm(length(mu),mu,sig),mu=mu,cps=tau[-c(1,m+1)]))
}

###simulate wave2
## N= segment length
## K=number of segments
## sig= residual sd
##output is data and mean function and changepoints

rwave2=function(N=150,K=10,sig=1){
if(K<2) K=2 ##minimum number of segments
tau=0:K*N

vec=c(0,(-1)^(2:K))
beta=(-cumsum(vec*(2^-5))+2^-6)*150/N ##slopes for each segment
alpha=-1/2
##calculate mu
m=length(tau)-1 ##number of segments
n=tau[m+1]
mu=rep(0,n);
for(i in 1:m){
    mu[(tau[i]+1):tau[i+1]]=alpha+beta[i]*(1:(tau[i+1]-tau[i]))
    alpha=mu[tau[i+1]]
    }

return(list(y=rnorm(length(mu),mu,sig),mu=mu,cps=tau[-c(1,m+1)]))

}

##simulate random
## nc is number of changespoints
## n is data length
## psi is sd of value of mean at each changepoint
## sig is residual sd

##output is data and mean function and changepoints

rRandom=function(n=1000,nc=9,psi=2,sig=1){
tau=(0:(nc+1))*as.integer(n/(nc+1))
m=length(tau)-1 ##number of segments
n=tau[m+1] #number of observations -- in case rounding
phi=rnorm(m+1,psi)
mu=rep(0,n)
for(i in 1:m){
    mu[(tau[i]+1):tau[i+1]]=phi[i]+(phi[i+1]-phi[i])*((tau[i]+1):tau[i+1] -tau[i])/(tau[i+1]-tau[i])
    }
  return(list(y=rnorm(length(mu),mu,sig),mu=mu,cps=tau[-c(1,m+1)]))

}

####################################################
# Functions to run NOT trend-filtering and CPOP
######################################################

library("not")
library("l1tf")
source("CPOP.R")
source("estimateCPOP.R")

##run CPOP
##input data, estimate of sigma (sigsquared=sig^2) and beta
CPOP.run=function(z,beta=NULL,sig=1,useC=T)
    {
        if(length(beta)!=1) beta=2*log(length(z)) ##BIC penalty
        ans<-CPOP(z,beta = beta,sigsquared = sig^2,useC=useC,useCprune=useC)
        CPS=ans[[2]]
	est=estimate.CPOP(z,CPS)
        
        return(list(cpt=CPS[-c(1,length(CPS))],fit=est$f))
    }



##runs the not funtion for pcwsLinContMean
##input data and truth
##output fit and changepoint locations
not.run=function(z,M=1e4,q.max=25){
##not output
w=not(z,M=M,method="not",contrast="pcwsLinContMean")
     fo <- features(w,q.max=q.max)
    return(list(cpt=fo$cpt,fit=predict(w,cpt=fo$c)))
}

##run trend-filtering
##input data; sigma; lambdas for grid search
##estimate based on SIC criteria
##output fit and changepoint locations
tf.run=function(z,sig=1,lambdas=exp(seq(log(10),log(1000),length=250)) )
    {
        M=length(lambdas)
        n=length(z)
        SIC=rep(0,M)
        for(i in 1:M){
            lambda<-lambdas[i]
            ans.tf<-l1tf(z,lambda) ##l1tf not there
            SIC[i]=sum((z-ans.tf)^2)
            dans<-round(diff(ans.tf),digits = 5 )
            CP<-c()
            for(ii in 1:(length(dans)-1)){
                if(dans[ii]!=dans[ii+1]){
                    CP<-c(CP,ii)
                }
            }
            SIC[i]=SIC[i]/sig^2+log(n)*length(CP)
        }
        k=which.min(SIC)
        ans.tf<-l1tf(z,lambdas[k]) ##l1tf not there
        #MSE=sum((mu-ans.tf)^2)
             dans<-round(diff(ans.tf),digits = 5 )
         CP<-c()
            for(ii in 1:(length(dans)-1)){
                if(dans[ii]!=dans[ii+1]){
                    CP<-c(CP,ii)
                }
            }
        
        return(list(cpt=CP,fit=ans.tf,lam=lambdas[k]))
        
    }

