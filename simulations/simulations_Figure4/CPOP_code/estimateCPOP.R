###
# function to estimate the underlying mean function when fitting change-in-slope model
#
# input: y -- data
#        tau -- vector of changepoints [can include 0 as first and n as last]
#
# output: phi -- value of mean at (0,changepoints,n)
#         f -- vector of mean function at each time-point
###

estimate.CPOP<-function(y,tau){
  if(length(y)==1){return(list(100000,tau,c(y,y)))}
  else{
  n<-length(y)
###add 0 and n to tau if needed##
  if(tau[1]!=0) tau=c(0,tau)
  if(max(tau)<n) tau=c(tau,n)
#####
  k<-length(tau)-1
  S1<-0
    for(i in 1:n){
      S1[i+1]<-S1[i]+y[i]
    }
  segS1<-c()
  for(i in 1:k){
    segS1[i]<-S1[tau[i+1]+1]-S1[tau[i]+1]
  }
  segS2<-c()
  seglength<-c()
  for(i in 1:k){
    seglength[i]<-tau[i+1]-tau[i]
    segS2[i]<-sum(y[(tau[i]+1):tau[i+1]]*(1:seglength[i]))
  }
  ####calculation of b####
  b<-c()
  b[1]<-(-1)*segS1[1]+(1/seglength[1])*segS2[1]
  if(k>1){
  for(i in 2:k){
    b[i]<-(-1)*segS1[i]+(1/seglength[i])*segS2[i]-(1/seglength[i-1])*segS2[i-1]
  }}
  b[k+1]<-(-1/seglength[k])*segS2[k]
  ####calculation of A####
  A<-matrix(0,nrow=k+1,ncol=k+1)
  A[1,1]<-(seglength[1]+1)*(2*seglength[1]+1)/(3*seglength[1])-2
  A[1,2]<-seglength[1]+1-(seglength[1]+1)*(2*seglength[1]+1)/(3*seglength[1])
  if(k>1){
  for(i in 2:k){
    A[i,i-1]<-seglength[i-1]+1-(seglength[i-1]+1)*(2*seglength[i-1]+1)/(3*seglength[i-1])
    A[i,i]<-(seglength[i]+1)*(2*seglength[i]+1)/(3*seglength[i])+(seglength[i-1]+1)*(2*seglength[i-1]+1)/(3*seglength[i-1])-2
    A[i,i+1]<-seglength[i]+1-(seglength[i]+1)*(2*seglength[i]+1)/(3*seglength[i])
  }}
  A[k+1,k]<-seglength[k]+1-(seglength[k]+1)*(2*seglength[k]+1)/(3*seglength[k])
  A[k+1,k+1]<-(seglength[k]+1)*(2*seglength[k]+1)/(3*seglength[k])
  A<-(-1/(2))*A
  ####Thomas Algorithm to solve A*theta=b####
  if(A[1,]==rep(0,k+1) && A[,1]==rep(0,k+1)){
    Ad<-A[-1,]
    Ad<-Ad[,-1]
    bd<-b[-1]
    theta1<--b[2]
  } else{
    Ad<-A
    bd<-b
  }
  theta<-thomas(Ad,bd)
  if(A[1,]==rep(0,k+1) && A[,1]==rep(0,k+1)){
    theta<-c(theta1,theta)
  }
 
  ###now estimate the mean function
  x0=tau[-length(tau)]
  y0=theta[-length(theta)]
  x1=tau[-1]
  y1=theta[-1]
  est<-c()
  for(i in 1:length(x0)){
        est[(x0[i]+1):x1[i]]<-(y1[i]-y0[i])/(x1[i]-x0[i])*(1:(x1[i]-x0[i]))+y0[i]
  }
      

  return(list(phi=theta,f=est))
}
  }

##thomas method uses in estimate.CPOP
thomas<-function(A,r){
  n<-length(r)
  gam<-c()
  rho<-c()
  x<-c()
  if(n==1){return(r/A)}
  else{
  gam[1]<-A[1,2]/A[1,1]
  rho[1]<-r[1]/A[1,1]
  if(n>2){
  for(i in 2:(n-1)){
    gam[i]<-A[i,i+1]/(A[i,i]-A[i,i-1]*gam[i-1])
    rho[i]<-(r[i]-A[i,i-1]*rho[i-1])/(A[i,i]-A[i,i-1]*gam[i-1])
  }}
  rho[n]<-(r[n]-A[n,n-1]*rho[n-1])/(A[n,n]-A[n,n-1]*gam[n-1])
  x[n]<-rho[n]
  for(i in (n-1):1){
    x[i]<-rho[i]-gam[i]*x[i+1]
  }
  return(x)
}}


