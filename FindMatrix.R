# This script transforms a long time series into a nx(core) matrix size, where n is the number of rows and core the number of columns.
# core should be set as the number of cores used to estimate the models.
# The purpose is to speed up the estimation processes by using parallel estimation.
# w is in-sample (training) length series (usually 1000 or 2500 financial returns days)
# reto is the raw series.

FindMatrix<-function(reto,core, w){
  
  library(gdata)
  library(primes)
  
  
  ma<-length(reto)-w
  
  repeat{
    
    if( ma %% core != 0){
      ma = ma-1 }
    if(ma %% core == 0){
      break
    }
  }
  
  rg<-ma/core
  mb<- ma+w
  
  pret<-(reto[(length(reto)- mb +1):(length(reto))])
  
  rea<-pret[(w+1):length(pret)]
  
  pret<-(reto[(length(reto)- mb +1):(length(reto)-1)])
  
  X<-as.numeric(reto[(length(reto)- mb +1):length(reto)]) #as.numeric
  
  Fu <- matrix(0, (w+rg), core)
  
  Fu[,1] <- X[1:(w+rg)]
  
  for (i in 2:core) {
    
       Fu[,i] = X[((i-1)*rg+1):(w+(i*rg))]
    
       }
 
  a<-tail(X)
  b<-tail(Fu[,core])
  
  resp<-identical(a,b)
  
  Fu<-Fu[1:(nrow(Fu)-1),]
  q<-Fu[(w:nrow(Fu)),(2:ncol(Fu))]
  q<-unmatrix(q, byrow=F)
  Q<-c(Fu[,1],q)
  
  plot(Q-pret)
  
  return(list(Fu=Fu, pret=pret, identical=resp, rea=rea))
  
}
