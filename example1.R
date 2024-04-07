rm(list=ls())
library(limSolve)
library(Rfast)
library(quantreg)
library(VGAM)
library(MASS)
library(parallel)
library(iterators)
library(foreach)
library(doParallel)

source("RDLSA.R")

N <- 1e5
pv <- 8
case_error <- 1
case_x <- 2
hete_level <- 1 # Possible values are 1 (the default), 5 or 10.
M <- 500
K_work <- c(10,20,500,100,200)
c_work <- length(K_work)
result_name <- c('CDMR1','CDMR2','CDMR3','RDLSA1','RDLSA2','RDLSA3')
result <- matrix(0,length(result_name)*c_work,pv+4)
result[,1] <- rep(K_work,each=length(result_name))
colnames(result) <- c('K',c(paste('x',c(0:pv))),'AEE','Cost')
rownames(result) <- rep(result_name,c_work)

index_resu <- 1
esti <- list()
for (k in 1:c_work) {
  K <- K_work[k]
  fit.cv <- chunk.f(n=N,K=K,partition = "uniform",type = "consecutive")
  
  ## parallel computation of M simulations
  c1 <- makeCluster(4)
  registerDoParallel(c1)
  a<-foreach(m = 1:M,  .combine=rbind, .packages = c("limSolve","Rfast","quantreg","VGAM")) %dopar% {
    para.esti(n=N,p=pv,case_x=case_x,hete_level=hete_level,case_error=case_error,fit.cv=fit.cv,K=K)
  }
  stopCluster(c1)
  
  ASE_RDLSA1 <- a[,1:(pv+3)]
  ASE_RDLSA2 <- a[,(pv+4):(2*pv+6)]
  ASE_CDMR1 <- a[,(2*pv+7):(3*pv+9)]
  ASE_CDMR2 <- a[,(3*pv+10):(4*pv+12)]
  ASE_CDMR3 <- a[,(4*pv+13):(5*pv+15)]
  ASE_RDLSA3 <- a[,(5*pv+16):(6*pv+18)]
  
  result[index_resu,-1] <- apply(ASE_CDMR1,2,mean)
  result[c(index_resu+1),-1] <- apply(ASE_CDMR2,2,mean)
  result[c(index_resu+2),-1] <- apply(ASE_CDMR3,2,mean)
  result[c(index_resu+3),-1] <- apply(ASE_RDLSA1,2,mean)
  result[c(index_resu+4),-1] <- apply(ASE_RDLSA2,2,mean)
  result[c(index_resu+5),-1] <- apply(ASE_RDLSA3,2,mean)
  index_resu <- index_resu + 6
  
  esti[[k]] <- a
  
}
result <-  result[,-(pv+4)]
result
