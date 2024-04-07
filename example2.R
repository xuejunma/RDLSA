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

source("RDLSA_SCAD.R")

N <- 1e5
ps <- 20
case_error <- 1
case_x <- 2
hete_level <- 1 # Possible values are 1 (the default), 5 or 10.
M <- 500
K_work <- c(10,20,500,100,200)
c_work <- length(K_work)
sel_rlt_name <- c('CDMR_SCAD','RDLSA_SCAD')
sel_rlt <- matrix(0,length(sel_rlt_name)*c_work,ps+6)
sel_rlt[,1] <- rep(K_work,each=length(sel_rlt_name))
colnames(sel_rlt) <- c('K',c(paste('x',c(0:ps))),'ASE','FP','FN','CM')
rownames(sel_rlt) <- rep(sel_rlt_name,c_work)

index_sel <- 1
sele <- list()
for (k in 1:c_work) {
  K <- K_work[k]
  fit.cv <- chunk.f(n=N,K=K,partition = "uniform",type = "consecutive")
  
  ## parallel computation of M simulations
  c2 <- makeCluster(4)
  registerDoParallel(c2)
  b<-foreach(m = 1:M,  .combine=rbind, .packages = c("limSolve","Rfast","quantreg","VGAM")) %dopar% {
    var_sele(n=N,p=ps,K=K,case_x=case_x,hete_level=hete_level,case_error= case_error,fit.cv=fit.cv)
  }
  stopCluster(c2)
  
  CDMR_sele <- b[,1:(ps+5)]
  RDLSA_sele <- b[,(ps+6):(2*ps+10)]
  
  sel_rlt[index_sel,-1] <- apply(CDMR_sele, 2, mean)
  sel_rlt[c(index_sel+1),-1] <- apply(RDLSA_sele, 2, mean)
  index_sel <- index_sel + 2
  sele[[k]] <- b
  
}
sel_rlt <-  sel_rlt[,-(ps+3)]
sel_rlt