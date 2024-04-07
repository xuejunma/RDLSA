

### Description:
###   Calculating the first/second order gradients.
###Arguments:
###   x: the covariate matrix in K-th machines/blocks; 
###   y: the response variable in K-th machines/blocks;
###   beta0: the estimators of model parameters;
###   h: the bandwidth.
grad12 <- function(xk,yk,beta0,h=0.5){
  nk <- length(yk) 
  resk <- yk - xk%*%beta0
  wk <- resk*dnorm(resk/h) / h^3
  
  ## the first order gradients
  grad1Qk <- t(xk) %*% wk /nk
  
  ## the second order gradients.
  grad2Qk <- t(xk)%*%(c(dnorm(resk/h)* resk^2 / h^5 - dnorm(resk/h)/h^3)*xk) / nk
  
  return(list(grad1Qk=grad1Qk,grad2Qk=grad2Qk))
}


### Description:
###   Used to divide data into K blocks.
### Arguments:
###   n: an integer giving the number of observations/samples to be split into groups;
###   K: split data into K block;
###   R: an integer giving the number of replications for repeated K-fold cross-validation.
###      This is ignored for for leave-one-out cross-validation and other non-random splits 
###      of the data;
###   partition: A string specifying the type of the number of observations/samples contained 
###              in each block. Possible values are "stochastic" (default) or "uniform";
###   type: a character string specifying the type of folds to be generated; 
###         Possible values are "random" (the default), "consecutive" or "interleaved".
chunk.f <- function(n=100,K=5,R = 1,partition=c("stochastic","uniform"),
                    type = c("random", "consecutive", "interleaved")){
  n <- round(rep(n, length.out=1))
  if(!isTRUE(n > 0)) stop("'n' must be positive")
  K <- round(rep(K, length.out=1))
  if(!isTRUE((K > 1) && K <= n)) stop("'K' outside allowable range")
  type <- if(K == n) "leave-one-out" else match.arg(type)
  if(type == "random") {
    # random K-fold splits with R replications
    R <- round(rep(R, length.out=1))
    if(!isTRUE(R > 0)) R <- 1
    subsets <- replicate(R, sample(n))
  }else {
    # leave-one-out CV or non-random splits, replication not meaningful
    R <- 1
    subsets <- as.matrix(seq_len(n))
  }
  partition <-match.arg(partition)
  if(partition == "stochastic"){
    #weight <- runif(K,0,1)
    weight <- rnorm(K,mean = 1, sd=0.1)
    repeat{
      if((weight[which.max(weight)]-weight[which.min(weight)])<0.4)break
      weight <- rnorm(K,mean = 1, sd=0.1)
    }
    weight <- weight/sum(weight)
    nk <- n*weight
    nk <- round(nk)
    if(sum(nk)>n){
      nk[which.max(nk)] <- nk[which.max(nk)] - (sum(nk)-n)
    }else{
      nk[which.min(nk)] <- nk[which.min(nk)] + (n-sum(nk))
    }
    
    which <- rep.int(seq_len(K),nk)
    
  }else{
    which <- rep(seq_len(K), length.out=n)
    which <- rep.int(seq_len(K), tabulate(which))
  }
  folds <- list(n=n, K=K, R=R, subsets=subsets, which=which)
  class(folds) <- "cvFolds"
  folds
}


### Description:
###   The optimal bandwidth h is selected by the grid search method.
### Arguments:
###   x: the covariate matrix; y: the response variable;
###   fit.cv: index of data for different machines;
###   K: the number of machines/blocks;
###   Intercept: whether the model contains an intercept term.Defaults to FALSE.
select_h <- function(x, y, fit.cv, K=20, Intercept=FALSE){
  if(Intercept==TRUE) x <- cbind(1, x)
  n <- length(y)
  p <- dim(x)[2]
  
  ####h begin
  ###################selection h   begin 
  rh <- rep(0, 60)
  res <- rep(0, 60)
  for (j in 1:60) {
    rh_vec <- NULL
    res_vec <- NULL
    ###############loop begin
    for(kk in 1:K){
      index.kk <- fit.cv$subsets[fit.cv$which==kk]
      resid0 <- quantreg::rq(y[index.kk]~x[index.kk, ]-1,tau=0.5)$residuals
      bdh <- 0.5 * sd(resid0) * 1.02 ^ (j-1)
      gh <- mean( (-dnorm(resid0 / bdh) * resid0 / bdh^3)^2   )
      fh <- mean( dnorm(resid0 / bdh) * resid0^2 / bdh^5 - dnorm(resid0 / bdh) /bdh^3 )
      rh_vec[kk] <- gh / fh^2 /(sd(resid0)) ^ 2
      res_vec[kk] <- sd(resid0)
    }
    ###############loop end
    rh[j] <- mean(rh_vec)
    res[j] <- mean(res_vec) 
  }
  if(order(rh)[1]==1){
    h <-  0.5 * res[(order(rh)[1])] * 1.02 ^ (order(rh)[1] )
  }else{
    h <-  0.5 * res[(order(rh)[1] - 1)] * 1.02 ^ (order(rh)[1] - 1)
  }
  return(list(rh_full=rh, res_full=res_vec, h.opt=h))
}


### Description:
###   MEM Algorithm for Modal Regression
### Arguments:
###   theta: initial value of MEM Algorithm;
###   x: the covariate matrix; y: the response variable; h: the bandwidth;
###   tolset: criteria for stopping algorithm convergence;
###   Intercept: whether the model contains an intercept term.Defaults to FALSE.
modreg.MEM <- function(theta, y, x, h, tolset=10^(-6), Intercept=FALSE) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(Intercept==TRUE) x <- cbind(1, x)
  tol <- 1
  iter <- 0
  #loop begin
  while (tol > tolset) {
    #if tol <=0.000001 stop
    #E-step
    w <- dnorm((y - x %*% theta) / h)
    ww <- c(w / sum(w))
    #M-step
    thetanew <-  solve(t(x) %*% (ww * x)) %*% t(x) %*% (ww * y)
    cha <- abs(thetanew - theta)
    tol <- sum(cha ^ 2)
    theta <- thetanew
    iter <- iter + 1
    if (iter >= 200)
      break
  }
  res <- y - x %*% theta
  #loop end
  return(list(
    iter = iter,
    theta = theta,
    tol = tol,
    res = res
  ))
}


### Description:
###   Functions for producing data
### Arguments:
###   n: number of observations;
###   p: dimension of observations;
###   K: the number of machines/blocks;
###   fit.cv: index of data for different machines;
###   case_beta: production methods for real model parameters. Possible values 
###              are "uniform" (default) or "var_sele"
###   case_x: the way covariates are generated. Default generation from a 
###             multivariate normal distribution;
###   hete_level: level of data heterogeneity across machines.Defaults to 1;
###   case_error: Type of error generation.Defaults to standard normal distribution.
data.gene <- function(n=10^4, p=8, K=5, fit.cv, case_beta="uniform", case_x=1, hete_level=1, case_error=1){
  switch (case_beta,
          uniform = beta0 <- runif((p+1), 0, 1) +1 ,
          var_sele = beta0 <- c(2,7,1.5,0,0,3,rep(0,p-5)) #p must be equal to or greater than 5
  )
  
  switch(case_x,
         { 
           ## X generate from a multivariate normal distribution
           ind <- 1:p
           sigmax  <- 0.5 ^abs(outer(ind, ind, "-"))
           X  <-Rfast::rmvnorm(n, rep(0, p), sigmax)     #sigmax:The covariance matrix of X
         },{ 
           ## For producing heterogeneous data
           rho_k <- runif(K,0.2,0.3)
           X <- matrix(0,n,p)
           ind <- 1:p
           nk <- tabulate(fit.cv$which)
           for (k in 1:K) {
             mean_k <- hete_level*(k/K)
             sigmax_k  <- rho_k[k] ^ abs(outer(ind, ind, "-"))
             index.k <- fit.cv$subsets[fit.cv$which==k]
             X[index.k,] <- MASS::mvrnorm(nk[k],rep(mean_k,p),sigmax_k)
             
           }
         }
  )
  
  switch(case_error,
         {   ## case_error 1:  normal mean 0
           Y<-  beta0[1] + X %*% beta0[-1] + rnorm(n)
         },{ ## case_error 2: t3
           Y<-  beta0[1] + X %*% beta0[-1] + rt(n, df=3)
         },{ ## case_error 3: standard Cauchy distribution
           Y<-  beta0[1] + X %*% beta0[-1] + rcauchy(n)
         },{ ## case_error 4：Laplace distribution
           Y<-  beta0[1] + X %*% beta0[-1] + rlaplace(n)
         })
  
  return(list(y=Y, x=X, beta=beta0))
}


### Description:
###   The robust distributed least squares approximation method.
### Arguments: 
###   x: the covariate matrix; y: the response variable; 
###   fit.cv: index of data for different machines;
###   K: the number of machines/blocks; h: the bandwidth;
###   tolset: criteria for stopping algorithm convergence; 
###   case_iter: fixed rounds of communication or making the RDLSA algorithm 
###               converge.Defaults to fixed rounds of communication;
###   estep: the nummber of rounds of communication.Defaults to three round; 
###   Intercept: whether the model contains an intercept term.Defaults to FALSE.
RDLSA <- function(x, y, fit.cv, K=10, h=0.5, tolset=10^(-6), case_iter=1, estep=3, Intercept=FALSE){
  time_begin_dmlsa <- Sys.time()
  if(Intercept==TRUE) x <- cbind(1, x)
  N <- length(y)
  p <- dim(x)[2]
  
  ## One round communication. 
  Sig_h <- matrix(0,p,p)
  m_h <- matrix(0,p,1)
  
  ##### calculating the second-order gradient of G_k
  for (k in 1:K) {
    index <- fit.cv$subsets[fit.cv$which==k] # Index of data on the K-th machine.
    yk <- y[index]
    xk <- x[index,]
    nk <- length(yk)
    ak <- nk/N
    
    # Initial value of modal regression on the k-th machine
    fitk <- lm(yk~xk-1)
    betak <- coef(fitk)
    
    # MR estimator on the k-th machine
    beta_kh <- modreg.MEM(theta=betak,y=yk,x=xk,h=h,tolset= tolset,Intercept=FALSE)$theta
    
    Sig_kh <- -(grad12(xk=xk,yk=yk,beta0=beta_kh,h=h)$grad2Qk)
    Sig_h <- Sig_h + ak*Sig_kh
    m_h <- m_h + ak * Sig_kh %*% beta_kh
  }
  beta_w <- limSolve::Solve(Sig_h) %*% m_h
  
  beta_one <- c(beta_w)   # The robust weighted least squares estimator(RWLSE) 
  time_end_one <- Sys.time()
  
  
  # Multi-round communication
  tol <- 1
  iter <- 1
  switch (case_iter,
          {# To obtain the estep-RWLSE
            time_end_estep <- rep(0,(estep-1))
            beta_estep <- matrix(0,p,(estep-1))
            for (i in 1:(estep-1)) {
              Sig_h2 <- matrix(0,p,p)
              m_h2 <- matrix(0,p,1)
              for (k2 in 1:K) {
                index <- fit.cv$subsets[fit.cv$which==k2]
                yk <- y[index]
                xk <- x[index,]
                nk <- length(yk)
                ak <- nk/N
                
                fit_x1 <- grad12(xk=xk, yk=yk, beta0=beta_w, h=h)
                grad2_xk_int <- -(fit_x1$grad2Qk)             
                grad1_xk_int <- -(fit_x1$grad1Qk) 
                beta_kh2 <- beta_w - limSolve::Solve(grad2_xk_int) %*% grad1_xk_int
                Sig_h2 <- Sig_h2 + ak*grad2_xk_int
                m_h2 <- m_h2 + ak*grad2_xk_int %*% beta_kh2
              }
              beta_new <- limSolve::Solve(Sig_h2) %*% m_h2
              beta_w <- beta_new
              iter<- iter + 1
              beta_estep[,i] <- beta_new
              time_end_estep[i] <- Sys.time()
            }
          },{# Making the DMLSA algorithm converge
            beta_estep <- NULL
            while (tol > tolset) {
              Sig_h2 <- matrix(0,p,p)
              m_h2 <- matrix(0,p,1)
              for (k2 in 1:K) {
                index <- fit.cv$subsets[fit.cv$which==k2]
                yk <- y[index]
                xk <- x[index,]
                nk <- length(yk)
                ak <- nk/N
                
                fit_x1 <- grad12(xk=xk, yk=yk, beta0=beta_w, h=h)
                grad2_xk_int <- -(fit_x1$grad2Qk)             
                grad1_xk_int <- -(fit_x1$grad1Qk) 
                beta_kh2 <- beta_w - limSolve::Solve(grad2_xk_int) %*% grad1_xk_int
                Sig_h2 <- Sig_h2 + ak*grad2_xk_int
                m_h2 <- m_h2 + ak*grad2_xk_int %*% beta_kh2
              }
              beta_new <- limSolve::Solve(Sig_h2) %*% m_h2
              
              tol <- sum((beta_new - beta_w) ^ 2)
              beta_w <- beta_new
              iter<- iter + 1
              beta_estep <- beta_new
              if(iter >=200) break
            }
            time_end_estep <- Sys.time()
          }
  )
  
  time_end <- Sys.time()
  time_one <- time_end_one - time_begin_dmlsa
  time_iter <- as.POSIXlt(time_end_estep) - time_begin_dmlsa
  time_all <- time_end - time_begin_dmlsa
  return(list(coef_one=beta_one,beta=beta_estep,coef=beta_w,iter=iter,h_opt=h,time_one=time_one,time_iter=time_iter,time_all=time_all))
}


### Description:
###   The communication-efficient distributed modal regression method.
### Arguments: 
###   x: the covariate matrix; y: the response variable; 
###   fit.cv: index of data for different machines;
###   K: the number of machines/blocks; h: the bandwidth;
###   tolset: criteria for stopping algorithm convergence; 
###   case_iter: fixed rounds of communication or making the RDLSA algorithm 
###               converge.Defaults to fixed rounds of communication;
###   round_num: the nummber of rounds of communication.Defaults to one round; 
###   Intercept: whether the model contains an intercept term.Defaults to FALSE.
CDMR <- function(x, y, fit.cv, K=10, h=0.5, tolset=10^(-6), case_iter=1, round_num=1, Intercept=FALSE){
  time_begin_cdmr <- Sys.time()
  if(Intercept==TRUE) x <- cbind(1, x)
  n <- length(y)
  p <- dim(x)[2]  
  
  ## Getting the index of the data in the first machine
  index.1 <- fit.cv$subsets[fit.cv$which==1]
  y1 <- y[index.1]
  x1 <- x[index.1, ]
  
  ## Getting the initial value
  fit0 <- quantreg::rq(y1~x1-1,tau=0.5)
  beta0 <- coef(fit0)
  
  switch (case_iter,
          {#### Fixed rounds of communication.
            for (i in 1:round_num) {
              
              ## calculating intermediate quantities for renewing
              fit_x1 <- grad12(xk=x1, yk=y1, beta0=beta0, h=h)
              grad2_x1_int <- fit_x1$grad2Qk             
              grad1_x1_int <- fit_x1$grad1Qk      
              grad1_full_int <- grad1_x1_int
              for(k in 2:K){
                index.k <- fit.cv$subsets[fit.cv$which==k]    
                yk <- y[index.k]               
                xk <- x[index.k, ]
                grad1_full_int <- grad12(xk=xk, yk=yk, beta0=beta0, h=h)$grad1Qk + grad1_full_int
              }
              grad1_full_int <- grad1_full_int / K           
              grad1_1_N_int <- grad1_x1_int - grad1_full_int
              
              ## updata beta
              betanew <-  beta0 - limSolve::Solve(grad2_x1_int) %*%  (grad1_x1_int - grad1_1_N_int )
              
              beta0 <- betanew
              iter <- i
            }
          },{#### Making the CDMR algorithm converge
            tol <- 1
            iter <- 0
            while (tol > tolset) {
              
              ## calculating intermediate quantities for renewing
              fit_x1 <- grad12(xk=x1, yk=y1, beta0=beta0, h=h)
              grad2_x1_int <- fit_x1$grad2Qk            
              grad1_x1_int <- fit_x1$grad1Qk      
              grad1_full_int <- grad1_x1_int
              for(k in 2:K){
                index.k <- fit.cv$subsets[fit.cv$which==k]    
                yk <- y[index.k]               
                xk <- x[index.k, ]
                grad1_full_int <- grad12(xk=xk, yk=yk, beta0=beta0, h=h)$grad1Qk + grad1_full_int
              }
              grad1_full_int <- grad1_full_int / K           
              grad1_1_N_int <- grad1_x1_int - grad1_full_int
              
              ## updata beta
              betanew <-  beta0 - limSolve::Solve(grad2_x1_int) %*%  (grad1_x1_int - grad1_1_N_int )
              
              tol <- sum((betanew - beta0) ^ 2)
              beta0 <- betanew
              iter<- iter + 1
              if(iter >=200) break
            }
          }
  )
  time_end_cdmr <- Sys.time()
  time_iter <- time_end_cdmr - time_begin_cdmr
  return(list(coefficients=beta0, h_opt=h,time_iter=time_iter,iter=iter))
}


### Description:
###   Estimators of parameters generated using both CDMR and RDLSA methods
### Arguments:
###   n: an integer giving the number of observations/samples;
###   p: dimension of observations;
###   K: the number of machines/blocks;
###   case_x: the way covariates are generated. Default generation from a 
###             multivariate normal distribution;
###   hete_level: level of data heterogeneity across machines.Defaults to 1;
###   case_error: Type of error generation.Defaults to standard normal distribution;
###   fit.cv: index of data for different machines.
para.esti <- function(n,p=8,K=10,case_x=1,hete_level=1,case_error=1,fit.cv){
  
  ASE_CDMR1 <- rep(0,p+3)
  ASE_CDMR2 <- rep(0,p+3)
  ASE_CDMR3 <- rep(0,p+3)
  ASE_RDLSA1 <- rep(0,p+3)
  ASE_RDLSA2 <- rep(0,p+3)
  ASE_RDLSA3 <- rep(0,p+3)
  
  data <- data.gene(n=n,p=p,K=K,fit.cv=fit.cv,case_beta="uniform",case_x=case_x,hete_level=hete_level,case_error=case_error)
  x <- data$x
  y <- data$y
  beta_true <- data$beta
  h.opt1 <- select_h(x=x,y=y,fit.cv=fit.cv,K=K,Intercept = TRUE)$h.opt
  
  fit_CDMR1 <- CDMR(x=x,y=y,fit.cv=fit.cv,K=K,h=h.opt1,case_iter=1,round_num=1,Intercept = TRUE)

  fit_CDMR2 <- CDMR(x=x,y=y,fit.cv=fit.cv,K=K,h=h.opt1,case_iter=1,round_num=2,Intercept = TRUE)

  fit_CDMR3 <- CDMR(x=x,y=y,fit.cv=fit.cv,K=K,h=h.opt1,case_iter=1,round_num=3,Intercept = TRUE)

  fit_RDLSA <- RDLSA(x=x,y=y,fit.cv=fit.cv,K=K,h=h.opt1,case_iter=1,estep=3,Intercept = TRUE)

  
  ASE_RDLSA1[1:(p+1)] <- abs(fit_RDLSA$coef_one - beta_true)
  ASE_RDLSA2[1:(p+1)] <- abs(fit_RDLSA$beta[,1] - beta_true)
  ASE_RDLSA3[1:(p+1)] <- abs(fit_RDLSA$beta[,2] - beta_true)
  ASE_CDMR1[1:(p+1)] <- abs(fit_CDMR1$coefficients - beta_true)
  ASE_CDMR2[1:(p+1)] <- abs(fit_CDMR2$coefficients - beta_true)
  ASE_CDMR3[1:(p+1)] <- abs(fit_CDMR3$coefficients - beta_true)

  ##  the average estimation error (AEE)
  ASE_RDLSA1[p+2] <- mean(abs(fit_RDLSA$coef_one - beta_true))
  ASE_RDLSA2[p+2] <- mean(abs(fit_RDLSA$beta[,1] - beta_true))
  ASE_RDLSA3[p+2] <- mean(abs(fit_RDLSA$beta[,2] - beta_true))
  ASE_CDMR1[p+2] <- mean(abs(fit_CDMR1$coefficients - beta_true))
  ASE_CDMR2[p+2] <- mean(abs(fit_CDMR2$coefficients - beta_true))
  ASE_CDMR3[p+2] <- mean(abs(fit_CDMR3$coefficients - beta_true))

  RDLSA_h.opt <- fit_RDLSA$h_opt
  
  return(c(ASE_RDLSA1,ASE_RDLSA2,ASE_CDMR1,ASE_CDMR2,ASE_CDMR3,ASE_RDLSA3,RDLSA_h.opt))

}

