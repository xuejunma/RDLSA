

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
###   K: the number of machines;
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
###   Generate a sequence of λ values
### Arguments:
###   x: the covariate matrix; y: the response variable;
###   n: number of observations/samples;
###   p: dimension of covariate vectors;
###   g: the number of λ values.
f_lam <- function(x,y,n=1000,p,g=100){
  lam_v <- NULL
  lam_max <- rep(0,p)
  for (j in 1:p) {
    lam_max[j] <- abs(2*t(x[,j])%*%y)/n
  }
  lam_1 <- 1/n
  lam_g <- lam_max[which.max(lam_max)]
  c <- (log(lam_g)-log(lam_1))/(g-1)
  #lam_v[1] <- lam_1
  lam_v[g] <- lam_g
  for (i in g:2) {
    lam_v[i-1] <- lam_v[i]*exp(-c)
  }
  return(lam_v)
}


### Description:
###   a indicator function, when |a| <= b it takes the value 1 and 0 otherwise.
Indicator <- function(v,lam){
  ifelse( abs(v) <= lam, 1, 0)
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
###   hete_level: Level of data heterogeneity across machines.Defaults to 1;
###   case_error: Type of error generation.Defaults to standard normal distribution.
data.gene <- function(n=10^4, p=8, K=10, fit.cv, case_beta="uniform", case_x=1, hete_level=1, case_error=1){
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
         {    ## case_error 1:  normal mean 0
           Y<- beta0[1] + X %*% beta0[-1] + rnorm(n)
         }, { ## case_error 2: t3
           Y<-  beta0[1] + X %*% beta0[-1] + rt(n, df=3)
         }, { ## case_error 3: standard Cauchy distribution
           Y<-  beta0[1] + X %*% beta0[-1] + rcauchy(n)
         },{  ## case_error 4：Laplace distribution
           Y<-  beta0[1] + X %*% beta0[-1] + rlaplace(n)
         })
  
  return(list(y=Y, x=X, beta=beta0))
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
###   The distributed variable selection procedure based on RDLSA method and SCAD technique 
### Arguments:
###   x: the covariate matrix; y: the response variable; 
###   fit.cv: index of data for different machines;
###   K: the number of machines/blocks; h: the bandwidth;
###   lam: a nonnegative tuning parameter;
###   thre: threshold for coefficient compression to zero;
###   Intercept: whether the model contains an intercept term.Defaults to FALSE.
RDLSA_SCAD <- function(x,y,K=10,fit.cv,h=0.5,lam=0.1,thre=0.001,Intercept=FALSE){
  if(Intercept==TRUE) {
    x <- cbind(1, x)
    p <- dim(x)[2]
    N <- dim(x)[1]
    x_colname <- colnames(x)
    if(is.null(x_colname)==TRUE){
      x_colname <- c('Intercept',paste('x',1:(p-1),sep = ""))
    }
  }else{
    p <- dim(x)[2]
    N <- dim(x)[1]
    x_colname <- colnames(x)
    if(is.null(x_colname)==TRUE){
      x_colname <- paste('x',1:p,sep = "")
    }
  }
  
  ####################################  Start of initial value
  Sig_h <- matrix(0,p,p)
  m_h <- matrix(0,p,1)
  
  ## calculating the second-order gradient of G_k 
  for (k in 1:K) {
    index <- fit.cv$subsets[fit.cv$which==k]
    yk <- y[index]
    xk <- x[index,]
    nk <- length(yk)
    ak <- nk/N
    fitk <- lm(yk~xk-1)
    betak <- coef(fitk)
    beta_kh <- modreg.MEM(theta=betak,y=yk,x=xk,h=h,tolset= 10^(-6),Intercept=FALSE)$theta
    
    Sig_kh <- -(grad12(xk=xk,yk=yk,beta0=beta_kh,h=h)$grad2Qk)
    Sig_h <- Sig_h + ak*Sig_kh
    m_h <- m_h + ak * Sig_kh %*% beta_kh
  }
  
  ## calculating the initial value
  beta_w <- limSolve::Solve(Sig_h) %*% (m_h)
  beta_one <- c(beta_w)
  ####################################  End  of initial value
  
  
  ################################################# Variable selection begin
  a<-3.7
  beta_0 <- beta_one
  
  ## First derivative of the SCAD penalty function
  grad1_p <- lam*(Indicator(v=beta_0,lam = lam) + pmax((a*lam-abs(beta_0)),0)
                  /((a-1)*lam) * (1-Indicator(v=beta_0,lam = lam)))
  
  U_beta <- diag(grad1_p/abs(beta_0))
  
  ## updata beta
  sbeta_w <- limSolve::Solve(Sig_h+U_beta) %*% (m_h)
  sbeta_one <- c(sbeta_w)
  
  ## coefficient compression
  beta_1thre <- sbeta_one
  beta_1thre[which(abs(sbeta_one) <= thre)] <- 0
  beta_1full <- c(beta_1thre)
  names(beta_1full) <- c(x_colname)
  ################################################# Variable selection end 
  
  
  ## calculating the value of DBIC
  beta_cha <- beta_1full - beta_one
  DBIC_my <- t(beta_cha)%*%Sig_h%*%beta_cha + log(N)*sum(beta_1thre!=0)
  
  
  return(list(beta_one = beta_one, sbeta = beta_1full, DBIC = DBIC_my))
}


### Description:
###   The distributed variable selection procedure based on CDMR method and SCAD technique 
### Arguments:
###   x: the covariate matrix; y: the response variable; 
###   fit.cv: index of data for different machines;
###   K: the number of machines/blocks; h: the bandwidth;
###   lam: a nonnegative tuning parameter;
###   thre: threshold for coefficient compression to zero;
###   Intercept: whether the model contains an intercept term.Defaults to FALSE.
CDMR_SCAD <- function(x,y,K=10,fit.cv,h=0.5,lam=0.1,thre=0.001,Intercept=FALSE){
  if(Intercept==TRUE) {
    x <- cbind(1, x)
    p <- dim(x)[2]
    n <- dim(x)[1]
    x_colname <- colnames(x)
    if(is.null(x_colname)==TRUE){
      x_colname <- c('Intercept',paste('x',1:(p-1),sep = ""))
    }
  }else{
    p <- dim(x)[2]
    n <- dim(x)[1]
    x_colname <- colnames(x)
    if(is.null(x_colname)==TRUE){
      x_colname <- paste('x',1:p,sep = "")
    }
  }
  
  ## Getting the index of the data in the first machine
  index.1 <- fit.cv$subsets[fit.cv$which==1]
  y1 <- y[index.1]
  x1 <- x[index.1, ]
  
  ## Getting the initial value
  beta0 <- coef(quantreg::rq(y1~x1-1,tau=0.5))
  betam <- beta0
  
  ####################################### Variable selection begin
  a <- 3.7
  
  ## calculating intermediate quantities for renewing
  fit_x1 <- grad12(xk=x1, yk=y1, beta0=betam, h=h)
  grad2_x1_int <- fit_x1$grad2Qk             
  grad1_x1_int <- fit_x1$grad1Qk      
  grad1_full_int <- grad1_x1_int
  for(k in 2:K){
    index.k <- fit.cv$subsets[fit.cv$which==k]    
    yk <- y[index.k]               
    xk <- x[index.k, ]
    grad1_full_int <- grad12(xk=xk, yk=yk, beta0=betam, h=h)$grad1Qk + grad1_full_int
  }
  grad1_full_int <- grad1_full_int / K           
  grad1_1_N_int <- grad1_x1_int - grad1_full_int

  
  for (i in 1:1) {
    ## update the gradients in the first machine
    grad12_update <- grad12(xk=x1, yk=y1, beta0=betam, h=h)
    grad2_x1_update <- grad12_update$grad2Qk          
    grad1_x1_update <- grad12_update$grad1Qk
    
    ## calculating the value of the first derivative of the SCAD penalty function
    grad1_p <- lam*(Indicator(v=betam,lam = lam) + pmax((a*lam-abs(betam)),0)
                    /((a-1)*lam) * (1-Indicator(v=betam,lam = lam)))
    
    ## calculating the value of U(beta)
    U_beta <- diag(grad1_p/abs(betam))
    
    ## updata beta
    beta_new <- betam - limSolve::Solve(grad2_x1_update + U_beta)%*%
      (U_beta%*%betam + (grad1_x1_update-grad1_1_N_int))
    
    betam <- beta_new
    
  }
  
  ## coefficient compression
  beta_thre <- betam
  beta_thre[which(abs(betam) <= thre)] <- 0
  beta_full <- c(beta_thre)
  names(beta_full) <- c(x_colname)
  ####################################### Variable selection end
  
  
  ## calculating the value of BIC
  res <- y1 - x1%*%beta_thre
  Q_1h <- mean(dnorm(res/h)/h)
  Q_Nh_til <- Q_1h - t(grad1_1_N_int)%*%beta_thre
  BIC_my <- -Q_Nh_til + log(n)/n*sum(beta_thre!=0)
  
  
  return(list(beta=beta_full,BIC=BIC_my))
  
}


### Description:
###   Variable selection use both CDMR_SCAD and RDLSA_SCAD methods
### Arguments:
###   n: number of observations;
###   p: dimension of observations;
###   K: the number of machines;
###   case_x: the way covariates are generated. Default generation from a 
###             multivariate normal distribution;
###   hete_level: Level of data heterogeneity across machines.Defaults to 1;
###   case_error: Type of error generation.Defaults to standard normal distribution;
###   fit.cv: index of data for different machines.
var_sele <- function(n,p,K=10,case_x=1,hete_level=1,case_error=1,fit.cv){
  ASE_CDMR_SCAD <- rep(0,p+5)
  ASE_RDLSA_SCAD <- rep(0,p+5)
  
  data <- data.gene(n=n,p=p,K=K,fit.cv=fit.cv,case_beta="var_sele",case_x=case_x,hete_level=hete_level,case_error=case_error)
  x <- data$x
  y <- data$y
  beta_true <- data$beta
  MT <- which(beta_true!=0) ## Number of non-zero elements in beta_true
  M0 <- which(beta_true==0) ## Number of zero elements in beta_true
  
  h.opt1 <- select_h(x=x,y=y,fit.cv=fit.cv,K=K,Intercept = TRUE)$h.opt
  lam_v <- f_lam(x=x,y=y,n=n,p=p)
  CDMR_BIC <- matrix(0,length(lam_v),(p+2))
  RDLSA_DBIC <- matrix(0,length(lam_v),(p+2))
  iter<- 1
  for (  i in lam_v) {
    fit_CDMR_SCAD <- CDMR_SCAD(x=x,y=y,fit.cv=fit.cv,K=K,h=h.opt1,lam=i,Intercept = TRUE)
    fit_RDLSA_SCAD <- RDLSA_SCAD(x=x,y=y,fit.cv=fit.cv,K=K,h=h.opt1,lam=i,Intercept = TRUE)
    CDMR_BIC[iter,] <- c(fit_CDMR_SCAD$beta,fit_CDMR_SCAD$BIC)
    RDLSA_DBIC[iter,] <- c(fit_RDLSA_SCAD$sbeta,fit_RDLSA_SCAD$DBIC)
    
    iter<- iter + 1
  }
  
  ## select the estimator corresponding to the smallest BIC
  CDMR_fit <- c(CDMR_BIC[order(CDMR_BIC[,(p+2)]),][1,(1:(p+1))])
  
  ## select the estimator corresponding to the smallest DBIC
  RDLSA_fit <- c(RDLSA_DBIC[order(RDLSA_DBIC[,(p+2)]),][1,(1:(p+1))])
  
  ## calculating the values of supp(θ) for CDMR_SCAD and RDLSA_SCAD methods
  CDMR_MThat <- which(CDMR_fit != 0)
  RDLSA_MThat <- which(RDLSA_fit != 0)
  
  ## Number of zero elements in the estimators of CDMR_SCAD and RDLSA_SCAD methods
  CDMR_M0hat <- which(CDMR_fit == 0)
  RDLSA_M0hat <- which(RDLSA_fit == 0)
  
  
  ASE_CDMR_SCAD[1:(p+1)] <- abs(CDMR_fit - beta_true)
  ASE_RDLSA_SCAD[1:(p+1)] <- abs(RDLSA_fit - beta_true)
  
  ## calculating the values of AEE for CDMR_SCAD and RDLSA_SCAD methods
  ASE_CDMR_SCAD[(p+2)] <- mean(abs(CDMR_fit - beta_true))
  ASE_RDLSA_SCAD[(p+2)] <- mean(abs(RDLSA_fit - beta_true))
  
  ## calculating the values of FP for CDMR_SCAD and RDLSA_SCAD methods
  ASE_CDMR_SCAD[(p+3)] <- length(CDMR_MThat) - length(intersect(CDMR_MThat,MT))
  ASE_RDLSA_SCAD[(p+3)] <- length(RDLSA_MThat) - length(intersect(RDLSA_MThat,MT))
  
  ## calculating the values of FN for CDMR_SCAD and RDLSA_SCAD methods
  ASE_CDMR_SCAD[(p+4)] <- length(CDMR_M0hat) - length(intersect(CDMR_M0hat,M0))
  ASE_RDLSA_SCAD[(p+4)] <- length(RDLSA_M0hat) - length(intersect(RDLSA_M0hat,M0))
  
  ## calculating the values of CM for CDMR_SCAD and RDLSA_SCAD methods
  ASE_CDMR_SCAD[(p+5)] <- ifelse( setequal(CDMR_MThat,MT), 1, 0)
  ASE_RDLSA_SCAD[(p+5)] <- ifelse( setequal(RDLSA_MThat,MT), 1, 0)
  
  return(c(ASE_CDMR_SCAD,ASE_RDLSA_SCAD))
}

