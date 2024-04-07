RDLSA_Code.R has eight functions: grad12, chunk.f, select_h, modreg.MEM,
data.gene, RDLSA, CDMR, para.esti.


grad12

Description
  Calculating the first/second order gradients.
  
Usage
  grad12(xk,yk,beta0,h=0.5)
  
Arguments:
  x: the covariate matrix in K-th machines/blocks; 
  y: the response variable in K-th machines/blocks;
  beta0: the estimators of model parameters;
  h: the bandwidth.


chunk.f

Description
  Used to divide data into K blocks.
  
Usage
  chunk.f(n=100,K=5,R=1,partition=c("stochastic","uniform"),
                    type = c("random", "consecutive", "interleaved"))
  
Arguments:
  n: an integer giving the number of observations/samples to be split into groups;
  K: split data into K block;
  R: an integer giving the number of replications for repeated K-fold cross-validation.
     This is ignored for for leave-one-out cross-validation and other non-random splits 
     of the data;
  partition: A string specifying the type of the number of samples contained 
             in each block. Possible values are "stochastic" (default) or "uniform";
  type: a character string specifying the type of folds to be generated; 
        Possible values are "random" (the default), "consecutive" or "interleaved".
        
        
select_h

Description:
  The optimal bandwidth h is selected by the grid search method.

Usage
  select_h(x, y, fit.cv, K=20, Intercept=FALSE)

Arguments:
  x: the covariate matrix; 
  y: the response variable;
  fit.cv: index of data for different machines;
  K: the number of machines/blocks;
  Intercept: whether the model contains an intercept term.Defaults to FALSE.


modreg.MEM

Description:
  MEM Algorithm for Modal Regression

Usage
  modreg.MEM(theta, y, x, h, tolset=10^(-6), Intercept=FALSE)

Arguments:
  theta: initial value of MEM Algorithm;
  x: the covariate matrix;
  y: the response variable; 
  h: the bandwidth;
  tolset: criteria for stopping algorithm convergence;
  Intercept: whether the model contains an intercept term.Defaults to FALSE.


data.gene

Description:
  Functions for producing data

Usage
  data.gene(n=10^4, p=8, K=10, fit.cv, case_beta="uniform", case_x=1,
            hete_level=1, case_error=1)

Arguments:
  n: number of observations;
  p: dimension of observations;
  K: the number of machines/blocks;
  fit.cv: index of data for different machines;
  case_beta: production methods for real model parameters. Possible values 
             are "uniform" (default) or "var_sele"
  case_x: the way covariates are generated. Default generation from a 
            multivariate normal distribution;
  hete_level: level of data heterogeneity across machines.Defaults to 1;
  case_error: Type of error generation.Defaults to standard normal distribution.


RDLSA

Description:
  The robust distributed least squares approximation method.

Usage
  RDLSA(x, y, fit.cv, K=10, h=0.5, tolset=10^(-6), case_iter=1,
        estep=3, Intercept=FALSE)

Arguments: 
  x: the covariate matrix;
  y: the response variable; 
  fit.cv: index of data for different machines;
  K: the number of machines/blocks;
  h: the bandwidth;
  tolset: criteria for stopping algorithm convergence; 
  case_iter: fixed rounds of communication or making the RDLSA algorithm 
              converge.Defaults to fixed rounds of communication;
  estep: the nummber of rounds of communication.Defaults to three round; 
  Intercept: whether the model contains an intercept term.Defaults to FALSE.


CDMR

Description:
  The communication-efficient distributed modal regression method.
  
Usage
  CDMR(x, y, fit.cv, K=10, h=0.5, tolset=10^(-6), case_iter=1,
        round_num=1, Intercept=FALSE)  

Arguments: 
  x: the covariate matrix;
  y: the response variable; 
  fit.cv: index of data for different machines;
  K: the number of machines/blocks;
  h: the bandwidth;
  tolset: criteria for stopping algorithm convergence; 
  case_iter: fixed rounds of communication or making the RDLSA algorithm 
              converge.Defaults to fixed rounds of communication;
  round_num: the nummber of rounds of communication.Defaults to one round; 
  Intercept: whether the model contains an intercept term.Defaults to FALSE.


para.esti

Description:
  Estimators of parameters generated using both CDMR and RDLSA methods
  
Usage  
  para.esti(n,p=8,K=10,case_x=1,hete_level=1,case_error=1,fit.cv)  

Arguments:
  n: an integer giving the number of observations/samples;
  p: dimension of observations;
  K: the number of machines/blocks;
  case_x: the way covariates are generated. Default generation from a 
            multivariate normal distribution;
  hete_level: level of data heterogeneity across machines.Defaults to 1;
  case_error: Type of error generation.Defaults to standard normal distribution;
  fit.cv: index of data for different machines.
