RDLSA_SCAD_Code.R has ten functions: grad12, chunk.f, select_h, f_lam, Indicator,
data.gene, modreg.MEM, RDLSA_SCAD, CDMR_SCAD, var_sele.


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


f_lam

Description:
  Generate a sequence of λ values
  
Usage
  f_lam(x,y,n=1000,p,g=100)

Arguments:
  x: the covariate matrix; y: the response variable;
  n: number of observations/samples;
  p: dimension of covariate vectors;
  g: the number of λ values.


Indicator

Description:
  a indicator function, when |a| <= b it takes the value 1 and 0 otherwise.
  
Usage
  Indicator(v,lam)
  
Arguments:
  v: a variable;
  lam: a variable.


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


RDLSA_SCAD

Description:
  The distributed variable selection procedure based on RDLSA method and SCAD technique 

Usage
  RDLSA_SCAD(x,y,K=10,fit.cv,h=0.5,lam=0.1,thre=0.001,Intercept=FALSE)

Arguments:
  x: the covariate matrix;
  y: the response variable; 
  fit.cv: index of data for different machines;
  K: the number of machines/blocks;
  h: the bandwidth;
  lam: a nonnegative tuning parameter;
  thre: threshold for coefficient compression to zero;
  Intercept: whether the model contains an intercept term.Defaults to FALSE.


CDMR_SCAD

Description:
  The distributed variable selection procedure based on CDMR method and SCAD technique

Usage
  CDMR_SCAD(x,y,K=10,fit.cv,h=0.5,lam=0.1,thre=0.001,Intercept=FALSE)

Arguments:
  x: the covariate matrix;
  y: the response variable; 
  fit.cv: index of data for different machines;
  K: the number of machines/blocks;
  h: the bandwidth;
  lam: a nonnegative tuning parameter;
  thre: threshold for coefficient compression to zero;
  Intercept: whether the model contains an intercept term.Defaults to FALSE.


var_sele

Description:
  Variable selection use both CDMR_SCAD and RDLSA_SCAD methods

Usage
  var_sele(n,p,K=10,case_x=1,hete_level=1,case_error=1,fit.cv)

Arguments:
  n: number of observations;
  p: dimension of observations;
  K: the number of machines;
  case_x: the way covariates are generated. Default generation from a 
            multivariate normal distribution;
  hete_level: Level of data heterogeneity across machines.Defaults to 1;
  case_error: Type of error generation.Defaults to standard normal distribution;
  fit.cv: index of data for different machines.
