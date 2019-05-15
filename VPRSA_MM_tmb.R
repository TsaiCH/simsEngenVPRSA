########################################################################################################################
#### The mixed-effects model (MM) approach for estimating VPRSA (cf. IML approach of Engen et al. 2002)
#### see eq. 1-4 and eq. 14
#### CAUTION: the space ahead of the line is needed for C code template before submitting to TMB compiler
########################################################################################################################
library(TMB)
library(methods)
engen.cpp <-"// Normal non-linear mixed model specified through sparse design matrices.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(rho);         // Observations
  DATA_VECTOR(lag);         // Observations
  DATA_SPARSE_MATRIX(Z);  // Random effect design matrix
  PARAMETER_VECTOR(u_a);    // Random effects vector
  PARAMETER_VECTOR(u_i);    // Random effects vector
  PARAMETER_VECTOR(u_d);    // Random effects vector
  PARAMETER(logasymp); // Fixed effects vector
  PARAMETER(loginit); // Fixed effects vector
  PARAMETER(logdelta); // Fixed effects vector
  PARAMETER(logsdu_a);      // Random effect standard deviations
  PARAMETER(logsdu_i);      // Random effect standard deviations
  PARAMETER(logsdu_d);      // Random effect standard deviations
  PARAMETER(logsd0);      // Measurement standard deviation
  // Distribution of random effect (u_a, u_i, u_d):
  Type ans = 0;
  ans -= dnorm(u_a, Type(0), exp(logsdu_a), true).sum();
  ans -= dnorm(u_i, Type(0), exp(logsdu_i), true).sum();
  ans -= dnorm(u_d, Type(0), exp(logsdu_d), true).sum();
  // Distribution of obs given random effects (x|u):
  vector<Type> Fasymp = exp(logasymp + Z*u_a);
  vector<Type> Fdelta = exp(logdelta + Z*u_d)*lag;
  vector<Type> Finit = exp(loginit + Z*u_i);
  vector<Type> y = (Finit-Fasymp)*exp(-Fdelta) + Fasymp;
  ans -= dnorm(rho, y, exp(logsd0), true).sum();
  return ans;
}"
writeLines(engen.cpp,con="engen.cpp")
compile("engen.cpp")
dyn.load(dynlib("engen"))
setwd('~/csv_sim_output/') # direct to the folder that stores csv output files of 100 simulation results for each scenario of community dynamics
csvs <- list.files()
parfix.all <- parrandom.all <- c()
for(i in 1:length(csvs)){
  print(csvs[i])
  fit.data.all <- read.csv(csvs[i])
  fit.data.all$reefID <- factor(as.character(fit.data.all$reefID))
  Z = model.matrix(~reefID-1,fit.data.all)
  Z <- as(Z,"dgTMatrix")
  simdata <- list(rho=fit.data.all$rho,lag=fit.data.all$lag,Z=Z)
  par <- list(u_a=rnorm(ncol(Z))*0,u_i=rnorm(ncol(Z))*0,u_d=rnorm(ncol(Z))*0,
              logasymp=-0.2,loginit=-0.05,logdelta=-0.1,
              logsdu_a=-3,logsdu_i=-5,logsdu_d=-3,logsd0=-2)
  obj <- MakeADFun(data=simdata,
                   parameters=par,
                   random=c("u_a","u_i","u_d"),
                   DLL="engen",
                   hessian=TRUE,
                   silent=TRUE)
  opt <- nlminb(obj$par,obj$fn,obj$gr,upper=c(-0.05,0,10,10,10,10,10),lower=c(-10,-0.05,-15,-15,-15,-15,-15))
  report <- sdreport(obj)
  parfix.all <- rbind(parfix.all, report$par.fixed)
  parrandom.all <- rbind(parrandom.all, report$par.random)
}
save(parfix.all,parrandom.all,file="~/VPRSA_MM_tmb_allsims_results.RData")
