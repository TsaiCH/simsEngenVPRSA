### This model is the baseline model for all other scenarios of community dynamics
tot.sim <- 1 # the simulaiton number (100 times in the manuscript)
source('SAD_pBoot_functions_zeta_final.R')
require(compiler)
require(poilog)
require(MASS)
cmp_pBootSAD_fn <- cmpfun(pBootSAD_fn)
cmp_PoiLogFit_fn <- cmpfun(PoiLogFit_fn)
#
print(Sys.time())
sigma.E <- seq(0.05,1,length.out=40) # magnitude in environmental variance varies across 40 reefs
nsp <- 100 #true species richness
fit.data.all <- A.allsim <- E.allsim <- D.allsim <- c()
for(k in 1:length(sigma.E)){
nrepeat <- 1
repeat{
  print(paste("Sim Reef=",k,sep=''))
  repeat
  {
    alpha <- 0 # no inter-specific density dependence
    B <- matrix(1,nsp,nsp)*alpha
    diag(B) <- 0.84  # strength of intra-specific density dependence = 1-0.84 = 0.16
    eig <- eigen(B)$value
    unit <- Re(eig)^2+Im(eig)^2
    if(all(unit<0.99))
    {
      break
    }
  }
  A <- rnorm(nsp,1.5,0.25)   # normally distributed intrinsic growth rate on log scale (intrinsic species differences)
  X <- rep(10,nsp)     # # initialise species abundances (doesn't matter for the stationary case)
  X.all <- E.all <- D.all <- c()
  for(i in 1:1000)
  {
    varSigma <- rep(0.5,nsp) # no species differences in environmental variance
    D <- matrix(rnorm(nsp,0,sqrt(0.5)),nsp,1)/sqrt(exp(X)) # Demographic stochasticity
    E <- matrix(mvrnorm(1,rep(0,nsp),diag(varSigma*sigma.E[k])),nsp,1) # Environmental stochasticity  
    X <- A+B%*%X+E+D
    E.all <- cbind(E.all,E)
    D.all <- cbind(D.all,D)
    X.all <- cbind(X.all,X)
  }
  rSAD_gompertz <- tail(t(X.all),11) # time-series time window = 11 years (comparable to LTMP data), simulated 'true' species abundance
  SAD_time <- cmp_pBootSAD_fn(rSAD_gompertz,totalN=1500)$SAD_time # Poisson sampling from simulated 'true' species abundance
  SAD_time[is.na(SAD_time)] <- 0 
  nspcheck <- apply(SAD_time,2,function(x){sum(length(which(x!=0)))})
  nrepeat <- nrepeat+1
  if(all(nspcheck>(nsp*0.4))){
    print(paste("Observed (sampled) richness > 40 for 11 years",sep=""))
    break 
  }
  else if(nrepeat>30){
    print(paste("Warning: nrepeat >30 so that R is terminated",sep=""))
    break
  }
  else{print(paste("Warning: sampling intensity should be increased !!! nrepeat = ",nrepeat,sep=""))}
}
  if(nrepeat>30){
    fitout <- cmp_PoiLogFit_fn(NA)
  }
  else{
    fitout <- cmp_PoiLogFit_fn(SAD_time)
  }
  A.allsim <- cbind(A.allsim,A)
  E.allsim <- cbind(E.allsim,as.numeric(E.all))
  D.allsim <- cbind(D.allsim,as.numeric(D.all))
  fit <- fitout$rho_list
  fit.data <- data.frame(rho=fit$rho,lag=fit$lag,reefID=paste("SimReef",k,sep=''),sigma=sigma.E[k],sig1=fit$sig1,sig2=fit$sig2,strue1=fit$s.true1,strue2=fit$s.true2,nsize=fit$nsize,mu1=fit$mu1,mu2=fit$mu2)
  fit.data.all <- rbind(fit.data.all,fit.data)
}
print(Sys.time())
#
filename <- paste0("~/simdata",tot.sim,".csv")
write.csv(fit.data.all,filename)
