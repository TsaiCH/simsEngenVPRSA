# varInterS scenarios of Gompertz-type communtiy dynamics
tot.sim <- 1 # simulation number (100 time in total for the manuscript)
source('SAD_pBoot_functions_zeta_final.R')
require(compiler)
require(poilog)
require(MASS)
cmp_pBootSAD_fn <- cmpfun(pBootSAD_fn )
cmp_PoiLogFit_fn <- cmpfun(PoiLogFit_fn)
#
print(Sys.time())
sigma.E <- seq(0.05,1,length.out=40) # magnitude in environmental variance varies across 40 reefs
ABE.list <- vector("list",40)
nsp <- 100 # "true" species richness
fit.data.all <- A.allsim <- E.allsim <- D.allsim <- c()
for(k in 1:length(sigma.E)){
nrepeat <- 1
repeat{
  print(paste("Sim Reef=",k,sep=''))
  repeat
  {
    alpha <- runif(nsp*(nsp-1)*0.5,0,0.002) # nsp-1 species in total will have ~60% strength of intra-specific density dependence of species i (~0.001*99/0.16)
    B <- matrix(1,nsp,nsp)
    B[upper.tri(B)] <- alpha
    B[lower.tri(B)] <- alpha
    diag(B) <- 0.84  # strength of intra-specific density dependence ~ 1-0.84 = 0.16 as per the baseline model
    eig <- eigen(B)$value
    unit <- Re(eig)^2+Im(eig)^2
    if(all(unit<0.99))
    {
      break
    }
  }
  A <- rnorm(nsp,1.5,0.25)  # normally distributed intrinsic growth rate on log scale (intrinsic species differences)
  ABE.list[[k]][[1]] <- A
  ABE.list[[k]][[2]] <- B
  varNhat <- var(solve(diag(rep(1,nrow(B)))-B)%*%A) # Variance of deterministic part equilibrium abundances (var of N*)
  varSigma <- rep(0.5,nsp)
  require(matrixcalc)
  # WARNING! the intensive computation of following environmental variance (varEnv) (see Ives et al. 2003 for the matrix solution)
  varcovEnv <- matrix(solve(diag(1,nrow(B)^2)-kronecker(B,B))%*%vec(diag(varSigma*sigma.E[k])),nsp,nsp)
  varEnv <- diag(varcovEnv) # species' environmental variances after adjusting density dependences
  varEnv.m <- mean(varEnv) # averaged (true) speceis' environmental variances 
  ABE.list[[k]][[3]] <- varcovEnv
  X <- rep(10,nsp)     # initialize species abundances (doesn't matter for the stationary case)
  X.all <- E.all <- D.all <- c()
  for(i in 1:1000)
  {
    varSigma <- rep(0.5,nsp) # no speceies differences in environmental variance as per the baseline model
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
    print(paste("Observed (sampled) richness > 40 for all years",sep=""))
    break 
  }
  else if(nrepeat>30){
    print(paste("Warning: nrepeat >30 so the loop is terminated",sep=""))
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
  fit.data <- data.frame(rho=fit$rho,lag=fit$lag,reefID=paste("SimReef",k,sep=''),sigma=sigma.E[k],sig1=fit$sig1,sig2=fit$sig2,strue1=fit$s.true1,strue2=fit$s.true2,nsize=fit$nsize,mu1=fit$mu1,mu2=fit$mu2,varNhat=varNhat,varEnvM=varEnv.m)
  fit.data.all <- rbind(fit.data.all,fit.data)
}
print(Sys.time())
#
filename <- paste0("~/simdata",tot.sim,".csv")
filename2 <- paste0("~/ABElist",tot.sim,".RData")
write.csv(fit.data.all, filename)
save(ABE.list, file=filename2)
