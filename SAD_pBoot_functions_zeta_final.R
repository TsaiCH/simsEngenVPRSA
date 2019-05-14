#####################################################################################################################
#### Function for generating Poisson sampling of species-abundance distributions (Connolly et al. 2009)
#### rSAD_gompertz: simulated "true" log-abundances of species over time
#### totalN: total number of individuals as Poisson sampling intensity
###################################################################################################################
pBootSAD_fn <- function(rSAD_gompertz,totalN=1500)
{
  SAD_time <- c()
  for(k in 1:NROW(rSAD_gompertz))
  {  
    x <- rSAD_gompertz[k,]
    rSAD <- exp(x)/sum(exp(x)) 
    nsp <- length(rSAD)
    rx <- rSAD
    splist <- c(1:nsp)
    SAD_sim <- numeric()
    N <- totalN
    for(i in 1:nsp)
    {
      SAD_sim[i] <- rpois(1,rx[i]*N)
    }
    SAD_time <- cbind(SAD_time,SAD_sim)
  }
  return(list(SAD_time=SAD_time))
}

###################################################################################################################
#### Function for estimating bivariate Poisson-lognormal species-abundance distributions by MLE (Engen et al. 2002)
#### input SADdat: species-abundance distribution over time (e.g., SAD_time from the above fn, nlag=11)
#### output rho_list: MLEs of bivariate Possion-lognormal species-abundance distributions
#### CAUTION: this function is tweek for only 11-yr time window
###################################################################################################################
PoiLogFit_fn <- function(SADdat)
{
  require(poilog) # Engen et al. 2002, 2011
  nlag <- ncol(SADdat) 
  parms1<-parms2<-parms3<-parms4<-parms5<-parms6<-parms7<-parms8<-parms9<-parms10<-c()         
  for(m in 1:nlag){
    for(n in 1:nlag){
      if(m<n)
      {
        print(paste('m=',m,', n=',n,sep=''))
        xy <- cbind(SADdat[,m],SADdat[,n])
        new.xy <- xy[rowSums(xy)!=0,]
        nsp.x <- length(which(new.xy[,1]!=0))
        nsp.y <- length(which(new.xy[,2]!=0))
        
        est <- try(bipoilogMLE(new.xy))
        if(NROW(est)!=1)
        {
          parms <- c(est$par,est$logLval,nrow(new.xy),nsp.x/est$p[1],nsp.y/est$p[2])
        }
        else
        {
          print(paste('BFGS convergence failed at:',' m=',m,', n=',n," ;Try SANN",sep=''))
          est2 <- try(bipoilogMLE(new.xy,method='SANN'))
          if(NROW(est2)!=1)
          {
            parms <- c(est2$par,est2$logLval,nrow(new.xy),nsp.x/est2$p[1],nsp.y/est2$p[2])
          }
          else
          {
            parms <- NA
          }          
        }
        if((n-m)==1){parms1 <- cbind(parms1,parms)}
        else if((n-m)==2){parms2 <- cbind(parms2,parms)}
        else if((n-m)==3){parms3 <- cbind(parms3,parms)}
        else if((n-m)==4){parms4 <- cbind(parms4,parms)}
        else if((n-m)==5){parms5 <- cbind(parms5,parms)}
        else if((n-m)==6){parms6 <- cbind(parms6,parms)}
        else if((n-m)==7){parms7 <- cbind(parms7,parms)}
        else if((n-m)==8){parms8 <- cbind(parms8,parms)}
        else if((n-m)==9){parms9 <- cbind(parms9,parms)}
        else if((n-m)==10){parms10 <- cbind(parms10,parms)}
      }        
    }
  }
  poilog.parms <- list(lag1=parms1,lag2=parms2,lag3=parms3,lag4=parms4,lag5=parms5,lag6=parms6,lag7=parms7,lag8=parms8,lag9=parms9,lag10=parms10)  
  # rho_list to collect rho and rho.boots from poilogMLE results 
  rho.all<-lag.all<-mu1.all<-mu2.all<-sig1.all<-sig2.all<-nsize.all<-s.true1.all<-s.true2.all<-AICc.all<-c()
  for(k in 1:10)
  {
    mu1 <- poilog.parms[[k]][1,]
    mu2 <- poilog.parms[[k]][2,]
    sig1 <- poilog.parms[[k]][3,]
    sig2 <- poilog.parms[[k]][4,]
    rho <- poilog.parms[[k]][5,]
    AICc <- 2*5*-2*poilog.parms[[k]][6,]+2*5*(5+1)/(poilog.parms[[k]][7,]-5-1)
    nsize <- poilog.parms[[k]][7,]
    s.true1 <- poilog.parms[[k]][8,]
    s.true2 <- poilog.parms[[k]][9,]
    lag <- rep(k,length(rho))     
    
    rho.all <- c(rho.all,rho)
    lag.all <- c(lag.all,lag)
    mu1.all <- c(mu1.all,mu1)
    mu2.all <- c(mu2.all,mu2)
    sig1.all <- c(sig1.all,sig1)
    sig2.all <- c(sig2.all,sig2)
    nsize.all <- c(nsize.all,nsize)
    s.true1.all <- c(s.true1.all,s.true1)
    s.true2.all <- c(s.true2.all,s.true2)
    AICc.all <- c(AICc.all,AICc)
  }
  rho_list <- data.frame(rho=rho.all,lag=lag.all,mu1=mu1.all,mu2=mu2.all,sig1=sig1.all,sig2=sig2.all,nsize=nsize.all,s.true1=s.true1.all,s.true2=s.true2.all,AICc=AICc.all)
  return(list(rho_list=rho_list,poilog_parms=poilog.parms))
}
#function end here



