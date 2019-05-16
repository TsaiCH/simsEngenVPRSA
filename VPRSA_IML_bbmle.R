########################################################################################################################
#### The individual maximum likelihood (IML) approach for estimating VPRSA (cf. proposed MM approach of the manuscript)
#### see eq. 1-4 and cf eq. 14 
########################################################################################################################
# call library for ML estimation
library(bbmle)
# likelihood function of VPRSA
LL <- function(loginit=log(0.95), logasymp=log(0.6), logdelta=log(0.1), logsig=log(0.01)) {
  -sum(stats::dnorm(rho, mean=(exp(loginit)-exp(logasymp))*exp(-exp(logdelta)*lag)+exp(logasymp), sd=exp(logsig), log=TRUE))
}
setwd('~/csv_sim_output/') # direct to the folder that stores csv output files of 100 simulation results for each scenario of community dynamics
csvs <- list.files()
MLEreefDATA.out <- c()
for(j in 1:length(csvs)){
  print(csvs[j])
  simID <- gsub(".csv","",gsub("simdata","",csvs[j]))
  reefDATA.out <- read.csv(csvs[j])[,-1]
  reefs <- unique(as.character(reefDATA.out$reefID))
  MLEreef.out <- c()
  for(i in 1:length(reefs)){
    dat <- reefDATA.out[reefDATA.out$reefID==reefs[i],]
    Esigma <- as.numeric(unique(dat$sigma))
    rho <- dat$rho
    lag <- dat$lag
    datt <- data.frame(cbind(rho,lag))
    fit <- mle2(rho~dnorm(mean=(exp(loginit)-exp(logasymp))*exp(-exp(logdelta)*lag)+exp(logasymp),sd=exp(logsig)),
                start=list(loginit=log(0.95),logasymp=log(0.8),logdelta=log(0.1),logsig=log(0.01)),
                data=datt,method="L-BFGS-B",upper=c(loginit=0, logasymp=0))
    res <- residuals(fit)
    est <- data.frame(t(exp(coef(fit))))
    colnames(est) <- c("init","asymp","delta","bbmlesig")
    MLE.out <- data.frame(est,reefID=reefs[i],simID=simID,Esigma=Esigma)
    MLEreef.out <- rbind(MLEreef.out,MLE.out)
  }
  MLEreefDATA.out <- rbind(MLEreefDATA.out,MLEreef.out)
}
write.csv(MLEreefDATA.out,"~/MLEreefDATA_scenarioName.csv",row.names=FALSE) # ouput the IML estimation of VPRSA for 100 simulation results of a given scenario








