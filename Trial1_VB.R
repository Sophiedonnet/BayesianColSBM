rm(list=ls())
#------------------ set wd
if(Sys.info()[[4]]=='sophie-Latitude-5310'){
  setwd('/home/sophie/WORK_LOCAL/RECHERCHE/TRAVAUX_DE_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}
if(Sys.info()[[4]]=="donnet-Precision-Tower-5810"){
  setwd('/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}

#--------------- packages and functions
library(sbm)
library(gtools)
library(DescTools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)



#seed  <- 14579
#set.seed(seed)
M = 6
K = 4


########### block proportions simul piColSBM avec classes vides

blockProp <-  rdirichlet(M,rep(1/(K-1),K))  #### emptying some blocks in certain netwokrs 
blockProp[1,] <- rep(1/K,K) 
print(blockProp)


########### block proportions simul piColSBM avec quasi même proba partout
# blockProp  
blockProp  <-  rdirichlet(M,rep(5*K,K))  #### all with the same proba
blockProp[1,] <- rep(1/K,K) 
 print(blockProp)




########### block proportions simul iidColSBM 
#blockProp <- list()
#blockProp <-  matrix(rdirichlet(1,rep(K,K)),byrow = TRUE,nrow = M, ncol = K) #### emptying some blocks in certain netwokrs 
#blockProp$col <- matrix(rdirichlet(1,rep(K,K)),byrow = TRUE,nrow = M, ncol = K)   #### emptying some blocks in certain netwokros
#print(blockProp)




############## connectivity matrix
means <- matrix(rbeta(K*K,1/1.1,1/1.1), K, K)  
connectParam <- list(mean = round(means,2))
or <- order(diag(connectParam$mean),decreasing = TRUE)
connectParam$mean <- connectParam$mean[or,or]


########### sizes of networks
nbNodes <- sample(10*c(2:6),M,replace = TRUE)


############################################################
##################################"""""""" SIMULATION
############################################################


mySampler <- lapply(1:M, function(m){sampleSimpleSBM(nbNodes = nbNodes[m], blockProp =  blockProp[m,],  connectParam,directed = TRUE)})
collecNetworks <- lapply(mySampler,function(l){Net.m <- l$networkData; diag(Net.m) <- 0; Net.m})


 
############################################################
############### init Par Frequentist
#########################################################"
initSimple <- lapply(collecNetworks,estimateSimpleSBM)
myRef = which.max(sum(t(sapply(initSimple,function(m){m$nbBlocks}))))
collecTau <- initCollecTau(initSimple, ref = myRef)

KEstim <- initSimple[[myRef]]$nbBlocks[1]


#######################################################"
########### Estim par VBEM
######################################################"

# Estim avec K K true


estimOptions <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-10,
                     valStopCritVB = 10^-10,epsTau = 10^-10)



################# iid
hyperparamPrior <- setHyperparamPrior(M,K, emissionDist  = 'bernoulli' , model ='iidColSBM')
resEstim_iid  <- VBEMColSBM(collecNetworks,hyperparamPrior,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='iidColSBM')
logLikMarg_iid <- computeLogLikMarg_VB(collecNetworks, resEstim_iid$collecTau,hyperparamPrior, emissionDist  = 'bernoulli' , model ='iidColSBM')

hyperparamPrior_picol <- setHyperparamPrior(M,KEstim, emissionDist  = 'bernoulli' , model ='piColSBM')
resEstim_picol  <- VBEMColSBM(collecNetworks,hyperparamPrior_picol,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='piColSBM')
logLikMarg_picol <- computeLogLikMarg_VB(collecNetworks, resEstim_picol$collecTau,hyperparamPrior_picol, emissionDist  = 'bernoulli' , model ='piColSBM')

res <- c(sum(logLikMarg_iid),sum(logLikMarg_picol))
names(res) = c('iidCol','piCol')
best_model <- names(res)[which.max(res)]
print(res)
print(best_model)

if(best_model=='iidCol'){resEstim <- resEstim_iid}
if(best_model=='piCol'){resEstim <- resEstim_picol}



postMeanEstim <- list()
postMeanEstim$connectParam <- resEstim$hyperparamPost$connectParam$alpha/(resEstim$hyperparamPost$connectParam$alpha + resEstim$hyperparamPost$connectParam$beta)
postMeanEstim$blockProp <- normByRow(resEstim$hyperparamPost$blockProp)

postMemberships <- lapply(resEstim$collecTau,function(tau){
  Z <- list()
  Z <- apply(tau,1,which.max)
  return(Z)
})


lapply(1:M,function(m){table(postMemberships[[m]],mySampler[[m]]$memberships)})
lapply(1:M,function(m){table(postMemberships[[m]],initSimple[[m]]$memberships)})




 
################ sepSBM 


logLikMarg_sep <- matrix(0,M,2)
for (m in 1:M){
  print(paste0('Network ',m))
  K_m <- initSimple[[m]]$nbBlocks
  tau_m <- initSimple[[m]]$probMemberships
  priorParam_m <- setHyperparamPrior(1, K_m,emissionDist  = 'bernoulli' , model ='iidColSBM')
  resEstim_m  <- VBEMColSBM(list(collecNetworks[[m]]),priorParam_m,list(tau_m),estimOptions,emissionDist = 'bernoulli',model ='iidColSBM')
  logLikMarg_sep[m,] <-  computeLogLikMarg_VB(list(collecNetworks[[m]]), list(tau_m),priorParam_m,emissionDist = 'bernoulli',model ='iidColSBM')
}



cbind(sum(logLikMarg),sum(logLikMarg_sep))
    

muPost <-   resEstim$hyperparamPost$connectParam$alpha/(resEstim$hyperparamPost$connectParam$beta + resEstim$hyperparamPost$connectParam$alpha)

############## Plot post 
par(mfrow = c(KEstim,KEstim))
for (k in 1:KEstim){
  for (l in 1:KEstim){
  curve(dbeta(x,resEstim$hyperparamPost$connectParam$alpha[k,l],resEstim$hyperparamPost$connectParam$beta[k,l]),ylab = 'post')
  abline(v=connectParam$mean[k,l],col='red')
  abline(v=muPost[k,l],col='green')
  }
}
  
