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
blockProp  <-  rdirichlet(M,rep(5*K,K))  #### emptying some blocks in certain netwokrs 
blockProp[1,] <- rep(1/K,K) 
# print(blockProp)




########### block proportions simul iidColSBM 
#blockProp <- list()
#blockProp <-  matrix(rdirichlet(1,rep(K,K)),byrow = TRUE,nrow = M, ncol = K) #### emptying some blocks in certain netwokrs 
#blockProp$col <- matrix(rdirichlet(1,rep(K,K)),byrow = TRUE,nrow = M, ncol = K)   #### emptying some blocks in certain netwokros
#print(blockProp)




############## connectivity matrix
means <- matrix(rbeta(K*K,1/1.1,1/1.1), K, K)  
connectParam <- list(mean = round(means,2))
oCol <- order(colSums(connectParam$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParam$mean[,]),decreasing = TRUE)
connectParam$mean <- connectParam$mean[oRow,oCol]


########### sizes of networks
nbNodes <- sample(10*c(6:10),M,replace = TRUE)


############################################################
##################################"""""""" SIMULATION
############################################################


mySampler <- lapply(1:M, function(m){sampleSimpleSBM(nbNodes = nbNodes[m], blockProp =  blockProp[m,],  connectParam,directed = TRUE)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})





############### init Tau

initSimple <- lapply(collecNetworks,estimateSimpleSBM)
myRef = which.max(sum(t(sapply(initSimple,function(m){m$nbBlocks}))))
collecTau <- initCollecTau(initSimple, ref = myRef)

KEstim <- initSimple[[myRef]]$nbBlocks[1]



 
###################### Estim avec K K true


estimOptions <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)



#   Making the groups correspond between networks









hyperparamPrior_iid <- setHyperparamPrior(M,K, emissionDist  = 'bernoulli' , model ='iidColSBM')
resEstim_iid  <- VBEMColSBM(collecNetworks,hyperparamPrior_iid,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='iidColSBM')
logLikMarg_iid <- computeLogLikMarg(collecNetworks, resEstim_iid$collecTau,hyperparamPrior_iid, emissionDist  = 'bernoulli' , model ='iidColSBM')

hyperparamPrior_picol <- setHyperparamPrior(M,KEstim,KEstim, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
resEstim_picol  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_picol,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
logLikMarg_picol <- computeLogLikMarg(collecNetworks, resEstim_picol$collecTau,hyperparamPrior_picol, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')

res <- c(sum(logLikMarg_iid),sum(logLikMarg_picol))
names(res) = c('iidCol','piCol')
best_model <- names(res)[which.max(res)]
print(res)
print(best_model)

if(best_model=='iidCol'){resEstim <- resEstim_iid}
if(best_model=='piCol'){resEstim <- resEstim_picol}



postMeanEstim <- list()
postMeanEstim$connectParam <- resEstim$hyperparamPost$connectParam$alpha/(resEstim$hyperparamPost$connectParam$alpha + resEstim$hyperparamPost$connectParam$beta)
postMeanEstim$blockProp <- lapply(resEstim$hyperparamPost$blockProp,function(l){normByRow(l)}) 


postMemberships <- lapply(resEstim$collecTau,function(tau){
  Z <- list()
  Z <- apply(tau,1,which.max)
  Z$col <- apply(tau$col,1,which.max)
  return(Z)
})


lapply(1:M,function(m){table(postMemberships[[m]],mySampler[[m]]$memberships)})
lapply(1:M,function(m){table(postMemberships[[m]]$col,mySampler[[m]]$memberships$col)})


lapply(1:M,function(m){table(postMemberships[[m]],initBipartite[[m]]$memberships)})
lapply(1:M,function(m){table(postMemberships[[m]]$col,initBipartite[[m]]$memberships$col)})

lapply(1:M,function(m){table(postMemberships[[m]],mySampler[[m]]$memberships)})
lapply(1:M,function(m){table(postMemberships[[m]]$col,mySampler[[m]]$memberships$col)})


 
################ sepSBM 


logLikMarg_sep <- matrix(0,M,2)
for (m in 1:M){
  print(paste0('Network ',m))
  K_m <- initBipartite[[m]]$nbBlocks[1]
  K_m <- initBipartite[[m]]$nbBlocks[2]
  tau_m <- initBipartite[[m]]$probMemberships
  priorParam_m <- setHyperparamPrior(1, K_m,K_m,emissionDist  = 'bernoulli' , model ='iidColBipartiteSBM')
  resEstim_m  <- VBEMBipartiteColSBM(list(collecNetworks[[m]]),priorParam_m,list(tau_m),estimOptions,emissionDist = 'bernoulli',model ='iidColBipartiteSBM')
  logLikMarg_sep[m,] <-  computeLogLikMarg(list(collecNetworks[[m]]), list(tau_m),priorParam_m,emissionDist = 'bernoulli',model ='iidColBipartiteSBM')
}



cbind(sum(logLikMarg),sum(logLikMarg_sep))
    

muPost <-   resEstim$postParam$connectParam$alpha/(resEstim$postParam$connectParam$beta + resEstim$postParam$connectParam$alpha)
oCol <- order(colSums(muPost),decreasing = TRUE)
oRow <- order(rowSums(muPost),decreasing = TRUE)
resEstim$postParam$connectParam$alpha <-resEstim$postParam$connectParam$alpha[oRow,oCol]
resEstim$postParam$connectParam$beta <-resEstim$postParam$connectParam$beta[oRow,oCol]

resEstim$postParam$blockProp <- resEstim$postParam$blockProp[,oRow]
resEstim$postParam$blockProp$cow <- resEstim$postParam$blockProp$cow[,oCol]



############## Plot post 
par(mfrow = c(KEstim,KEstim))
for (k in 1:K){
  for (l in 1:K){
  curve(dbeta(x,resEstim$postParam$connectParam$alpha[k,l],resEstim$postParam$connectParam$beta[k,l]),ylab = 'post')
  abline(v=connectParam$mean[k,l],col='red')
  abline(v=muPost[k,l],col='green')
  }
}
  