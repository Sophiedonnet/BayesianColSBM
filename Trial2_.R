rm(list=ls())
library(sbm)
source('utils.R')
source('VBEMBipartiteColSBM.R')
source('initializationCollecTau.R')
library(gtools)


###########################################################################################
############# simulation avec iid colSBM mais certains très petits réseaux. 
#############################################################################################


M = 6
KRow = 4
KCol = 3




########### block proportions simul iidColSBM 
blockProp <- list()
blockProp$row <-  matrix(rdirichlet(1,rep(KRow,KRow)),byrow = TRUE,nrow = M, ncol = KRow) #### emptying some blocks in certain netwokrs 
blockProp$col <- matrix(rdirichlet(1,rep(KCol,KCol)),byrow = TRUE,nrow = M, ncol = KCol)   #### emptying some blocks in certain netwokros
print(blockProp)




############## connectivity matrix
means <- matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  
connectParam <- list(mean = round(means,2))
oCol <- order(colSums(connectParam$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParam$mean[,]),decreasing = TRUE)
connectParam$mean <- connectParam$mean[oRow,oCol]


########### sizes of networks
nbNodes <- matrix(sample(10*c(6:10),M*2,replace = TRUE),M,2)
nbNodes[5,] = c(30,20)
nbNodes[6,] = c(20,30)


############################################################
##################################"""""""" SIMULATION
############################################################


mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockProp$row[m,],col=blockProp$col[m,]),  connectParam)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})





############### init Tau

initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef = which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau <- initCollecTau(initBipartite,ref = myRef)

KRowEstim <- initBipartite[[myRef]]$nbBlocks[1]
KColEstim <- initBipartite[[myRef]]$nbBlocks[2]



###################### Estim avec Kcol Krow true


estimOptions <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)






hyperparamPrior_iid <- setHyperparamPrior(M,KRowEstim,KColEstim, emissionDist  = 'bernoulli' , model ='iidColBipartiteSBM')
resEstim_iid  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_iid,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='iidColBipartiteSBM')
logLikMarg_iid <- computeLogLikMarg(collecNetworks, resEstim_iid$collecTau,hyperparamPrior_iid, emissionDist  = 'bernoulli' , model ='iidColBipartiteSBM')

hyperparamPrior_picol <- setHyperparamPrior(M,KRowEstim,KColEstim, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
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
postMeanEstim$connectParam <- resEstim$hyperparamPost$connectParam$alpha/(resEstim$postParam$connectParam$alpha + resEstim$postParam$connectParam$beta)
postMeanEstim$blockProp <- lapply(resEstim$hyperparamPost$blockProp,function(l){normByRow(l)}) 


postMemberships <- lapply(resEstim$collecTau,function(tau){
  Z <- list()
  Z$row <- apply(tau$row,1,which.max)
  Z$col <- apply(tau$col,1,which.max)
  return(Z)
  })


lapply(1:M,function(m){table(postMemberships[[m]]$row,mySampler[[m]]$memberships$row)})
lapply(1:M,function(m){table(postMemberships[[m]]$col,mySampler[[m]]$memberships$col)})


lapply(1:M,function(m){table(postMemberships[[m]]$row,initBipartite[[m]]$memberships$row)})
lapply(1:M,function(m){table(postMemberships[[m]]$col,initBipartite[[m]]$memberships$col)})




 

   

muPost <-   resEstim$postParam$connectParam$alpha/(resEstim$postParam$connectParam$beta + resEstim$postParam$connectParam$alpha)
oCol <- order(colSums(muPost),decreasing = TRUE)
oRow <- order(rowSums(muPost),decreasing = TRUE)
resEstim$postParam$connectParam$alpha <-resEstim$postParam$connectParam$alpha[oRow,oCol]
resEstim$postParam$connectParam$beta <-resEstim$postParam$connectParam$beta[oRow,oCol]

resEstim$postParam$blockProp$row <- resEstim$postParam$blockProp$row[,oRow]
resEstim$postParam$blockProp$cow <- resEstim$postParam$blockProp$cow[,oCol]



############## Plot post 
par(mfrow = c(KRow,KCol))
for (k in 1:KRow){
  for (l in 1:KCol){
  curve(dbeta(x,resEstim$postParam$connectParam$alpha[k,l],resEstim$postParam$connectParam$beta[k,l]),ylab = 'post')
  abline(v=connectParam$mean[k,l],col='red')
  abline(v=muPost[k,l],col='green')
  }
}
  