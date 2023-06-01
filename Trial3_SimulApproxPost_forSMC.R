rm(list=ls())
library(sbm)
library(gtools)

source('Functions/utils.R')
source('Functions/VBEMBipartiteColSBM.R')
source('Functions/initializationCollecTau.R')
source('Functions/Functions-SMC.R')
source('Functions/Functions-colBipartiteSBM.R')
source('Functions/Functions-checks.R')
source('Functions/Function-MCMC-colBipartiteSBM.R')

emissionDist = 'bernoulli'
model = 'piColBipartiteSBM'
model = 'iidColBipartiteSBM'
###########################################################################################
############# simulation avec iid colSBM mais certains très petits réseaux. 
#############################################################################################


M = 1
KRow = 4
KCol = 3
#---  block proportions simul iidColSBM 
blockProp <- list()
blockProp$row <-  matrix(rdirichlet(1,rep(KRow,KRow)),byrow = TRUE,nrow = M, ncol = KRow) #### not emptying some blocks in certain netwokrs 
blockProp$col <- matrix(rdirichlet(1,rep(KCol,KCol)),byrow = TRUE,nrow = M, ncol = KCol)   #### not  emptying some blocks in certain netwokros
blockPropTrue <- blockProp
rm(blockProp)

#-----  connectivity matrix
means <- matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  
connectParam <- list(mean = round(means,2))
oCol <- order(colSums(connectParam$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParam$mean[,]),decreasing = TRUE)
connectParam$mean <- connectParam$mean[oRow,oCol]

connectParamTrue <- connectParam
rm(connectParam)
#-------  sizes of networks
nbNodes <- matrix(sample(10*c(6:10),M*2,replace = TRUE),M,2)
nbNodes[5,] = c(30,20)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row[m,],col=blockPropTrue$col[m,]),  connectParamTrue)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})

 
#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef <- which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]


#---------------  Prior distribution

hyperparamPrior <- setHyperparamPrior(M,KRow,KCol, emissionDist, model)
hyperparam <- hyperparamPrior
MC = 3000
HSamplePrior <- rParamZ(MC, hyperparam, emissionDist, model ,nbNodes)
checkSample(HSamplePrior,MC,M,KRow,KCol,hyperparam,emissionDist,model)


#------------- variational estim 

estimOptions <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)
rm(nbNodes) 
nbNodes <-  t(sapply(initBipartite, function(sbm){sbm$nbNodes}))
resEstimVBEM  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior,collecTau_init,estimOptions, emissionDist, model)

##################################"" TEST rParamZ
hyperparamPost <- resEstimVBEM$hyperparamPost
HSample <- rParamZ(MC, hyperparam = hyperparamPost, emissionDist, model ,nbNodes)
checkSample(HSample,MC,M,KRow,KCol,hyperparamPost,emissionDist,model)

hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau
HSample <- rParamZ(MC, hyperparam = hyperparamApproxPost , emissionDist, model ,nbNodes)
checkSample(HSample,MC,M,KRow,KCol,hyperparamApproxPost,emissionDist,model)



##################################

mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)
L <- likelihood(mydata,HSample,emissionDist  = 'bernoulli')
L <- likelihood(mydata,HSamplePrior,emissionDist  = 'bernoulli')

LPriorCP <- logDistConnectParam(HSample,MC,hyperparam,emissionDist)

LPriorBP <- logDirichletBlockProp(HSample,MC,hyperparam,model)
LPriorZ <- logMultinomZ(HSample, M, MC, hyperparamPrior,model)

LPrior <- logJointParamZ(HSample,M,MC,hyperparamPrior,emissionDist,model)
LPost <- logJointParamZ(HSample,M,MC,hyperparamApproxPost,emissionDist,model)

########################################"
mc <- 2
H.mc <- list()
H.mc$Z <- lapply(1:M, function(m){list(row = HSample$ZSample[[m]]$row[,,mc],col = HSample$ZSample[[m]]$col[,,mc])})
H.mc$connectParam <- HSample$connectParamSample[,,mc]
if (model == 'iidColBipartiteSBM'){
  H.mc$blockProp <- list(row = HSample$blockPropSample$row[,mc],col = HSample$blockPropSample$col[,mc])
}
if (model == 'piColBipartiteSBM'){
  H.mc$blockProp <- list(row = HSample$blockPropSample$row[,,mc],col = HSample$blockPropSample$col[,,mc])
}
H.mc.init <- H.mc


  
###################### Estim avec Kcol Krow true

paramsMCMC = list(nbIterMCMC = 10000)

H.mc.init$Z <-lapply(1:M,function(m){mySampler[[m]]$indMemberships})
H.mc$blockProp <- list(row = blockPropTrue$row[1,],col = blockPropTrue$col[1,])
H.mc$connectParam <- connectParamTrue

paramsMCMC$opEchan = list(connectParam = FALSE,blockProp = FALSE, Z = TRUE)



resMCMC <- MCMCKernel(data = mydata, H.mc.init, alpha.t = 1, hyperparApproxPost, hyperparamPrior, emissionDist = 'bernoulli', model = 'iidColBipartiteSBM',paramsMCMC, opSave=TRUE)

k  = sample(1:KRow,1)
l = sample(1:KCol,1)

par(mfrow=c(2,2))
burnin  = 5000 
extr <- burnin:paramsMCMC$nbIterMCMC
plot(density(resMCMC$seqConnectParam[k,l,extr]),main='alpha')  
curve(dbeta(x,hyperparamApproxPost$connectParam$alpha[k,l],hyperparamApproxPost$connectParam$beta[k,l]),col='red',add=TRUE)
abline(v = connectParamTrue$mean[k,l])
plot(density(resMCMC$seqBlockProp$row[k,extr]),main='pi_k')  
curve(dbeta(x,hyperparamApproxPost$blockProp$row[k], sum(hyperparamApproxPost$blockProp$row)-hyperparamApproxPost$blockProp$row[k]),col='red',add=TRUE)
abline(v = blockPropTrue$row[1,k])
plot(density(resMCMC$seqBlockProp$col[l,extr]),main='rho_l')  
curve(dbeta(x,hyperparamApproxPost$blockProp$col[l], sum(hyperparamApproxPost$blockProp$col)-hyperparamApproxPost$blockProp$col[l]),col='red',add=TRUE)
abline(v = blockPropTrue$col[1,l])


