rm(list=ls())
library(sbm)
library(gtools)

source('Functions/utils.R')
source('Functions/VBEMBipartiteColSBM.R')
source('Functions/initializationCollecTau.R')
source('Functions/Functions-SMC.R')
source('Functions/Functions-colBipartiteSBM.R')
source('Functions/Functions-checks.R')


emissionDist = 'bernoulli'
model = 'piColBipartiteSBM'
model = 'iidColBipartiteSBM'
###########################################################################################
############# simulation avec iid colSBM mais certains très petits réseaux. 
#############################################################################################


M = 6
KRow = 4
KCol = 3
#---  block proportions simul iidColSBM 
blockProp <- list()
blockProp$row <-  matrix(rdirichlet(1,rep(KRow,KRow)),byrow = TRUE,nrow = M, ncol = KRow) #### emptying some blocks in certain netwokrs 
blockProp$col <- matrix(rdirichlet(1,rep(KCol,KCol)),byrow = TRUE,nrow = M, ncol = KCol)   #### emptying some blocks in certain netwokros
print(blockProp)

#-----  connectivity matrix
means <- matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  
connectParam <- list(mean = round(means,2))
oCol <- order(colSums(connectParam$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParam$mean[,]),decreasing = TRUE)
connectParam$mean <- connectParam$mean[oRow,oCol]


#-------  sizes of networks
nbNodes <- matrix(sample(10*c(6:10),M*2,replace = TRUE),M,2)
nbNodes[5,] = c(30,20)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockProp$row[m,],col=blockProp$col[m,]),  connectParam)})
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

hyperparamPost <- resEstimVBEM$hyperparamPost
hyperparamPost$collecTau <- resEstimVBEM$collecTau
HSample <- rParamZ(MC, hyperparam = hyperparamPost , emissionDist, model ,nbNodes)
checkSample(HSample,MC,M,KRow,KCol,hyperparamPost,emissionDist,model)



##################################

data <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)
L <- likelihood(data,HSample,emissionDist  = 'bernoulli')
L <- likelihood(data,HSamplePrior,emissionDist  = 'bernoulli')

LPriorCP <- logDist.connectParam(HSample,MC,hyperparam,emissionDist)

LPriorBP <- logDirichlet.blockProp(HSample,MC,hyperparam,model)
LPriorZ <- logMultinom.Z(HSample, M, MC, hyperparamPrior,model)

LPrior <- logJointParamZ(HSample,M,MC,hyperparamPrior,emissionDist,model)
LPost <- logJointParamZ(HSample,M,MC,hyperparamPost,emissionDist,model)

 
###################### Estim avec Kcol Krow true

mc <- sample(1:MC,1)
H.mc = list()
H.mc$connectParam <- HSample$connectParamSample[,,mc]
if(model == 'piColBipartiteSBM'){
  H.mc$blockProp <- list(row = HSample$blockPropSample$row[,,mc], col = HSample$blockPropSample$col[,,mc])}
if(model == 'iidColBipartiteSBM'){
  H.mc$blockProp <- list(row = HSample$blockPropSample$row[,mc], col = HSample$blockPropSample$col[,mc])}

H.mc$Z<- lapply(1:M,function(m){list(row = HSample$ZSample[[m]]$row[,,mc],col = HSample$ZSample[[m]]$col[,,mc])})


estimOptions <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)







hyperparamPrior_picol <- setHyperparamPrior(M,KRowEstim,KColEstim, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
resEstim_picol  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_picol,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
logLikMarg_picol <- computeLogLikMarg(collecNetworks, resEstim_picol$collecTau,hyperparamPrior_picol, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
