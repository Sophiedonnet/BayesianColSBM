rm(list=ls())
#------------------ set wd
if(Sys.info()[[4]]=='sophie-Latitude-5310'){
  setwd('/home/sophie/WORK_LOCAL/RECHERCHE/TRAVAUX_DE_RECHERCHE/Barbillon-Chabert/BayesianBipartiteColSBM')
}
if(Sys.info()[[4]]=="donnet-Precision-Tower-5810"){
  setwd('/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianBipartiteColSBM')
}

#--------------- packages and functions
library(sbm)
library(gtools)
library(DescTools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)

###############################################################
#######################        " TEST MI  functions 
###############################################################


MC = 1000
n  = 10
matZ  = matrix(0,n,MC)
for(i in 1:n){
  K = sample(1:5,1)
  matZ[i,] <- t(rmultinom(MC,1,rdirichlet(1,rep(1,K))))%*%matrix(1:K,ncol=1) 
}
mutualInformation(matZ)  
MutInf(matZ[1,], matZ[2,], base = 2)


###############################################################
########  TEST onf my models MODELS 
###############################################################

emissionDist = 'bernoulli'
model = 'iidColBipartiteSBM'


###########################################################################################
############# simulation 
#############################################################################################
M = 3
KRow = 4
KCol = 3


#-----  block proportions simul piColSBM avec classes vides
mySeed = sample(1:10000,1)
whereSaveData <- paste0(getwd(),'/Simu/Data/',model,'/data_simu_',mySeed,'.Rdata')
print(mySeed)
set.seed(mySeed)
blockPropTrue <- list()
blockPropTrue$row <-  rdirichlet(M,rep(1/1.2,KRow))  #### emptying some blocks in certain netwokrs
blockPropTrue$col <- rdirichlet(M,rep(1/1.2,KCol))   #### emptying some blocks in certain netwokros
blockPropTrue$row[1,] <- rep(1/KRow,KRow)
blockPropTrue$col[1,] <- rep(1/KCol,KCol)

print(blockPropTrue)

#-----  connectivity matrix
connectParamTrue <- list(mean = round( matrix(rbeta(KRow*KCol,1/1.01,1/1.01), KRow, KCol)  ,9))
oCol <- order(colSums(connectParamTrue$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParamTrue$mean[,]),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[oRow,oCol]

if(model=="iidColBipartiteSBM"){
  blockPropTrue$row <- blockPropTrue$row[1,]
  blockPropTrue$col <- blockPropTrue$col[1,]
}
  
  

#-------  sizes of networks
nbNodes <- matrix(sample(15:35,M*2,replace = TRUE),M,2)
#---------  SIMULATION
if (model =='iidColBipartiteSBM'){
  mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row,col=blockPropTrue$col),  connectParamTrue)})
}else{
  mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row[m,],col=blockPropTrue$col[m,]),  connectParamTrue)})
}
collecNetworks <- lapply(mySampler,function(l){l$networkData})
mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)


##############################################################
#---------- VEM
##############################################################
 


#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef <- which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]
print(c(KRow,KCol))


###################################################### 
#---------------  Prior distribution
######################################################
hyperparamPrior <- setHyperparamPrior(M,KRow,KCol, emissionDist, model)
HSample <- rParamZ(MC=1000, hyperparamPrior, emissionDist, model,nbNodes);
ZSample <-  HSample$ZSample
MIPrior <- myMutualInformationZ(ZSample,1:5)
MIPrior

###########################################################################################
#------------------ VBEM  + SMC -VBEM 
###########################################################################################


#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                         maxIterVE = 100,
                         valStopCritVE = 10^-5,
                         valStopCritVB = 10^-5)

resEstimVBEM  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)
LogLikMarg_VB <- sum(computeLogLikMarg_VB(collecNetworks,resEstimVBEM$collecTau,hyperparamPrior,emissionDist,model))



whereSaveResEstimVB <- paste0(getwd(),'/Simu/Res/',model,'/VB/resVB_',mySeed,'.Rdata')
#save(initBipartite,estimOptionsVBEM, resEstimVBEM_piCol,estimOptionsVBEM,hyperparamPrior_piCol,file=whereSaveResEstimVB)

#------------------ Set ApproxPost
hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau

HSample <- rParamZ(MC=10000, hyperparamApproxPost, emissionDist, model,nbNodes);
ZSample <-  HSample$ZSample
MIApproxPost <- myMutualInformationZ(HSample$ZSample,1:5)
MIApproxPost

estimOptionsSMC = list()
estimOptionsSMC$paramsMCMC <- list(nbIterMCMC=10)  
estimOptionsSMC$MC <- 1000
estimOptionsSMC$ESS.rate <- 0.9
estimOptionsSMC$cESS.rate <- 0.9
estimOptionsSMC$opSave <- FALSE
estimOptionsSMC$op.parallel  <- list(os =  .Platform$OS.type, mc.cores = 1); 
estimOptionsSMC$op.print<- TRUE
estimOptionsSMC$NB.iter.max  <- Inf # Inf
estimOptionsSMC$op.SMC.classic <- FALSE

resSMC_VB <- SMCColBipartiteSBM(data = mydata,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC)

List.matZ.0 <- transfZsampleIntoMatrix(resSMC_VB$HSample.0$ZSample)
List.matZ <- transfZsampleIntoMatrix(resSMC_VB$HSample_end$ZSample)


