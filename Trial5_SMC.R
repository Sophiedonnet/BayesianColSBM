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

#--------------------------------------- 
pathtoresSimu <- paste0(getwd(),'/Simu/Res_SIMU_2023806_22')


emissionDist = 'bernoulli'
model = 'piColSBM'

mySeed <- sample(1:999,1)
#mySeed <- 460
set.seed(mySeed)
###########################################################################################
############# simulation 
#############################################################################################
M = 1
K  = 4
 
#-----  block proportions simul iidColSBM 
blockPropTrue <- list()
if(model=='iidColSBM'){
  blockPropTrue  <-  matrix(rdirichlet(1,rep(K,K)),byrow = TRUE,nrow = M, ncol = K) #### not emptying some blocks in certain netwokrs 
}

#-----  block proportions simul piColSBM avec classes vides
if(model=='piColSBM'){
  blockPropTrue <-  rdirichlet(M,rep(1/(K-1),K))  #### emptying some blocks in certain netwokrs 
  blockPropTrue[1,] <- rep(1/K,K) 
}

#-----  connectivity matrix
connectParamTrue <- list(mean = round( matrix(rbeta(K*K,1/1.1,1/1.1), K, K)  ,9))
diag(connectParamTrue$mean) <- seq(0.9,0.1,len=K)
or <- order(diag(connectParamTrue$mean),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[or,or]
#-------  sizes of networks
nbNodes <- sample(10*c(6:10),M,replace = TRUE)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleSimpleSBM(nbNodes = nbNodes[m], blockProp =  blockPropTrue[m,],  connectParamTrue,directed = TRUE)})

propNoise = 0.1

collecNetworks <- lapply(mySampler,function(l){Net.m <- l$networkData; diag(Net.m) <- 0; Net.m})
noisyCollecNetworks <- noiseSampling(collecNetworks,propNoise)



mydata <- list(collecNetworks = noisyCollecNetworks, M= M, nbNodes = nbNodes)

names_file_data <- paste0(pathtoresSimu,'/mySimu.Rdata')
save(mydata,file = names_file_data )



###########################################################################################
#------------------ VBEM 
###########################################################################################
#--------------- init CollecTau
initSimple <- lapply(noisyCollecNetworks,estimateSimpleSBM)
myRef = which.max(sum(t(sapply(initSimple,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initSimple, ref = myRef)
KEstim <- initSimple[[myRef]]$nbBlocks[1]
print(c(K,KEstim))

#---------------  Prior distribution

hyperparamPrior <- setHyperparamPrior(M,K, emissionDist, model)


############################################### 
########### ORACLE 
#################################################
S <- vapply(1:M,function(m){t(mySampler[[m]]$indMemberships) %*% collecNetworks[[m]] %*% mySampler[[m]]$indMemberships},matrix(0,K,K))
S <- apply(S,c(1,2),sum)
N <- vapply(1:M,function(m){t(mySampler[[m]]$indMemberships) %*% matrix(1,nbNodes[m],nbNodes[m]) %*% mySampler[[m]]$indMemberships},matrix(0,K,K))
N <- apply(N,c(1,2),sum)
connectParamOracle <- S/N
(1-propNoise)*connectParamTrue$mean  + propNoise*(1-connectParamTrue$mean)
############################################################
######################## VB ###############################
#############################################################
#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                         maxIterVE = 100,
                         valStopCritVE = 10^-10,
                         valStopCritVB = 10^-10,
                         epsTau = 10^-10)
nbNodes <-  t(sapply(initSimple, function(sbm){sbm$nbNodes}))
resEstimVBEM  <- VBEMColSBM(mydata,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)

#------------------ Set ApproxPost

hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau

names_file_VB <- paste0(pathtoresSimu,'/myResVBEM.Rdata')
save(hyperparamApproxPost, hyperparamPrior, file = names_file_VB )


###########################################################################################
#------------------  SMC 
###########################################################################################
estimOptionsSMC = list()
estimOptionsSMC$paramsMCMC <- list(nbIterMCMC=5)  
estimOptionsSMC$MC <- 1000
estimOptionsSMC$ESS.rate <- 0.9
estimOptionsSMC$cESS.rate <- 0.9
estimOptionsSMC$opSave <- FALSE
estimOptionsSMC$op.parallel  <- list(os =  .Platform$OS.type, mc.cores = 4); 
estimOptionsSMC$op.print<- TRUE
estimOptionsSMC$NB.iter.max  <- Inf # Inf
estimOptionsSMC$op.SMC.classic <- FALSE
resSMC <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC)

names_file_SMC <- paste0(pathtoresSimu,'/myResSMC.Rdata')
save(resSMC,estimOptionsSMC , file = names_file_SMC )

estimOptionsSMC$op.SMC.classic <- TRUE
rm(hyperparamApproxPost)
resSMC_classic <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)
names_file_SMC_classic <- paste0(pathtoresSimu,'/myResSMC_classic.Rdata')
save(resSMC,estimOptionsSMC , file = names_file_SMC_classic )

 
###################################################
################ PLOT Post 
#######################################################"
load(file =  paste0(pathtoresSimu,'/myResVBEM.Rdata'))

############## rho_h
par(mfrow=c(1,1))
plot(resSMC_classic$alpha.vec,type='l')
lines(resSMC$alpha.vec,col='red')

#############  # plot sum alpha_{kl}
par(mfrow=c(1,1))
plot(density(apply(resSMC_classic$HSample_end$connectParamSample,c(3),sum),weights=resSMC_classic$W.end),main='alpha',col='green')
lines(density(apply(resSMC$HSample_end$connectParamSample,c(3),sum),weights=resSMC$W.end),col='magenta')
simu_ApproxPost <- rParamZ(10000,hyperparam = hyperparamApproxPost,emissionDist = 'bernoulli',model = model,nbNodes = nbNodes)
lines(density(apply(simu_ApproxPost$connectParamSample,c(3),sum)),col='red')

##################################################### 

par(mfrow=c(1,1))
for (k in 1:K){
  if (k ==1){plot(density(resSMC_classic$HSample_end$connectParamSample[k,k,],weights=resSMC_classic$W.end),main='alphakk',col='green')
  }else{
    lines(density(resSMC_classic$HSample_end$connectParamSample[k,k,],weights=resSMC_classic$W.end),col='green')
  }  
  curve(dbeta(x,hyperparamApproxPost$blockProp[k], sum(hyperparamApproxPost$blockProp)-hyperparamApproxPost$blockProp[k]),col='red',add=TRUE)
  lines(density(resSMC$HSample_end$connectParamSample[k,k,],weights=resSMC$W.end),col='magenta')
}


plot(density(apply(resSMC_classic$HSample_end$connectParamSample,c(3),sum),weights=resSMC_classic$W.end),main='alpha',col='green')
lines(density(apply(resSMC$HSample_end$connectParamSample,c(3),sum),weights=resSMC$W.end),col='magenta')
simu_ApproxPost <- rParamZ(10000,hyperparam = hyperparamApproxPost,emissionDist = 'bernoulli',model = model,nbNodes = nbNodes)
lines(density(apply(simu_ApproxPost$connectParamSample,c(3),sum)),col='red')


####################################################### 

###################

par(mfrow=c(1,1))
for (k in 1:K){
  if (k ==1){plot(density(resSMC_classic$HSample_end$blockPropSample[1,k,],weights=resSMC_classic$W.end),col='green',main='pi_k',xlim = c(0,1))
    }else{
    lines(density(resSMC_classic$HSample_end$blockPropSample[1,k,],weights=resSMC_classic$W.end),col='green')
    }  
    #curve(dbeta(x,hyperparamApproxPost$blockProp[k], sum(hyperparamApproxPost$blockProp)-hyperparamApproxPost$blockProp[k]),col='red',add=TRUE)
    lines(density(resSMC$HSample_end$blockPropSample[1,k,],weights=resSMC$W.end),col='magenta')
}
  
################################### 
#####################################
##################  Proba d'être dans le même clusterpar VB
############################################

ZZ_VB <- hyperparamApproxPost$collecTau[[1]]%*%t(hyperparamApproxPost$collecTau[[1]])

ZSample <- lapply(1:M,function(m){resSMC$HSample_end$ZSample[[m]]})
postZSample <- transfZsampleIntoMatrix(ZSample)
ZZ_SMC  = matrix(0,nbNodes[1],nbNodes[1])
for (i in 1:nbNodes[1]){
  for (j in 1:nbNodes[1]){
    ZZ_SMC [i,j] = sum(resSMC$W.end*(postZSample[[1]][i,]== postZSample[[1]][j,])) 
  }
}

ZSample <- lapply(1:M,function(m){resSMC_classic$HSample_end$ZSample[[m]][,,extr]})
postZSample <- transfZsampleIntoMatrix(ZSample)
ZZ_SMC_classic  = matrix(0,nbNodes[1],nbNodes[1])
for (i in 1:nbNodes[1]){
  for (j in 1:nbNodes[1]){
    ZZ_SMC_classic [i,j] = sum(resSMC_classic$W.end*(postZSample[[1]][i,]== postZSample[[1]][j,])) 
  }
}


plotMyMatrix(ZZ_VB)
plotMyMatrix(ZZ_SMC-ZZ_VB,plotOptions = list(legend = TRUE))
plotMyMatrix(ZZ_SMC-ZZ_SMC_classic,plotOptions = list(legend = TRUE))
plot(ZZ_SMC,ZZ_SMC_classic); abline(a=0,b = 1)
plot(ZZ_VB,ZZ_SMC);abline(a=0,b = 1)
hist(ZZ_VB)
hist(ZZ_SMC)
hist(ZZ_SMC_classic)

