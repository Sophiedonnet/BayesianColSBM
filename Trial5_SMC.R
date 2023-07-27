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
pathtoresSimu <- paste0(getwd(),'/Simu/Res_SIMU_2023_06_26')


emissionDist = 'bernoulli'
model = 'piColSBM'
mySeed <- sample(1:999,1)
#mySeed <- 460
set.seed(mySeed)
###########################################################################################
############# simulation 
#############################################################################################
names_file_data <- paste0(pathtoresSimu,'/mySimu.Rdata')
simu = FALSE
if(simu){
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
  
  save(mydata,file = names_file_data )
}else{
  load(names_file_data)  
}


###########################################################################################
#------------------ VBEM 
###########################################################################################

names_file_VB <- paste0(pathtoresSimu,'/myResVBEM.Rdata')
estim = TRUE
#--------------- init CollecTau
if(estim){
  initSimple <- lapply(mydata$collecNetworks,estimateSimpleSBM)
  myRef = which.max(sum(t(sapply(initSimple,function(m){m$nbBlocks}))))
  collecTau_init <- initCollecTau(initSimple, ref = myRef)
  KEstim <- initSimple[[myRef]]$nbBlocks[1]
  print(c(KEstim))
  
  #---------------  Prior distribution
  
  hyperparamPrior <- setHyperparamPrior(mydata$M,KEstim, emissionDist, model)
  
  
  ############################################### 
  ########### ORACLE 
  #################################################
  # S <- vapply(1:M,function(m){t(mySampler[[m]]$indMemberships) %*% collecNetworks[[m]] %*% mySampler[[m]]$indMemberships},matrix(0,K,K))
  # S <- apply(S,c(1,2),sum)
  # N <- vapply(1:M,function(m){t(mySampler[[m]]$indMemberships) %*% matrix(1,nbNodes[m],nbNodes[m]) %*% mySampler[[m]]$indMemberships},matrix(0,K,K))
  # N <- apply(N,c(1,2),sum)
  # connectParamOracle <- S/N
  # (1-propNoise)*connectParamTrue$mean  + propNoise*(1-connectParamTrue$mean)
  # ############################################################
  ######################## VB ###############################
  #############################################################
  #------------- variational estim 
  
  estimOptionsVBEM <- list(maxIterVB = 3000,
                           maxIterVE = 100,
                           valStopCritVE = 10^-10,
                           valStopCritVB = 10^-10,
                           epsTau = 10^-10)
  nbNodes <-  t(sapply(initSimple, function(sbm){sbm$nbNodes}))
  resEstimVBEM  <- VBEMColSBM(mydata,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)
  
  #------------------ Set ApproxPost
  
  hyperparamApproxPost <- resEstimVBEM$hyperparamPost
  hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau
  save(hyperparamApproxPost, hyperparamPrior, file = names_file_VB )
}else{
  load(names_file_VB)
}



###########################################################################################
#------------------  SMC 
###########################################################################################
names_file_SMC_classic <- paste0(pathtoresSimu,'/myResSMC_classic.Rdata')
names_file_SMC <- paste0(pathtoresSimu,'/myResSMC.Rdata')

if(estim){
  estimOptionsSMC = list()
  estimOptionsSMC$paramsMCMC <- list(nbIterMCMC=5)  
  estimOptionsSMC$MC <- 2000
  estimOptionsSMC$ESS.rate <- 0.9
  estimOptionsSMC$cESS.rate <- 0.9
  estimOptionsSMC$opSave <- FALSE
  estimOptionsSMC$op.parallel  <- list(os =  .Platform$OS.type, mc.cores = 4); 
  estimOptionsSMC$op.print<- TRUE
  estimOptionsSMC$NB.iter.max  <- Inf # Inf
  estimOptionsSMC$op.SMC.classic <- FALSE
  resSMC <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC)
  
  save(resSMC,estimOptionsSMC , file = names_file_SMC )
  
  estimOptionsSMC$op.SMC.classic <- TRUE
  rm(hyperparamApproxPost)
  resSMC_classic <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)
  save(resSMC_classic,estimOptionsSMC , file = names_file_SMC_classic )
}else{
  load(names_file_SMC_classic)
  load(names_file_SMC)
}


###################################################
################ PLOT Post 
#######################################################"
load(names_file_VB)
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
K  = KEstim
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

ZSample <- lapply(1:mydata$M,function(m){resSMC$HSample_end$ZSample[[m]]})
postZSample <- transfZsampleIntoMatrix(ZSample)
ZZ_SMC  = matrix(0,mydata$nbNodes[1],mydata$nbNodes[1])
for (i in 1:mydata$nbNodes[1]){
  for (j in 1:mydata$nbNodes[1]){
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
plotMyMatrix(ZZ_SMC)
plotMyMatrix(ZZ_SMC-ZZ_VB,plotOptions = list(legend = TRUE))
plotMyMatrix(ZZ_SMC-ZZ_SMC_classic,plotOptions = list(legend = TRUE))

ZZ_DF <- cbind(c(ZZ_VB),c(ZZ_SMC),c(ZZ_SMC_classic))
ZZ_DF <- as.data.frame(ZZ_DF)
names(ZZ_DF) <- c('VB', 'SMC_VB','SMC_Classic')
library(ggplot2)
whereToSave <- "/home/sophie/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2023/2023_06_MC/plots/"
ggplot(ZZ_DF,aes(x=VB,y=SMC_VB)) + geom_point(size=4) + geom_abline(intercept=0, slope = 1) #+ ggtitle("Co-clustering probability") 
ggsave(paste0(whereToSave,'/Cocluster2.png'), width = 15, height = 15, units = "cm")

ggplot(ZZ_DF,aes(x=SMC_Classic,y=SMC_VB)) + geom_point(size=4) + geom_abline(intercept=0, slope = 1)# +  ggtitle("Co-clustering probability") 
ggsave(paste0(whereToSave,'/Cocluster1.png'), width = 15, height = 15, units = "cm")


plot(ZZ_SMC,ZZ_SMC_classic); abline(a=0,b = 1)
plot(ZZ_VB,ZZ_SMC_classic);abline(a=0,b = 1)
hist(ZZ_VB)
hist(ZZ_SMC)
hist(ZZ_SMC_classic)

hist(entropyBernoulli(ZZ_VB))
hist(entropyBernoulli(ZZ_SMC))

