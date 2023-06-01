rm(list=ls())
library(sbm)
library(gtools)

source('Functions/utils.R')
source('Functions/VBEMBipartiteColSBM.R')
source('Functions/initializationCollecTau.R')
source('Functions/Functions-SMC-colBipartiteSBM.R')
source('Functions/Functions-colBipartiteSBM.R')
source('Functions/Functions-checks.R')
source('Functions/Function-MCMC-colBipartiteSBM.R')

emissionDist = 'bernoulli'
model = 'piColBipartiteSBM'
model = 'iidColBipartiteSBM'

mySeed <- sample(1:999,1)
#mySeed <- 460
set.seed(mySeed)
###########################################################################################
############# simulation 
#############################################################################################
M = 5
KRow = 4
KCol = 3
#-----  block proportions simul iidColSBM 
blockPropTrue <- list()
if(model=='iidColBipartiteSBM'){
  blockPropTrue$row <-  matrix(rdirichlet(1,rep(KRow,KRow)),byrow = TRUE,nrow = M, ncol = KRow) #### not emptying some blocks in certain netwokrs 
  blockPropTrue$col <- matrix(rdirichlet(1,rep(KCol,KCol)),byrow = TRUE,nrow = M, ncol = KCol)   #### not  emptying some blocks in certain netwokros
}

#-----  block proportions simul piColSBM avec classes vides
if(model=='piColBipartiteSBM'){
  blockPropTrue$row <-  rdirichlet(M,rep(1/(KRow-1),KRow))  #### emptying some blocks in certain netwokrs 
  blockPropTrue$col <- rdirichlet(M,rep(1/(KCol-1),KCol))   #### emptying some blocks in certain netwokros
  blockPropTrue$row[1,] <- rep(1/KRow,KRow) 
  blockPropTrue$col[1,] <- rep(1/KCol,KCol)
}

#-----  connectivity matrix
connectParamTrue <- list(mean = round( matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  ,9))
oCol <- order(colSums(connectParamTrue$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParamTrue$mean[,]),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[oRow,oCol]
#-------  sizes of networks
nbNodes <- matrix(sample(10*c(6:10),M*2,replace = TRUE),M,2)
#nbNodes[1,] = c(150,200)
#nbNodes[5,] = c(30,20)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row[m,],col=blockPropTrue$col[m,]),  connectParamTrue)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})

save(collecNetworks,file='myTrialData.Rdata')
mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)

###########################################################################################
#------------------ VBEM 
###########################################################################################

#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef <- which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]

print(c(KRow,KCol))

#---------------  Prior distribution

hyperparamPrior <- setHyperparamPrior(M,KRow,KCol, emissionDist, model)

#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)

resEstimVBEM  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)

#------------------ Set ApproxPost
hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau



###########################################################################################
#------------------  SMC 
###########################################################################################
estimOptionsSMC = list()
estimOptionsSMC$paramsMCMC <- list(nbIterMCMC=2)  
estimOptionsSMC$MC <- 1000
estimOptionsSMC$ESS.rate <- 0.9
estimOptionsSMC$cESS.rate <- 0.9
estimOptionsSMC$opSave <- TRUE
estimOptionsSMC$op.parallel  <- list(os =  .Platform$OS.type, mc.cores = 1); 
estimOptionsSMC$op.print<- TRUE
estimOptionsSMC$NB.iter.max  <- Inf # Inf
estimOptionsSMC$op.SMC.classic <- FALSE
estimOptionsSMC$op.SMC.classic <- TRUE
resSMC <- SMCColBipartiteSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)


save(mydata,resSMC,hyperparamApproxPost,file='myTrialSMCResults.Rdata')

estimOptionsSMC$op.SMC.classic <- TRUE
resSCM_classic <- SMCColBipartiteSBM(data = mydata,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC)

plot(resSMC$alpha.vec,type='l')




###################### Estim avec Kcol Krow true

par(mfrow=c(KRow,KCol))
for (k in 1:KRow){
  for (l in 1:KCol){
 #   if ((k ==1) & (l==1)){
      plot(density(resMCMC$seqConnectParam[k,l,extr]),main='alpha',xlim=c(0,1))
#    }else{
#      lines(density(resMCMC$seqConnectParam[k,l,extr]))  
#      }
    curve(dbeta(x,hyperparamApproxPost$connectParam$alpha[k,l],hyperparamApproxPost$connectParam$beta[k,l]),col='red',add=TRUE)
    abline(v = connectParamTrue$mean[k,l])
  }
}

#------------------------------------ 
if(model=="iidColBipartiteSBM"){
  par(mfrow=c(1,1))
  for (k in 1:KRow){
    if (k ==1){plot(density(resMCMC$seqBlockProp$row[k,extr]),main='pi_k',xlim = c(0,1))
      }else{
        lines(density(resMCMC$seqBlockProp$row[k,extr]))
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp$row[k], sum(hyperparamApproxPost$blockProp$row)-hyperparamApproxPost$blockProp$row[k]),col='red',add=TRUE)
    abline(v = blockPropTrue$row[1,k])
  }
  
  for (l in 1:KCol){
    if (l ==1){plot(density(resMCMC$seqBlockProp$col[l,extr]),main='rho_k',xlim = c(0,1))
    }else{
      lines(density(resMCMC$seqBlockProp$col[l,extr]))
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp$col[l], sum(hyperparamApproxPost$blockProp$col)-hyperparamApproxPost$blockProp$col[l]),col='red',add=TRUE)
    abline(v = blockPropTrue$col[1,l])
  }
  
}

if(model=="piColBipartiteSBM"){
  par(mfrow=c(floor(M/2)+1,2))
  for (m in 1:M){
    for (k in 1:KRow){
      if (k ==1){plot(density(resMCMC$seqBlockProp$row[m,k,extr]),main='pi_k',xlim = c(0,1))
    }else{
      lines(density(resMCMC$seqBlockProp$row[m,k,extr]))
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp$row[m,k], sum(hyperparamApproxPost$blockProp$row[m,])-hyperparamApproxPost$blockProp$row[m,k]),col='red',add=TRUE)
    abline(v = blockPropTrue$row[m,k])
    }
  }  
  par(mfrow=c(floor(M/2)+1,2))
  for (m in 1:M){
    for (k in 1:KCol){
      if (k ==1){plot(density(resMCMC$seqBlockProp$col[m,k,extr]),main='pi_k',xlim = c(0,1))
      }else{
        lines(density(resMCMC$seqBlockProp$col[m,k,extr]))
      }  
      curve(dbeta(x,hyperparamApproxPost$blockProp$col[m,k], sum(hyperparamApproxPost$blockProp$col[m,])-hyperparamApproxPost$blockProp$col[m,k]),col='red',add=TRUE)
      abline(v = blockPropTrue$col[m,k])
    }
  }
}




