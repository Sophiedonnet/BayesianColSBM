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
or <- order(diag(connectParamTrue$mean),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[or,or]
#-------  sizes of networks
nbNodes <- sample(10*c(6:10),M,replace = TRUE)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleSimpleSBM(nbNodes = nbNodes[m], blockProp =  blockPropTrue[m,],  connectParamTrue,directed = TRUE)})



collecNetworks <- lapply(mySampler,function(l){Net.m <- l$networkData; diag(Net.m) <- 0; Net.m})
noisyCollecNetworks <- noiseSampling(collecNetworks,propNoise = 0.2)



mydata <- list(collecNetworks = noisyCollecNetworks, M= M, nbNodes = nbNodes)


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

#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                         maxIterVE = 100,
                         valStopCritVE = 10^-10,
                         valStopCritVB = 10^-10,
                         epsTau = 10^-3)
nbNodes <-  t(sapply(initSimple, function(sbm){sbm$nbNodes}))
resEstimVBEM  <- VBEMColSBM(mydata,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)

#------------------ Set ApproxPost

hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau




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

#estimOptionsSMC$op.SMC.classic <- TRUE
#resSMC_classic <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)
#save(mydata,resSMC,hyperparamApproxPost,file='myTrialSMCResults.Rdata')

estimOptionsSMC$op.SMC.classic <- FALSE
resSMC <- resSCM <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC)

plot(resSMC$alpha.vec,type='l')




###################### Estim avec  K true

par(mfrow=c(K,K))
for (k in 1:K){
  for (l in 1:K){
 #   if ((k ==1) & (l==1)){
      plot(density(resSMC$HSample_end$connectParamSample[k,l,],weights=resSMC$W.end),main='alpha',xlim=c(0,1),col='green')
#    }else{
#      lines(density(resMCMC$seqConnectParam[k,l,extr]))  
#      }
    curve(dbeta(x,hyperparamApproxPost$connectParam$alpha[k,l],hyperparamApproxPost$connectParam$beta[k,l]),col='red',add=TRUE)
    abline(v = connectParamTrue$mean[k,l])
  }
}

#------------------------------------ 
if(model=="iidColSBM"){
  par(mfrow=c(1,1))
  for (k in 1:K){
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

if(model=="piColSBM"){
  par(mfrow=c(floor(M/2)+1,2))
  for (m in 1:M){
    for (k in 1:K){
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




