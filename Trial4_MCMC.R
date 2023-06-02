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
#model = 'iidColBipartiteSBM'

###########################################################################################
############# simulation avec iid colSBM mais certains très petits réseaux. 
#############################################################################################


M = 4
K  = 4
 
#-----  block proportions simul iidColSBM 
#blockPropTrue <- list()
#blockPropTrue <-  matrix(rdirichlet(1,rep(K,K)),byrow = TRUE,nrow = M, ncol = K) #### not emptying some blocks in certain netwokrs 
#blockPropTrue$col <- matrix(rdirichlet(1,rep(KCol,KCol)),byrow = TRUE,nrow = M, ncol = KCol)   #### not  emptying some blocks in certain netwokros

#-----  block proportions simul piColSBM avec classes vides
blockPropTrue <-  rdirichlet(M,rep(1/(K-1),K))  #### emptying some blocks in certain netwokrs 
blockPropTrue[1,] <- rep(1/K,K) 
print(blockPropTrue)

#-----  connectivity matrix
connectParamTrue <- list(mean = round(matrix(rbeta(K*K,1/1,1/1), K, K)  ,2))
or <- order(diag(connectParamTrue$mean),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[or,or]
#-------  sizes of networks
nbNodes <- sample(10*c(5:10),M,replace = TRUE)
#nbNodes[1,] = c(150,200)
#nbNodes[5,] = c(30,20)

#---------  SIMULATION

mySampler <- lapply(1:M, function(m){sampleSimpleSBM(nbNodes = nbNodes[m], blockProp =  blockPropTrue[m,],  connectParamTrue,directed = TRUE)})
collecNetworks <- lapply(mySampler,function(l){Net.m <- l$networkData; diag(Net.m) <- 0; Net.m})

#--------------- init CollecTau
initSimple <- lapply(collecNetworks,estimateSimpleSBM)
myRef = which.max(sum(t(sapply(initSimple,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initSimple, ref = myRef)
KEstim <- initSimple[[myRef]]$nbBlocks[1]

print(c(KEstim))

#---------------  Prior distribution

hyperparamPrior <- setHyperparamPrior(M,K, emissionDist, model)

#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-10,
                     valStopCritVB = 10^-10)
nbNodes <-  t(sapply(initSimple, function(sbm){sbm$nbNodes}))
resEstimVBEM  <- VBEMColSBM(collecNetworks,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)

hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau




##################################
mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)


HSample <- rParamZ(200, hyperparam = hyperparamPrior , emissionDist, model ,nbNodes)
mc <- 1
H.mc.sim <- list(connectParam = HSample$connectParamSample[,,mc])
H.mc.sim$Z <- lapply(1:M, function(m){HSample$ZSample[[m]][,,mc]})
if (model == 'iidColSBM'){
  H.mc.sim$blockProp <- HSample$blockPropSample[,mc]
}
if (model == 'piColSBM'){
  H.mc.sim$blockProp <-  HSample$blockPropSample[,,mc]
}

  
###################### Estim avec Kcol K true

paramsMCMC = list(nbIterMCMC = 10000)
H.mc.init <- H.mc.sim
paramsMCMC$opEchan = list(connectParam = TRUE, blockProp = TRUE, Z = TRUE)
if (!paramsMCMC$opEchan$connectParam){H.mc.init$connectParam <- connectParamTrue$mean}
if (!paramsMCMC$opEchan$blockProp){H.mc.init$blockProp <- blockPropTrue}
if(!paramsMCMC$opEchan$Z){for (m in 1:M){H.mc.init$Z[[m]] <-mySampler[[m]]$indMemberships}}


resMCMC <- MCMCKernel(mydata, H.mc.init, alpha.t = 1, hyperparamPrior,hyperparamApproxPost = NULL, emissionDist = 'bernoulli', model =  'piColSBM',paramsMCMC, opSave= TRUE)

burnin  = 2500 
extr <- burnin:paramsMCMC$nbIterMCMC

m = sample(1:M,1)
apply(resMCMC$seqZ[[m]][,,extr],c(1,2),mean)[1:10,]
resMCMC$seqZ[[m]][,,1][1:10,]
hyperparamApproxPost$collecTau[[m]][1:10,]


apply(resMCMC$seqZ[[m]]$col[,,extr],c(1,2),mean)[1:10,]
resMCMC$seqZ[[m]]$col[,,1][1:10,]
hyperparamApproxPost$collecTau[[m]]$col[1:10,]

par(mfrow=c(K,KCol))
for (k in 1:K){
  for (l in 1:KCol){
 #   if ((k ==1) & (l==1)){
      plot(density(resMCMC$seqConnectParam[k,l,extr]),main='alpha',xlim=c(0,1),col='green')
#    }else{
#      lines(density(resMCMC$seqConnectParam[k,l,extr]))  
#      }
    curve(dbeta(x,hyperparamApproxPost$connectParam$alpha[k,l],hyperparamApproxPost$connectParam$beta[k,l]),col='red',add=TRUE)
    abline(v = connectParamTrue$mean)
  }
}

#------------------------------------ 
if(model=="iidColBipartiteSBM"){
  par(mfrow=c(1,1))
  for (k in 1:K){
    if (k ==1){plot(density(resMCMC$seqBlockProp[k,extr]),main='pi_k',xlim = c(0,1),col='green')
      }else{
        lines(density(resMCMC$seqBlockProp[k,extr]),col='green')
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp[k], sum(hyperparamApproxPost$blockProp)-hyperparamApproxPost$blockProp[k]),col='red',add=TRUE)
    abline(v = blockPropTrue[1,k])
  }
  
  for (l in 1:KCol){
    if (l ==1){plot(density(resMCMC$seqBlockProp$col[l,extr]),main='rho_k',xlim = c(0,1),col='green')
    }else{
      lines(density(resMCMC$seqBlockProp$col[l,extr]),col='green')
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp$col[l], sum(hyperparamApproxPost$blockProp$col)-hyperparamApproxPost$blockProp$col[l]),col='red',add=TRUE)
    abline(v = blockPropTrue$col[1,l])
  }
  
}

if(model=="piColBipartiteSBM"){
  par(mfrow=c(floor(M/2)+1,2))
  for (m in 1:M){
    for (k in 1:K){
      if (k ==1){plot(density(resMCMC$seqBlockProp[m,k,extr]),main='pi_k',xlim = c(0,1),col='green')
    }else{
      lines(density(resMCMC$seqBlockProp[m,k,extr]),col='green')
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp[m,k], sum(hyperparamApproxPost$blockProp[m,])-hyperparamApproxPost$blockProp[m,k]),col='red',add=TRUE)
    abline(v = blockPropTrue[m,k])
    }
  }  
  par(mfrow=c(floor(M/2)+1,2))
  for (m in 1:M){
    for (k in 1:KCol){
      if (k ==1){plot(density(resMCMC$seqBlockProp$col[m,k,extr]),main='pi_k',xlim = c(0,1),col='green')
      }else{
        lines(density(resMCMC$seqBlockProp$col[m,k,extr]),col='green')
      }  
      curve(dbeta(x,hyperparamApproxPost$blockProp$col[m,k], sum(hyperparamApproxPost$blockProp$col[m,])-hyperparamApproxPost$blockProp$col[m,k]),col='red',add=TRUE)
      abline(v = blockPropTrue$col[m,k])
    }
  }
}




