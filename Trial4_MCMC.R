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




emissionDist = 'bernoulli'
model = 'piColBipartiteSBM'
#model = 'iidColBipartiteSBM'

###########################################################################################
############# simulation avec iid colSBM mais certains très petits réseaux. 
#############################################################################################


M = 4
KRow = 4
KCol = 3
#-----  block proportions simul iidColSBM 
#blockPropTrue <- list()
#blockPropTrue$row <-  matrix(rdirichlet(1,rep(KRow,KRow)),byrow = TRUE,nrow = M, ncol = KRow) #### not emptying some blocks in certain netwokrs 
#blockPropTrue$col <- matrix(rdirichlet(1,rep(KCol,KCol)),byrow = TRUE,nrow = M, ncol = KCol)   #### not  emptying some blocks in certain netwokros

#-----  block proportions simul piColSBM avec classes vides
blockPropTrue <- list()
blockPropTrue$row <-  rdirichlet(M,rep(1/(KRow-1),KRow))  #### emptying some blocks in certain netwokrs 
blockPropTrue$col <- rdirichlet(M,rep(1/(KCol-1),KCol))   #### emptying some blocks in certain netwokros
blockPropTrue$row[1,] <- rep(1/KRow,KRow) 
blockPropTrue$col[1,] <- rep(1/KCol,KCol)
print(blockPropTrue)

#-----  connectivity matrix
connectParamTrue <- list(mean = round(matrix(rbeta(KRow*KCol,1/1,1/1), KRow, KCol)  ,2))
oCol <- order(colSums(connectParamTrue$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParamTrue$mean[,]),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[oRow,oCol]
#-------  sizes of networks
nbNodes <- matrix(sample(10*c(5:10),M*2,replace = TRUE),M,2)
#nbNodes[1,] = c(150,200)
#nbNodes[5,] = c(30,20)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row[m,],col=blockPropTrue$col[m,]),  connectParamTrue)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})


#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)

for(m in 1:M){

  b.m <- which( (initBipartite[[m]]$storedModels[,3]==KRow) & (initBipartite[[m]]$storedModels[,4]==KCol))
  print(c(m,b.m))
  if (length(b.m)>0){initBipartite[[m]]$setModel(b.m)}
}
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
nbNodes <-  t(sapply(initBipartite, function(sbm){sbm$nbNodes}))
resEstimVBEM  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)

hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau




##################################
mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)


HSample <- rParamZ(2, hyperparam = hyperparamPrior , emissionDist, model ,nbNodes)
mc <- 1
H.mc.sim <- list(connectParam = HSample$connectParamSample[,,mc])
H.mc.sim$Z <- lapply(1:M, function(m){list(row = HSample$ZSample[[m]]$row[,,mc],col = HSample$ZSample[[m]]$col[,,mc])})
if (model == 'iidColBipartiteSBM'){
  H.mc.sim$blockProp <- list(row = HSample$blockPropSample$row[,mc],col = HSample$blockPropSample$col[,mc])
}
if (model == 'piColBipartiteSBM'){
  H.mc.sim$blockProp <- list(row = HSample$blockPropSample$row[,,mc],col = HSample$blockPropSample$col[,,mc])
}

  
###################### Estim avec Kcol Krow true

paramsMCMC = list(nbIterMCMC = 10000)
H.mc.init <- H.mc.sim
paramsMCMC$opEchan = list(connectParam = TRUE, blockProp = TRUE, ZRow = TRUE, ZCol = TRUE)
if (!paramsMCMC$opEchan$connectParam){H.mc.init$connectParam <- connectParamTrue$mean}
if (!paramsMCMC$opEchan$blockProp){H.mc.init$blockProp <- list(row = blockPropTrue$row[1,],col = blockPropTrue$col[1,])}
if(!paramsMCMC$opEchan$ZRow){for (m in 1:M){H.mc.init$Z[[m]]$row <-mySampler[[m]]$indMemberships$row}}
if(!paramsMCMC$opEchan$ZCol){for (m in 1:M){H.mc.init$Z[[m]]$col <-mySampler[[m]]$indMemberships$col}}


resMCMC <- MCMCKernel(mydata, H.mc.init, alpha.t = 1, hyperparamPrior,hyperparamApproxPost = NULL, emissionDist = 'bernoulli', model =  'piColBipartiteSBM',paramsMCMC, opSave= TRUE)

burnin  = 2500 
extr <- burnin:paramsMCMC$nbIterMCMC

m = sample(1:M,1)
apply(resMCMC$seqZ[[m]]$row[,,extr],c(1,2),mean)[1:10,]
resMCMC$seqZ[[m]]$row[,,1][1:10,]
hyperparamApproxPost$collecTau[[m]]$row[1:10,]


apply(resMCMC$seqZ[[m]]$col[,,extr],c(1,2),mean)[1:10,]
resMCMC$seqZ[[m]]$col[,,1][1:10,]
hyperparamApproxPost$collecTau[[m]]$col[1:10,]

par(mfrow=c(KRow,KCol))
for (k in 1:KRow){
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
  for (k in 1:KRow){
    if (k ==1){plot(density(resMCMC$seqBlockProp$row[k,extr]),main='pi_k',xlim = c(0,1),col='green')
      }else{
        lines(density(resMCMC$seqBlockProp$row[k,extr]),col='green')
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp$row[k], sum(hyperparamApproxPost$blockProp$row)-hyperparamApproxPost$blockProp$row[k]),col='red',add=TRUE)
    abline(v = blockPropTrue$row[1,k])
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
    for (k in 1:KRow){
      if (k ==1){plot(density(resMCMC$seqBlockProp$row[m,k,extr]),main='pi_k',xlim = c(0,1),col='green')
    }else{
      lines(density(resMCMC$seqBlockProp$row[m,k,extr]),col='green')
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp$row[m,k], sum(hyperparamApproxPost$blockProp$row[m,])-hyperparamApproxPost$blockProp$row[m,k]),col='red',add=TRUE)
    abline(v = blockPropTrue$row[m,k])
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




