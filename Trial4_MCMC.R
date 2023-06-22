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


M = 2
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
nbNodes[1] = c(40)
nbNodes[2] = 30
#nbNodes[5,] = c(30,20)

#---------  SIMULATION

mySampler <- lapply(1:M, function(m){sampleSimpleSBM(nbNodes = nbNodes[m], blockProp =  blockPropTrue[m,],  connectParamTrue,directed = TRUE)})
collecNetworks <- lapply(mySampler,function(l){Net.m <- l$networkData; diag(Net.m) <- 0; Net.m})
propNoise = 0.2
noisyCollecNetworks <- noiseSampling(collecNetworks,propNoise)
rm(collecNetworks)
mydata <- list(collecNetworks = noisyCollecNetworks, M= M, nbNodes = nbNodes)


#######################################################
#--------------- init CollecTau
##########################################################
initSimple <- lapply(mydata$collecNetworks,estimateSimpleSBM)
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
resEstimVBEM  <- VBEMColSBM(mydata,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)

hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau




##################################


HSample <- rParamZ(1, hyperparam = hyperparamApproxPost, emissionDist, model ,nbNodes)
mc <- 1
H.mc.sim <- list(connectParam = HSample$connectParamSample[,,mc])
H.mc.sim$Z <- lapply(1:M, function(m){HSample$ZSample[[m]][,,mc]})
if (model == 'iidColSBM'){
  H.mc.sim$blockProp <- HSample$blockPropSample[,1]
}
if (model == 'piColSBM'){
  H.mc.sim$blockProp <-  HSample$blockPropSample[,,1]
  if(M==1){
    H.mc.sim$blockProp <-  matrix(HSample$blockPropSample[,,1],nrow=1)
  }
}

  
###################### Estim avec Kcol K true

paramsMCMC = list(nbIterMCMC = 5000)
H.mc.init <- H.mc.sim
paramsMCMC$opEchan = list(connectParam = TRUE, blockProp = TRUE, Z = TRUE)
if (!paramsMCMC$opEchan$connectParam){H.mc.init$connectParam <- connectParamTrue$mean}
if (!paramsMCMC$opEchan$blockProp){H.mc.init$blockProp <- blockPropTrue}
if(!paramsMCMC$opEchan$Z){
  for (m in 1:M){
    Z.m <- mySampler[[m]]$memberships
    n.m <- length(Z.m)
    H.mc.init$Z[[m]] = matrix(0,nrow = n.m,ncol = K) 
    for (i in 1:n.m){H.mc.init$Z[[m]][i,Z.m[i]] = 1}
  }
}
paramsMCMC$opPrint = TRUE
alpha.t = 1
data = mydata
emissionDist = 'bernoulli'
model =  'piColSBM'
opSave = TRUE

resMCMC <- MCMCKernel(mydata, H.mc.init, alpha.t = 1, hyperparamPrior,hyperparamApproxPost = NULL, emissionDist = 'bernoulli', model =  'piColSBM',paramsMCMC, opSave= TRUE)

burnin  = 1000 
extr <- seq(burnin,paramsMCMC$nbIterMCMC,by=5)

m = sample(1:M,1)
apply(resMCMC$seqZ[[m]][,,extr],c(1,2),mean)[1:5,]
hyperparamApproxPost$collecTau[[m]][1:5,]

#####################################################################""" 
#######################################   PLot
###################################################################### 


par(mfrow=c(K,K))
for (k in 1:K){
  for (l in 1:K){
    plot(resMCMC$seqConnectParam[k,l,extr],main='alpha',col='green',type='l')
    abline(h = (1-propNoise)*connectParamTrue$mean[k,l] + propNoise*(1-connectParamTrue$mean[k,l]))
  }
}


par(mfrow=c(K,K))
for (k in 1:K){
  for (l in 1:K){
    plot(density(resMCMC$seqConnectParam[k,l,extr]),main=c(k,l),xlim=c(0,1),col='green')
    curve(dbeta(x,hyperparamApproxPost$connectParam$alpha[k,l],hyperparamApproxPost$connectParam$beta[k,l]),col='red',add=TRUE)
    abline(v = (1-propNoise)*connectParamTrue$mean[k,l] + propNoise*(1-connectParamTrue$mean[k,l]))
  }
}

#------------------------------------ 
if(model=="iidColSBM"){
  par(mfrow=c(1,1))
  for (k in 1:K){
    if (k ==1){plot(density(resMCMC$seqBlockProp[k,extr]),main='pi_k',xlim = c(0,1),col='green')
      }else{
        lines(density(resMCMC$seqBlockProp[k,extr]),col='green')
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp[k], sum(hyperparamApproxPost$blockProp)-hyperparamApproxPost$blockProp[k]),col='red',add=TRUE)
    abline(v = blockPropTrue[1,k])
  }
}

if(model=="piColSBM"){
  par(mfrow=c(floor(M/2)+1,2))
  if(M==1){par(mfrow=c(1,1))}
  for (m in 1:M){
    for (k in 1:K){
      if (k ==1){plot(density(resMCMC$seqBlockProp[m,k,extr]),main='pi_k',xlim = c(0,1),col='green')
    }else{
      lines(density(resMCMC$seqBlockProp[m,k,extr]),col='green',lty=2,lwd=2)
    }  
    curve(dbeta(x,hyperparamApproxPost$blockProp[m,k], sum(hyperparamApproxPost$blockProp[m,])-hyperparamApproxPost$blockProp[m,k]),col='red',add=TRUE)
    abline(v = blockPropTrue[m,k])
    }
  }  
}


#####################################
##################  Proba d'être dans le même clusterpar VB

for (m in 1:M){
  plotMyMatrix(hyperparamApproxPost$collecTau[[m]] %*%t( hyperparamApproxPost$collecTau[[m]]))
}
VZZ <- hyperparamApproxPost$collecTau[[1]]%*%t(hyperparamApproxPost$collecTau[[2]])

ZSample <- lapply(1:M,function(m){resMCMC$seqZ[[m]][,,extr]})
postZSample <- transfZsampleIntoMatrix(ZSample)


VZZ <- hyperparamApproxPost$collecTau[[1]]%*%t(hyperparamApproxPost$collecTau[[1]])

UZZ  = matrix(0,nbNodes[1],nbNodes[1])
for (i in 1:nbNodes[1]){
  for (j in 1:nbNodes[1]){
  UZZ[i,j] = mean(postZSample[[1]][i,]== postZSample[[1]][j,]) 
  }
}


UZZ  = matrix(0,nbNodes[2],nbNodes[2])
for (i in 1:nbNodes[2]){
  for (j in 1:nbNodes[2]){
    UZZ[i,j] = mean(postZSample[[2]][i,]== postZSample[[2]][j,]) 
  }
}
VZZ <- hyperparamApproxPost$collecTau[[2]]%*%t(hyperparamApproxPost$collecTau[[2]])
plotMyMatrix(UZZ)
plotMyMatrix(VZZ)

# 1/2 = bordel absolu. 
mean(abs(UZZ-1/2)) # MCMC 
mean(abs(VZZ-1/2)) # VB 
hist(abs(UZZ-1/2))
hist(abs(VZZ-1/2))

# 
EntropyBernoulli = function(p){
  -p * log2(p)  - (1-p)*log2 (1-p)
}



     