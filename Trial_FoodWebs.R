rm(list=ls())
if(Sys.info()[[4]]=='sophie-Latitude-5310'){
  setwd('/home/sophie/WORK_LOCAL/RECHERCHE/TRAVAUX_DE_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}
if(Sys.info()[[4]]=="donnet-Precision-Tower-5810"){
  setwd('/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}
############################ Files and data
library(sbm)
library(gtools)
library(DescTools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)

#------------------------- all data
pathtodata <- paste0(getwd(),'/dataThomson/data')
my_files <- list.files(pathtodata)
food_webs <- lapply(
  X = seq_along(my_files),
  FUN = function(i) {
    df <- read.csv(file = paste0(pathtodata,"/",my_files[i]), header = TRUE, row.names = 1)
    A <- as.matrix(df)
    return(list(
      net = A,
      nr = nrow(A),
      nc = ncol(A),
      dens = mean(A),
      id = stringr::str_sub(my_files[i], 1, -5))
    )
  }
)
site_names <- c("Martins(M)", "Cooper(NC)", "Herlzier(NC)", "Venlaw(NC)",
                "Berwick(NZ)", "NorthCol(NZ)", "Powder(NZ)", "Trib(NZ)" )
names(food_webs)<- site_names

#########################################################"
#---------------------- What to do 
#########################################################"

net <- 'simple' # collec
net <- 'collec'
estim <- TRUE # 'No' # "Simple",'iidColSbm'
model <- 'piColSBM'

#------------------ 

if(net=='simple'){
  pathtores <- paste0(getwd(),'/dataThomson/res_estim_simple')
  extr = 3
}else{
  if(model=='iidColSbm'){
  pathtores <- paste0(getwd(),'/dataThomson/res_estim_iidColSbm')
  }
  if(model=='piColSbm'){
    pathtores <- paste0(getwd(),'/dataThomson/res_estim_piColSbm')
  }
  extr = 1:3
}
#---------------- collec networks
collecNetworks  = list()
nbNodes <- rep(0,length(extr))
for (i in 1:length(extr)){
  name_i <- site_names[extr[i]]
  Li <- food_webs[[extr[i]]]
  Mat <- Li$net
  nr <- Li$nr
  U <- order(rownames(Mat))
  Mat <- Mat[U,U]
  collecNetworks[[i]] <- Mat
  nbNodes[i]= nr
}
names(collecNetworks) <- site_names[extr]
mydata <- list(collecNetworks  = collecNetworks, M = length(extr), nbNodes = nbNodes)

#################### 
### Model 
########################

emissionDist = 'bernoulli'





############################################################
######################## VB ###############################
#############################################################
if(net == 'simple'){
  names_files_res <- paste0(pathtores,'/res_estim_',site_names[i],'_VB.Rdata')
}
if(net=='collec'){
  names_files_res <- paste0(pathtores,'/res_estim_VB.Rdata')
}


if(estim){
  
  ###########################################################################################
  #------------------ VBEM 
  ###########################################################################################
  #--------------- init CollecTau
  initSimple <- lapply(collecNetworks,estimateSimpleSBM)
  myRef = which.max(sum(t(sapply(initSimple,function(m){m$nbBlocks}))))
  collecTau_init <- initCollecTau(initSimple, ref = myRef)
  KEstim <- initSimple[[myRef]]$nbBlocks[1]
  print(c(KEstim))
  
  #---------------  Prior distribution
  
  hyperparamPrior <- setHyperparamPrior(mydata$M,KEstim, emissionDist, model)
  
  
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

  save(hyperparamPrior,collecTau_init,hyperparamApproxPost,names_files_res,file = names_files_res )
}else{
  load(file = names_files_res)
}

w <- apply(hyperparamApproxPost$collecTau[[1]],1,which.max)
mu_post <- hyperparamApproxPost$connectParam$alpha/(hyperparamApproxPost$connectParam$alpha+hyperparamApproxPost$connectParam$beta)

o <- order(rowSums(mu_post))
plotMyMatrix(mu_post[o,o])

mymat <- 1*(mu_post>0.004)

plotMyMatrix(mymat)
o <- c(4,,3,1,5)
plotMyMatrix(mu_post[o,o])

P <- permutations(5, 5)
score <- sum(mu_post*lower.tri(mu_post))
best_o <- 1:5 
test_score <- rep(0,nrow(P))
test_score[1] <-score
for(p in 2:nrow(P)){
  print(p)
  o <- P[p,]
  M <- 1 * (mu_post[o,o]>0.01)
  plotMyMatrix(M)
  test_score[p]  <- mean(M*lower.tri(M,diag = FALSE))
}

my_o <- order(test_score,decreasing = TRUE)
list_PLot <- list()
i= 0
for (i in 1:10){
  i = i+1 
  u <- P[my_o[i],]
  plotMyMatrix(mu_post[u,u])
  
}

###########################################################################################
#------------------  SMC 
###########################################################################################


if(net == 'simple'){
  names_files_res <- paste0(pathtores,'/res_estim_',site_names[i],'_SMCVB.Rdata')
}else{
  names_files_res <- paste0(pathtores,'/res_estim_SMCVB.Rdata')
}
if(estim){
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
  save(estimOptionsSMC, hyperparamPrior,hyperparamApproxPost,resSMC,names_files_res,file = names_files_res )
}else{
  load(file=names_files_res)
}
plot(resSMC$alpha.vec,type='l')

######################################



#####################################
##################  Proba d'être dans le même clusterpar VB
############################################
M = 1
ZZ_VB <- lapply(1:M,function(m){hyperparamApproxPost$collecTau[[m]]%*%t(hyperparamApproxPost$collecTau[[m]])})

ZSample <- lapply(1:M,function(m){resSMC$HSample_end$ZSample[[m]]})
postZSample <- transfZsampleIntoMatrix(ZSample)
ZZ_SMC <- lapply(1:M,function(m){
  MZZ <- matrix(0,nbNodes[m],nbNodes[m])
  for (i in 1:nbNodes[m]){
    for (j in 1:nbNodes[m]){
      MZZ [i,j] = sum(resSMC$W.end*(postZSample[[m]][i,]== postZSample[[m]][j,])) 
    }
  }
  return(MZZ)
}
)



#ZSample <- lapply(1:M,function(m){resSMC_classic$HSample_end$ZSample[[m]][,,extr]})
#postZSample <- transfZsampleIntoMatrix(ZSample)
#ZZ_SMC_classic  = matrix(0,nbNodes[1],nbNodes[1])
#for (i in 1:nbNodes[1]){
#  for (j in 1:nbNodes[1]){
#    ZZ_SMC_classic [i,j] = sum(resSMC_classic$W.end[extr]*(postZSample[[1]][i,]== postZSample[[1]][j,])) 
#  }
#}

UVB <- ZZ_VB[[1]][upper.tri(ZZ_VB[[1]])]
USMC <- ZZ_SMC[[1]][upper.tri(ZZ_SMC[[1]])]
ZZ_DF <- cbind(UVB,USMC)
names(ZZ_DF) <- c('VB','VBSMC')
ggplot(ZZ_DF) + geom_points(aes(x=VB,y=VBSMC))
par(mfrow=c(1,2))
for (m in 1:M){
  plot(ZZ_VB[[m]][w],ZZ_SMC[[m]][w])
  plot(ZZ_VB[[m]],ZZ_SMC[[m]])
  
  abline(a=0,b=1)
  
  hist(ZZ_SMC[[m]])
  hist(ZZ_VB[[m]],add= TRUE,col='blue')
  
  
  mean( abs(ZZ_SMC[[m]]-1/2))
  mean( abs(ZZ_VB[[m]]-1/2))
  
  hist(abs(1/2- ZZ_VB[[m]]),add= TRUE,col='blue')
  
  
}
ZZ_all <- lower.tri(ZZ_SMC)*ZZ_SMC + upper.tri(ZZ_VB)*ZZ_VB
plotMyMatrix(ZZ_all)
plotMyMatrix(ZZ_VB)
plotMyMatrix(ZZ_SMC)

hist(ZZ_VB)
hist(ZZ_SMC,add = TRUE,col="blue")
###############################
K <- KEstim
par(mfrow=c(1,1))
for (k in 1:K){
  # 
 #if((k ==1)){
      plot(density(apply(resSMC_classic$HSample_end$connectParamSample,c(3),sum),weights=resSMC_classic$W.end),main='alpha',col='green')
  #}
  #else{
  #  lines(density(resSMC$HSample_end$connectParamSample[k,k,],weights=resSMC$W.end),col='green')
  #}
  curve(dbeta(x,hyperparamApproxPost$connectParam$alpha[k,k],hyperparamApproxPost$connectParam$beta[k,k]),col='red',add=TRUE)
  lines(density(apply(resSMC$HSample_end$connectParamSample,c(3),sum),weights=resSMC$W.end),col='magenta')
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


plot(resSMC$HSample_end$blockPropSample)



mean(abs(ZZ_VB-1/2))
mean(abs(ZZ_SMC-1/2))
mean(abs(ZZ_SMC_classic-1/2))



