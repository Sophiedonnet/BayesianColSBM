rm(list=ls())
if(Sys.info()[[4]]=='sophie-Latitude-5310'){
  setwd('/home/sophie/WORK_LOCAL/RECHERCHE/TRAVAUX_DE_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}
if(Sys.info()[[4]]=="donnet-Precision-Tower-5810"){
  setwd('/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}
#--------------------------------- 
library(sbm)
library(gtools)
library(DescTools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)


#------------------ 
pathtodata <- paste0(getwd(),'/dataThomson')
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

#for (i in 1:length(site_names)){

i = 1 
name_i <- site_names[i]
Li <- food_webs[[i]]
Mat <- Li$net
nr <- Li$nr
U <- order(rownames(Mat))
Mat <- Mat[U,U]

#################### 
### Model 
########################

emissionDist = 'bernoulli'
model = 'piColSBM'
collecNetworks <- list(Mat)
M = 1
mydata <- list(collecNetworks  = collecNetworks, M= 1, nbNodes = nr)

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

hyperparamPrior <- setHyperparamPrior(M,KEstim, emissionDist, model)



############################################################
######################## VB ###############################
#############################################################
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

names_files_res <- paste0(pathtodata,'/res_estim_',name_i,'_VB.Rdata')
save(hyperparamPrior,collecTau_init,hyperparamApproxPost,names_files_res,file = names_files_res )




###########################################################################################
#------------------  SMC 
###########################################################################################
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

names_files_res <- paste0(pathtodata,'/res_estim_',name_i,'_SMCVB.Rdata')
save(estimOptionsSMC, hyperparamPrior,hyperparamApproxPost,resSMC,names_files_res,file = names_files_res )

plot(resSMC$alpha.vec,type='l')

rm(hyperparamApproxPost)
estimOptionsSMC$op.SMC.classic <- TRUE
resSMC_classic <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)
names_files_res <- paste0(pathtodata,'/res_estim_',name_i,'_SMCClassic.Rdata')
save(estimOptionsSMC, hyperparamPrior,resSMC_classic,names_files_res,file = names_files_res )



######################################



#####################################
##################  Proba d'être dans le même clusterpar VB
############################################
load( paste0(pathtodata,'/res_estim_',name_i,'_VB.Rdata'))
for (m in 1:M){
  plotMyMatrix(hyperparamApproxPost$collecTau[[m]] %*%t( hyperparamApproxPost$collecTau[[m]]))
}

extr <- which(apply(resSMC_classic$HSample_end$connectParamSample,c(3),sum)<4)
length(extr)

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
    ZZ_SMC_classic [i,j] = sum(resSMC_classic$W.end[extr]*(postZSample[[1]][i,]== postZSample[[1]][j,])) 
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



