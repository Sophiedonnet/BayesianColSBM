rm(list=ls())
#------------------ set wd
if(Sys.info()[[4]]=='sophie-Latitude-5310'){
  setwd('/home/sophie/WORK_LOCAL/RECHERCHE/TRAVAUX_DE_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}
if(Sys.info()[[4]]=="donnet-Precision-Tower-5810"){
  setwd('/home/donnet/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianColSBM')
}
pathtodata <- paste0(getwd(),'/dataThomson')
#--------------------------------- 
library(sbm)
library(gtools)
library(DescTools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)

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

plot(resSMC$alpha.vec,type='l')


estimOptionsSMC$op.SMC.classic <- TRUE
resSMC_classic <- SMCColSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)
save(mydata,resSMC,hyperparamApproxPost,resSMC_classic,file='resultsHerlzier.Rdata')
######################################



#####################################
##################  Proba d'être dans le même clusterpar VB
############################################ 
for (m in 1:M){
  plotMyMatrix(hyperparamApproxPost$collecTau[[m]] %*%t( hyperparamApproxPost$collecTau[[m]]))
}



VZZ <- hyperparamApproxPost$collecTau[[1]]%*%t(hyperparamApproxPost$collecTau[[1]])
ZSample <- lapply(1:M,function(m){resSMC$HSample_end$ZSample[[m]]})
postZSample <- transfZsampleIntoMatrix(ZSample)
UZZ  = matrix(0,nbNodes[1],nbNodes[1])
for (i in 1:nbNodes[1]){
  for (j in 1:nbNodes[1]){
    UZZ[i,j] = sum(resSMC$W.end*(postZSample[[1]][i,]== postZSample[[1]][j,])) 
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

mean(abs(UZZ-1/2))
mean(abs(VZZ-1/2))

EntropyBernoulli = function(p){
  -p * log2(p)  - (1-p)*log2 (1-p)
}

