
##################### divide by its sum each row (to obtain probabilities)
normByRow <- function(Mat){
  U <- t(scale(t(Mat), center = FALSE, scale = colSums(t(Mat))))
  attr(U,"scaled:scale")<- NULL
  return(U)
}

##################### IN VE step, to obtain the tau  
fromBtoTau <- function(B,eps = 10^-5){
  B <- B - matrix(apply(B,1,max),nrow = nrow(B),ncol = ncol(B),byrow = FALSE)
  temp <- exp(B)
  temp2 <- normByRow(temp)
  temp2[temp2 < eps] <- eps
  temp2[temp2 > (1 - eps)] <- 1 - eps
  Tau <- normByRow(temp2)
  attr(Tau,"scaled:scale")<- NULL
  return(Tau)
}

##################### Difference between tau and tauOld
distTau  <- function(collecTau,collecTauOld)
{
  M <- length(collecTau)
  vdis <- sapply(1:M,function(m){
    d  <- sqrt(sum(as.vector(collecTau[[m]] - collecTauOld[[m]])^2))
    return(d)
  })
  return(sum(vdis))
}

##################### Difference between paramPostOld and paramPost

disthyperparamPost <- function(hyperparamPost,hyperparamPostOld){
  
  d1 <- sqrt(sum(as.vector((hyperparamPost$connectParam$alpha - hyperparamPostOld$connectParam$alpha)^2)))
  d2 <- sqrt(sum(as.vector((hyperparamPost$connectParam$beta -  hyperparamPostOld$connectParam$beta )^2)))
  d3 <- sqrt(sum(as.vector((hyperparamPost$blockProp -  hyperparamPostOld$blockProp )^2)))
  return(d1 + d2 + d3 )
}


####################### Reorder groups row / cols. 
reorderBlocks  <- function(hyperparamPost, collecTau ,model){
  
  M <- length(collecTau)
  meanPost <- hyperparamPost$connectParam$alpha / (hyperparamPost$connectParam$alpha + hyperparamPost$connectParam$beta)
  alphaDiag <- diag(meanPost); 
  ord <- order(alphaDiag,decreasing = TRUE)
  hyperparamPost$connectParam <- lapply(hyperparamPost$connectParam,function(u){u[ord,ord]})
  if(model == 'iidColSBM'){  
    hyperparamPost$blockProp  <- hyperparamPost$blockProp[ord]
  }
  if(model == 'piColSBM'){
    hyperparamPost$blockProp <- hyperparamPost$blockProp[,ord]
    if(M==1){
      hyperparamPost$blockProp  <- matrix(hyperparamPost$blockProp ,nrow=1)
    }
  }
  
  
  collecTau <- lapply(collecTau ,function(tau){
    tauNew <- tau
    tauNew <-tau[,ord]
    return(tauNew)}
    )
  return(list(hyperparamPost = hyperparamPost,collecTau  = collecTau))
  
  
}


#######################""
mylBetaFunction = function(alpha){
   sum(lgamma(alpha))-lgamma(sum(alpha))
}


################################## set prior param
setHyperparamPrior <- function(M,K,emissionDist, model){
  # M  : number of  networks
  # K  : number of blocks in row
  # Kcol  : number of blocks in col
  # emissionDistr : poisson or bernoulli
  
  hyperparamPrior <- list()
  
  
  
  hyperparamPrior$blockProp  
  if (model == 'iidColSBM'){
    hyperparamPrior$blockProp  = rep(1/K,K) ### Jeffreys dirichlet prior the pi and rho
}
  if (model == 'piColSBM'){
    hyperparamPrior$blockProp = matrix(1/K,M,K) ### Jeffreys dirichlet prior the pi^m and rho^m
  }
  
  hyperparamPrior$connectParam<- list()
  if (emissionDist =='bernoulli'){
    hyperparamPrior$connectParam$alpha <- matrix(1,K,K) ### Uniform prior on the alpha_{kl}
    hyperparamPrior$connectParam$beta <- matrix(1,K,K)
  }
  if (emissionDist =='poisson'){
    hyperparamPrior$connectParam$alpha <- matrix(1,K,K) ### Exp de param 1/100
    hyperparamPrior$connectParam$beta <- matrix(1/100,K,K)
  }
  return(hyperparamPrior)
}

xlogx <- function(x){
  res = x
  res[x>0] <- x[x>0]*log(x[x>0])
}
entropyBernoulli = function(p){
  
  xlogx(p) + xlogx(1-p)
  
}


