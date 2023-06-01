VBEMColSBM = function(collecNetworks,hyperparamPrior,collecTau,estimOptions, emissionDist, model ){

  
  if(is.null(estimOptions)){
    estimOptions <- list(maxIterVB = 100,
                        maxIterVE = 100,
                        valStopCritVE = 10^-5,
                        valStopCritVB = 10^-5)
  } 

  
  
  
  M <- length(collecNetworks)
  nbNodes <- sapply(collecNetworks,function(X)(nrow(X)))
  K <-nrow(hyperparamPrior$connectParam$alpha)
  
  #-------------------- initialisation
  hyperparamPost <- Mstep(collecNetworks,M, nbNodes,collecTau,hyperparamPrior, emissionDist, model)
  noConvergence <- 0;
  iterVB <- 0
  stopVB <- 0 
  
  #-------------- RUN VBEM
  while (iterVB < estimOptions$maxIterVB & stopVB == 0){
  
    iterVB  <- iterVB  + 1
    hyperparamPostOld <- hyperparamPost
  
    ##--------------- VE step
    res_Estep <- Estep(collecNetworks,M, nbNodes,K,collecTau,hyperparamPost,estimOptions,emissionDist, model)
    collecTau <- res_Estep$collecTau
    #----- end  VB M step
    
    ##--------------- VB step
    hyperparamPost <- Mstep(collecNetworks,M, nbNodes,collecTau,hyperparamPrior, emissionDist, model)
    #----- end  VB M step
    
    ##-------------- stop criteria
    if (disthyperparamPost(hyperparamPost,hyperparamPostOld) < estimOptions$valStopCritVB) {stopVB <- 1}
    print(iterVB)
  }
  
  return(reorderBlocks(hyperparamPost, collecTau, model))
    

}


############################################################
####################" MSTEP from CollecTau to hyperparamPost
###############################################################
Mstep <- function(collecNetworks,M, nbNodes,collecTau,hyperparamPrior,emissionDist, model){
  
  
  hyperparamPost <- hyperparamPrior
  S <- lapply(1:M,function(m){t(collecTau[[m]]) %*% matrix(1,nbNodes[m],nbNodes[m]) %*% collecTau[[m]]})
  SY <- lapply(1:M,function(m){t(collecTau[[m]]) %*% collecNetworks[[m]] %*% collecTau[[m]]})
  
  
  hyperparamPost$connectParam$alpha  <- hyperparamPrior$connectParam$alpha +  Reduce(`+`, SY)
  hyperparamPost$connectParam$beta  <- hyperparamPrior$connectParam$beta +   Reduce(`+`, S) - (emissionDist == 'bernoulli') * Reduce(`+`, SY)
  
  Rtau <-  t(sapply(1:M,function(m){colSums(collecTau[[m]])}))

  if(model == 'iidColSBM'){  
    hyperparamPost$blockProp  <- hyperparamPrior$blockProp + apply(Rtau,2,sum)
  }
  if(model == 'piColSBM'){
    hyperparamPost$blockProp   <- hyperparamPrior$blockProp  + Rtau
  }
  return(hyperparamPost)
}

############################################################
####################" MSTEP from CollecTau to hyperparamPost
###############################################################
Estep <- function(collecNetworks,M, nNodes,K,collecTau,hyperparamPost,estimOptions,emissionDist, model){
  
  iterVE <- 0
  stopVE <- 0
  
  while ((iterVE < estimOptions$maxIterVE) & (stopVE == 0)){
    collecTauOld <- collecTau #useful ?
    if(K  == 1){for (m in 1:M){collecTau[[m]] <-  matrix(1,ncol  = 1,nrow = nRow[m])}}
    if (K>1) {
    ## useful quantities
      DiGAlpha <- digamma(hyperparamPost$connectParam$alpha) ### useful for Poisson and Bernoulli
      if(emissionDist == 'bernoulli'){
        DiGBeta <- digamma(hyperparamPost$connectParam$beta) 
        DiGAlphaBeta  <- digamma(hyperparamPost$connectParam$beta + hyperparamPost$connectParam$alpha) 
      }
      
      for(m in 1:M){
        if(model == 'iidColSBM'){
          DiblockProp_m <- digamma(hyperparamPost$blockProp) - digamma(sum(hyperparamPost$blockProp))
        }
        if(model == 'piColSBM'){
          DiblockProp_m <- digamma(hyperparamPost$blockProp[m,]) - digamma(sum(hyperparamPost$blockProp[m,]))
        }
        l3_m  <- matrix(DiblockProp_m,nrow = nbNodes[m],ncol = K,byrow = TRUE)
        if(emissionDist == 'poisson'){
            lY_m  <- collecNetworks[[m]] %*% tcrossprod(collecTau[[m]],log(hyperparamPost$connectParam$beta) + DiGAlpha)  
            lY_m <- lY_m - matrix(1,nbNBodes[m], nbNodes[m])  %*% tcrossprod(collecTau[[m]],hyperparamPost$connectParam$alpha/hyperparamPost$connectParam$beta)
        }
        if(emissionDist == 'bernoulli'){
          lY_m  <- collecNetworks[[m]] %*% tcrossprod(collecTau[[m]],DiGAlpha- DiGAlphaBeta)  +  (1-collecNetworks[[m]]) %*% tcrossprod(collecTau[[m]],DiGBeta- DiGAlphaBeta)
        }
        
        collecTau[[m]] <- fromBtoTau(lY_m + l3_m) 
      }

    deltaTau <- distTau(collecTau,collecTauOld)
    if (deltaTau < estimOptions$valStopCritVE) {stopVE <- 1}
    iterVE <- iterVE + 1
    noConvergence <- 1*(iterVE  == estimOptions$maxIterVE)
  }
  
  return(list(collecTau = collecTau,noConvergence  = noConvergence))
  }
  
  

############################## log marg likelihood
computeLogLikMarg_VB <- function(collecNetworks, collecTau, hyperparamPrior, emissionDist, model){
  
  M <- length(collecNetworks)
  collecZMAP <- lapply(collecTau,function(tau.m){
    indZ <- tau.m
    indZ <- t(sapply(1:nrow(tau.m),function(i){u <- 0*tau.m[i,]; u[which.max(tau.m[i,])]=1; return(u)}))
    if(min(dim(indZ))==1){indZ = matrix(indZ,ncol=1)}
    return(indZ)
  })
  
  if(emissionDist =='bernoulli'){
    S1 <- lapply(1:M,function(m){t(collecZMAP[[m]])%*% collecNetworks[[m]] %*% collecZMAP[[m]]})
    S0 <- lapply(1:M,function(m){t(collecZMAP[[m]])%*% (1-collecNetworks[[m]]) %*% collecZMAP[[m]]})
    ahat <- Reduce('+',S1)
    bhat <- Reduce('+',S0)
    lprobYZ <- sum(lbeta(ahat+hyperparamPrior$connectParam$alpha,bhat+hyperparamPrior$connectParam$beta) -lbeta(hyperparamPrior$connectParam$alpha,hyperparamPrior$connectParam$beta))
  }
  
  if(model=='poisson'){
    lprobYZ  = NA
  }
  
  if(model == 'piColSBM'){
    lprobZ <- Reduce("+", 
                     lapply(1:M,
                            function(m){
                              dPost<- apply(collecZMAP[[m]],2,sum) 
                              L   <- mylBetaFunction(dPost + hyperparamPrior$blockProp[m,]) - mylBetaFunction(hyperparamPrior$blockProp[m,])
                              return(L)
                            }
                     )
    )
  }
  if(model == 'iidColSBM'){
    dPost <- rowSums(sapply(1:M,function(m){apply(collecZMAP[[m]],2,sum)}))
    L  <- mylBetaFunction(dPost  + hyperparamPrior$blockProp) - mylBetaFunction(hyperparamPrior$blockProp )
    lprobZ <-  L
  }
  
  res <- c(lprobYZ,lprobZ)
  names(res) =c('lprobYZ','lprobZ')
  return(res)
}

