


#############################################################################
# MCMC kernel
#############################################################################


MCMCKernel <- function(data, H.mc.init, alpha.t, hyperparamPrior,hyperparamApproxPost, emissionDist, model,paramsMCMC, opSave=FALSE){
  ### 
  # Data 
  # H.mc.init : valeurs initiales des variables latentes et paramètres
  # alpha.t : param?tres de pond?ration post, post approch?e
  # hyperparamApproxPost : Parametres de la loi a priori
  # hyperparamPrior : parametres de la loi a posteriori approch?e
  # emissionDist : poisson or bernoulli
  # model:  piColBipartiteSBM or iidColBipartiteSBM
  # nbIterMCMC :  Number of it?rations
  
  ### variables latentes et param?tres
  if( !(model %in% c('piColSBM','iidColSBM'))){stop('Mispecified model')}
  if( !(emissionDist %in% c('poisson','bernoulli'))){stop('Mispecified emission distribution')}
  ##### 
  M <- data$M
  collecNetworks <- data$collecNetworks
  nbNodes <- data$nbNodes
  K <- nrow(hyperparamPrior$connectParam$alpha)
  H.mc <- H.mc.init
  
  ###### Normal MCMC
  if (is.null(hyperparamApproxPost)){
    hyperparamApproxPost <- hyperparamPrior
    hyperparamApproxPost$collecTau <- vector('list',data$M)
    for (m in 1:M){
      if(model=='iidColBipartiteSBM'){
        pi.m <- hyperparamPrior$blockProp
      }
      if (model=='piColBipartiteSBM'){
        pi.m <- hyperparamPrior$blockProp[m,]
      }
    hyperparamApproxPost$collecTau[[m]] <- matrix(pi.m ,nbNodes[m],K,byrow = TRUE)
    }
  }
  B <- paramsMCMC$nbIterMCMC
  if(is.null(paramsMCMC$opEchan)){
    opEchan <- list(connectParam = TRUE, blockProp = TRUE, Z = TRUE)
  }else{
    opEchan <- paramsMCMC$opEchan
  }
  
  
 
  # save
  if (opSave){
    seqConnectParam = array(0,c(K,K,1+B))
    seqConnectParam[,,1] = H.mc$connectParam
    seqZ  = lapply(1:M,function(m){array(0,c(nbNodes[m],K,1+B))})
    for (m in 1:M){
      seqZ[[m]][,,1] = H.mc$Z[[m]]
    }
    if(model=='iidColBipartiteSBM'){
      seqBlockProp = matrix(0,K,1+B)
      seqBlockProp[,1] = H.mc$blockProp
    }
    if(model=='piColBipartiteSBM'){
      seqBlockProp = array(0,c(M,K,1+ B))
      seqBlockProp[,,1] = H.mc$blockProp
    }
  }
    
    
  
  for (iterMCMC in 2:(1+B)){
  
    ################################################
    ############# simulation of connectParam
    ################################################
    
    if (opEchan$connectParam == 1) {
      
      S  <- lapply(1:M,function(m){t(H.mc$Z[[m]]) %*% matrix(1,nbNodes[m],nbNodes[m]) %*% H.mc$Z[[m]]})
      SY  <- lapply(1:M,function(m){t(H.mc$Z[[m]]) %*% collecNetworks[[m]] %*% H.mc$Z[[m]]})
      alpha_sim <- alpha.t*(hyperparamPrior$connectParam$alpha +  Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$alpha
      beta_sim <- alpha.t*(hyperparamPrior$connectParam$beta +  Reduce(`+`, S) - (emissionDist == 'bernoulli') * Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$beta
      if(emissionDist == 'bernoulli'){H.mc$connectParam <- matrix(rbeta(K^2,alpha_sim,beta_sim),K,K)}
      if(emissionDist == 'poisson'){H.mc$connectParam <- matrix(rgamma(K^2,alpha_sim,beta_sim),K,K)}
      
      if(opSave){seqConnectParam[,,iterMCMC] <- H.mc$connectParam}
    }
  
    ################################################
    ############# simulation of blockParam
    ################################################
    
    if (opEchan$blockProp == 1) {
      
      S  <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]])}))

      if(model=='iidColBipartiteSBM'){
        e <- alpha.t*(colSums(S) + hyperparamPrior$blockProp) + (1 - alpha.t) * hyperparamApproxPost$blockProp
        H.mc$blockProp <- c(rdirichlet(1,e))

        if(opSave){
          seqBlockProp[,iterMCMC] = H.mc$blockProp 
        }
      }
      if(model=='piColBipartiteSBM'){
        e <- alpha.t*(S + hyperparamPrior$blockProp) + (1 - alpha.t) * hyperparamApproxPost$blockProp
        H.mc$blockProp<- t(sapply(1:M,function(m){rdirichlet(1,e[m,])}))
         if(opSave){
          seqBlockProp[,,iterMCMC] = H.mc$blockProp 
       
        }
      }
    }
    ################################################
    ############# simulation of Z 
    ################################################
    #browser
    orderNetworks <-  sample(1:M,M,replace = FALSE)
    for (m in orderNetworks){ # for all networks
      
      
      ############ Zrow
      if (opEchan$ZRow){  
        if(model=='iidColBipartiteSBM'){piRow.m <- H.mc$blockProp}
        if(model=='piColBipartiteSBM'){piRow.m <- H.mc$blockProp[m,]}
        logZRow <-  matrix(log(piRow.m),nrow=nbNodes[m,1],ncol=KRow,byrow = TRUE)  
        
        if(emissionDist == 'bernoulli'){
          logYRow <- collecNetworks[[m]]%*%H.mc$Z[[m]]$col %*%t(log(H.mc$connectParam)) +   (1-collecNetworks[[m]])%*%H.mc$Z[[m]]$col %*%t(log(1-H.mc$connectParam))
        }
        if(emissionDist == 'poisson'){
          logYRow <- collecNetworks[[m]]%*%H.mc$Z[[m]]$col %*%t(log(H.mc$connectParam))  - matrix(1,nbNodes[m,1],nbNodes[m,2]) %*%  H.mc$Z[[m]]$col %*% t(H.mc$connectParam)
        }
        logProbZRow  <- alpha.t * (logYRow + logZRow) + (1-alpha.t)*log(hyperparamApproxPost$collecTau[[m]] )
        ProbZRow <- fromBtoTau(logProbZRow, eps = 10^-10)
        H.mc$Z[[m]] <- t(sapply(1:nbNodes[m,1],function(i){rmultinom(1,size=1,prob = ProbZRow[i,])}))
        if(opSave){
          seqZ[[m]][,,iterMCMC] = H.mc$Z[[m]]
        }
      }
      
    }
    
      
      
  }
  res <- H.mc
  if(opSave){res = list(H.mc = H.mc,seqConnectParam = seqConnectParam,seqBlockProp = seqBlockProp, seqZ = seqZ)}
  return(res)
}
      
# cat('Fin MCM.Kernel')
