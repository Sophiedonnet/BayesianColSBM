


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
  
  if(is.null(paramsMCMC$opPrint)){paramsMCMC$opPrint = TRUE}
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
    hyperparamApproxPost$collecTau <- vector('list',M)
    for (m in 1:M){
      if(model=='iidColSBM'){
        pi.m <- hyperparamPrior$blockProp
      }
      if (model=='piColSBM'){
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
    seqConnectParam = array(0,c(K,K,B))
    seqConnectParam[,,1] = H.mc$connectParam
    seqZ  = lapply(1:M,function(m){array(0,c(nbNodes[m],K,B))})
    #for (m in 1:M){
    #  seqZ[[m]][,,1] = H.mc$Z[[m]]
    #}
    if(model=='iidColSBM'){
      seqBlockProp = matrix(0,K,B)
      seqBlockProp[,1] = H.mc$blockProp
    }
    if(model=='piColSBM'){
      seqBlockProp = array(0,c(M,K, B))
      seqBlockProp[,,1] = H.mc$blockProp
    }
  }
  
  
  
  for (iterMCMC in 2:B){
    
    if ((paramsMCMC$opPrint) & (iterMCMC%%100==0)){ 
      print(iterMCMC)
    }
    ################################################
    ############# simulation of Z 
    ################################################
    #browser
    orderNetworks <-  sample(1:M,M,replace = FALSE)
    if (opEchan$Z){ 
      for (m in orderNetworks){ # for all networks
        
        
        n.m <- nbNodes[m]
        orderNodes.m <- sample(1:n.m,n.m)
        if(model=='iidColSBM'){log.pi.m <- log(H.mc$blockProp)}
        if(model=='piColSBM'){log.pi.m <- log(H.mc$blockProp[m,])}
        
        log.tau.m <- log(hyperparamApproxPost$collecTau[[m]])
        Y.m <- collecNetworks[[m]]
        IY.m <- 1-Y.m
        diag(IY.m) <- 0
        for (i in orderNodes.m){ ### for any node
          #print(i)
          if(emissionDist == 'bernoulli'){
            logYRow.i.1 <- matrix(Y.m[i,],nrow = 1) %*%H.mc$Z[[m]] %*%t(log(H.mc$connectParam)) +   matrix(IY.m[i,],nrow=1)%*%H.mc$Z[[m]] %*%t(log(1-H.mc$connectParam))
            logYRow.i.2 <- matrix(Y.m[,i],nrow = 1) %*%H.mc$Z[[m]] %*%log(H.mc$connectParam) +   matrix(IY.m[,i],nrow=1)%*%H.mc$Z[[m]] %*%log(1-H.mc$connectParam)
            logYRow.i <- logYRow.i.1 + logYRow.i.2
          }
          logProbZRow.i  <- alpha.t * (logYRow.i + log.pi.m) + (1-alpha.t)*log.tau.m[i,]
          ProbZRow.i <- c(fromBtoTau(matrix(logProbZRow.i,nrow = 1), eps = 10^-10))
          H.mc$Z[[m]][i,] <- c(rmultinom(1,size=1,prob = ProbZRow.i))
        } #end boucle sur nodes
      }# end boucle sur m 
    }# end if echanZ  = TRUE
    ################################################
    ############# simulation of connectParam
    ################################################
    
    if (opEchan$connectParam) {
      
      S  <- lapply(1:M,function(m){t(H.mc$Z[[m]]) %*% matrix(1,nbNodes[m],nbNodes[m]) %*% H.mc$Z[[m]]})
      SY  <- lapply(1:M,function(m){t(H.mc$Z[[m]]) %*% collecNetworks[[m]] %*% H.mc$Z[[m]]})
      alpha_sim <- alpha.t*(hyperparamPrior$connectParam$alpha +  Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$alpha
      beta_sim <- alpha.t*(hyperparamPrior$connectParam$beta +  Reduce(`+`, S) - (emissionDist == 'bernoulli') * Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$beta
      if(emissionDist == 'bernoulli'){H.mc$connectParam <- matrix(rbeta(K^2,alpha_sim,beta_sim),K,K)}
      if(emissionDist == 'poisson'){H.mc$connectParam <- matrix(rgamma(K^2,alpha_sim,beta_sim),K,K)}
      
    }
    
    ################################################
    ############# simulation of blockParam
    ################################################
    
    if (opEchan$blockProp) {
      
      S  <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]])}))
      
      if(model=='iidColSBM'){
        e <- alpha.t*(colSums(S) + hyperparamPrior$blockProp) + (1 - alpha.t) * hyperparamApproxPost$blockProp
        H.mc$blockProp <- c(rdirichlet(1,e))
      }
      if(model=='piColSBM'){
        e <- alpha.t*(S + hyperparamPrior$blockProp) + (1 - alpha.t) * hyperparamApproxPost$blockProp
        H.mc$blockProp<- t(sapply(1:M,function(m){rdirichlet(1,e[m,])}))
      }
    }
    
    
    
    #############################################""   
    if(opSave){
      seqConnectParam[,,iterMCMC] <- H.mc$connectParam
      if(model=='iidColSBM'){seqBlockProp[,iterMCMC] = H.mc$blockProp}
      if(model=='piColSBM'){seqBlockProp[,,iterMCMC] = H.mc$blockProp}
      for (m in 1:M){seqZ[[m]][,,iterMCMC] = H.mc$Z[[m]]}
    }
    
  }# end boucle sur iterMCMC
  
  res <- H.mc
  if(opSave){res = list(H.mc = H.mc,seqConnectParam = seqConnectParam,seqBlockProp = seqBlockProp, seqZ = seqZ)}
  return(res)
}

# cat('Fin MCM.Kernel')
