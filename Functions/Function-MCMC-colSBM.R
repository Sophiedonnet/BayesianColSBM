


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
  
  M <- data$M
  collecNetworks <- data$collecNetworks
  nbNodes <- data$nbNodes
  KRow <- nrow(hyperparamPrior$connectParam$alpha)
  KCol <- ncol(hyperparamPrior$connectParam$alpha)
  H.mc <- H.mc.init
  
  ###### Normal MCMC
  if (is.null(hyperparamApproxPost)){
    hyperparamApproxPost <- hyperparamPrior
    hyperparamApproxPost$collecTau <- vector('list',data$M)
    for (m in 1:M){
      if(model=='iidColBipartiteSBM'){
        pi.row.m <- hyperparamPrior$blockProp$row
        pi.col.m <- hyperparamPrior$blockProp$col
      }
      if (model=='piColBipartiteSBM'){
        pi.row.m <- hyperparamPrior$blockProp$row[m,]
        pi.col.m <- hyperparamPrior$blockProp$col[m,]
      }
    hyperparamApproxPost$collecTau[[m]]$row <- matrix(pi.row.m ,nbNodes[m,1],KRow,byrow = TRUE)
    hyperparamApproxPost$collecTau[[m]]$col <- matrix(pi.col.m ,nbNodes[m,2],KCol,byrow = TRUE)
    }
  }
  B <- paramsMCMC$nbIterMCMC
  if(is.null(paramsMCMC$opEchan)){
    opEchan <- list(connectParam = TRUE, blockProp = TRUE, ZRow = TRUE, ZCol = TRUE)
  }else{
    opEchan <- paramsMCMC$opEchan
  }
  
  
 
  # save
  if (opSave){
    seqConnectParam = array(0,c(KRow,KCol,1+B))
    seqConnectParam[,,1] = H.mc$connectParam
    seqZ  = lapply(1:M,function(m){list(row=array(0,c(nbNodes[m,1],KRow,1+B)),col=array(0,c(nbNodes[m,2],KCol,1+B)))})
    for (m in 1:M){
      seqZ[[m]]$row[,,1] = H.mc$Z[[m]]$row
      seqZ[[m]]$col[,,1] = H.mc$Z[[m]]$col
    }
    if(model=='iidColBipartiteSBM'){
      seqBlockProp = list(row =  matrix(0,KRow,1+B), col =  matrix(0,KCol,1+B))
      seqBlockProp$row[,1] = H.mc$blockProp$row
      seqBlockProp$col[,1] = H.mc$blockProp$col
    }
    if(model=='piColBipartiteSBM'){
      seqBlockProp = list(row =  array(0,c(M,KRow,1+ B)), col =  array(0,c(M,KCol,1+B)))
      seqBlockProp$row[,,1] = H.mc$blockProp$row
      seqBlockProp$col[,,1] = H.mc$blockProp$col
    }
  }
    
    
  
  for (iterMCMC in 2:(1+B)){
  
    ################################################
    ############# simulation of connectParam
    ################################################
    
    if (opEchan$connectParam == 1) {
      
      S  <- lapply(1:M,function(m){t(H.mc$Z[[m]]$row) %*% matrix(1,nbNodes[m,1],nbNodes[m,2]) %*% H.mc$Z[[m]]$col})
      SY  <- lapply(1:M,function(m){t(H.mc$Z[[m]]$row) %*% collecNetworks[[m]] %*% H.mc$Z[[m]]$col})
      alpha_sim <- alpha.t*(hyperparamPrior$connectParam$alpha +  Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$alpha
      beta_sim <- alpha.t*(hyperparamPrior$connectParam$beta +  Reduce(`+`, S) - (emissionDist == 'bernoulli') * Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$beta
      if(emissionDist == 'bernoulli'){H.mc$connectParam <- matrix(rbeta(KRow*KCol,alpha_sim,beta_sim),KRow,KCol)}
      if(emissionDist == 'poisson'){H.mc$connectParam <- matrix(rgamma(KRow*KCol,alpha_sim,beta_sim),KRow,KCol)}
      
      if(opSave){seqConnectParam[,,iterMCMC] <- H.mc$connectParam}
    }
  
    ################################################
    ############# simulation of blockParam
    ################################################
    
    if (opEchan$blockProp == 1) {
      
      S.row <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]]$row)}))
      S.col <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]]$col)}))
      
      if(model=='iidColBipartiteSBM'){
        e_row <- alpha.t*(colSums(S.row) + hyperparamPrior$blockProp$row) + (1 - alpha.t) * hyperparamApproxPost$blockProp$row
        H.mc$blockProp$row <- c(rdirichlet(1,e_row))
        e_col <- alpha.t*(colSums(S.col) + hyperparamPrior$blockProp$col) + (1 - alpha.t) * hyperparamApproxPost$blockProp$col
        H.mc$blockProp$col <- c(rdirichlet(1,e_col))
  
        if(opSave){
          seqBlockProp$row[,iterMCMC] = H.mc$blockProp$row 
          seqBlockProp$col[,iterMCMC] = H.mc$blockProp$col 
        }
      }
      if(model=='piColBipartiteSBM'){
        e_row <- alpha.t*(S.row + hyperparamPrior$blockProp$row) + (1 - alpha.t) * hyperparamApproxPost$blockProp$row
        H.mc$blockProp$row <- t(sapply(1:M,function(m){rdirichlet(1,e_row[m,])}))
        e_col <- alpha.t*(S.col + hyperparamPrior$blockProp$col) + (1 - alpha.t) * hyperparamApproxPost$blockProp$col
        H.mc$blockProp$col <- t(sapply(1:M,function(m){rdirichlet(1,e_col[m,])}))
        if(opSave){
          seqBlockProp$row[,,iterMCMC] = H.mc$blockProp$row 
          seqBlockProp$col[,,iterMCMC] = H.mc$blockProp$col 
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
        if(model=='iidColBipartiteSBM'){piRow.m <- H.mc$blockProp$row}
        if(model=='piColBipartiteSBM'){piRow.m <- H.mc$blockProp$row[m,]}
        logZRow <-  matrix(log(piRow.m),nrow=nbNodes[m,1],ncol=KRow,byrow = TRUE)  
        
        if(emissionDist == 'bernoulli'){
          logYRow <- collecNetworks[[m]]%*%H.mc$Z[[m]]$col %*%t(log(H.mc$connectParam)) +   (1-collecNetworks[[m]])%*%H.mc$Z[[m]]$col %*%t(log(1-H.mc$connectParam))
        }
        if(emissionDist == 'poisson'){
          logYRow <- collecNetworks[[m]]%*%H.mc$Z[[m]]$col %*%t(log(H.mc$connectParam))  - matrix(1,nbNodes[m,1],nbNodes[m,2]) %*%  H.mc$Z[[m]]$col %*% t(H.mc$connectParam)
        }
        logProbZRow  <- alpha.t * (logYRow + logZRow) + (1-alpha.t)*log(hyperparamApproxPost$collecTau[[m]]$row )
        ProbZRow <- fromBtoTau(logProbZRow, eps = 10^-10)
        H.mc$Z[[m]]$row <- t(sapply(1:nbNodes[m,1],function(i){rmultinom(1,size=1,prob = ProbZRow[i,])}))
        if(opSave){
          seqZ[[m]]$row[,,iterMCMC] = H.mc$Z[[m]]$row
        }
      }
      
        ############ Zcol 
      if (opEchan$ZCol) {  
        if(model=='iidColBipartiteSBM'){piCol.m <- H.mc$blockProp$col}
        if(model=='piColBipartiteSBM'){piCol.m <- H.mc$blockProp$col[m,]}
        logZCol <- matrix(log(piCol.m),nrow=nbNodes[m,2],ncol=KCol,byrow = TRUE) 
        if(emissionDist == 'bernoulli'){
          logYCol <- t(collecNetworks[[m]])%*%  H.mc$Z[[m]]$row  %*%log(H.mc$connectParam)  +   t(1-collecNetworks[[m]])%*%  H.mc$Z[[m]]$row  %*%log(1-H.mc$connectParam)
        }
        if(emissionDist == 'poisson'){
          logYCol <- t(collecNetworks[[m]])%*%  H.mc$Z[[m]]$row  %*%log(H.mc$connectParam)  - matrix(1,nbNodes[m,1],nbNodes[m,2]) %*%  H.mc$Z[[m]]$col %*% t(H.mc$connectParam)
        }
        logProbZCol <-  alpha.t * (logZCol + logYCol) + (1-alpha.t)*log(hyperparamApproxPost$collecTau[[m]]$col )
        ProbZCol <- fromBtoTau(logProbZCol, eps = 10^-10)
        H.mc$Z[[m]]$col <- t(sapply(1:nbNodes[m,2],function(j){rmultinom(1,size=1,prob = ProbZCol[j,])}))
        if(opSave){
        seqZ[[m]]$col[,,iterMCMC] = H.mc$Z[[m]]$col
        }
        
      }
    }
    
      
      
  }
  res <- H.mc
  if(opSave){res = list(H.mc = H.mc,seqConnectParam = seqConnectParam,seqBlockProp = seqBlockProp, seqZ = seqZ)}
  return(res)
}
      
# cat('Fin MCM.Kernel')
