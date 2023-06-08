
#############################################################################
# Sampling in the prior distribution
#############################################################################
rParamZ <- function(MC, hyperparam, emissionDist, model ,nbNodes){
  
  if( !(model %in% c('piColSBM','iidColSBM'))){stop('Mispecified model')}
  if( !(emissionDist %in% c('poisson','bernoulli'))){stop('Mispecified emission distribution')}
  
  K  <- nrow(hyperparam$connectParam$alpha)
  M <- length(nbNodes)
  
  collecTau <- hyperparam$collecTau
  
  
  #--------- Simul of connectParamSample = array(dim = c(K, K,MC))
  if(emissionDist == 'bernoulli'){
    connectParamSample <- vapply(1:MC,function(mc){matrix(rbeta(K*K,hyperparam$connectParam$alpha,hyperparam$connectParam$beta),K,K)},matrix(0,K,K))
  }
  if(emissionDist == 'poisson'){
    connectParamSample <- vapply(1:MC,function(mc){matrix(rgamma(K*K,hyperparam$connectParam$alpha,hyperparam$connectParam$beta),K,K)},matrix(0,K,K))
  }
  
  #--------- Simul of blockPropSample
  
  if(model == 'iidColSBM'){
    #--------- Simul of blockPropSample  = array(dim = c(K,MC))
    blockPropSample <- vapply(1:MC,function(mc){rdirichlet(1,hyperparam$blockProp)},rep(0,K))
  }
  
  if(model == 'piColSBM'){
    #--------- Simul of blockPropSample  = array(dim = c(K,M,MC))
      blockPropSample <- vapply(1:MC,function(mc){t(sapply(1:M,function(m){rdirichlet(1,hyperparam$blockProp[m,])}))},matrix(0,M,K))
    }
  
  
  # Sampling Z
  
  if(is.null(collecTau)){
   
    
    ZSample <- lapply(1:M,function(m){
      n.m <- nbNodes[m]
      Z.m <- vapply(1:MC, function(mc){
        if(model == 'piColSBM'){pi.m <- blockPropSample[m,,mc]}
        if(model == 'iidColSBM'){pi.m <- blockPropSample[,mc]}
        t(rmultinom(n.m,size=1,prob = pi.m))},
        matrix(0,n.m,K)
        )
      return(Z.m)
    })
  }else{
    ZSample <- lapply(1:M,function(m){
     
      n.m <- nbNodes[m]
      Z.m <- array(0,c(n.m,K,MC))
      for (i in 1:n.m){Z.m[i,,] = rmultinom(MC,size=1,prob = collecTau[[m]][i,])}
      return(Z.m)
    })
  }
  
  HSample <- list(connectParamSample = connectParamSample, blockPropSample = blockPropSample,  ZSample = ZSample)
  return(HSample)
}



#############################################################################
# COND log lik
#############################################################################
condLogLik = function(data, H.mc,emissionDist){
  

  if( !(emissionDist %in% c('poisson','bernoulli'))){stop('Mispecified emission distribution')}
  
  
  M <- data$M; 
  collecNetworks <- data$collecNetworks
  if(emissionDist == 'bernoulli'){
      v <- sapply(1:M,function(m){
        par.m <- H.mc$Z[[m]]%*% H.mc$connectParam %*% t(H.mc$Z[[m]])
        res.m  <- sum(dbinom(collecNetworks[[m]],1,par.m,log = TRUE))
        }
        )
  }
  if(emissionDist == 'poisson'){
    v <- sapply(1:M,function(m){
      par.m <- H.mc$Z[[m]]%*% H.mc$connectParam %*% t(H.mc$Z[[m]])
      res.m  <- sum(dpois(collecNetworks[[m]],par.m,log = TRUE))
    }
    )
  }
  res <- sum(v)
  return(res)
}

#############################################################################
#  COND log lik for a sample H.sample
#############################################################################

likelihood <- function(data, HSample,emissionDist){
  # cat('likelihood ')
  if( !(emissionDist %in% c('poisson','bernoulli'))){stop('Mispecified emission distribution')}
  
  MC <- dim(HSample$connectParamSample)[3]
 
  #------------------------------------
  condloglik.sample = vapply(1:MC, function(mc){
    H.mc <- list(connectParam = HSample$connectParamSample[,,mc])
    H.mc$Z <- lapply(HSample$ZSample,function(Zm){Zm[,,mc]})
    return(condLogLik(data,H.mc,emissionDist))
    }
  ,1)
  #------------------------------------
  return(condloglik.sample) 
}




#################################################################################
# PRIOR log density
#################################################################################


#-----------------------------------------------------------------------
logDistConnectParam <- function(HSample,MC, hyperparam,emissionDist){
  # 
  if( !(emissionDist %in% c('poisson','bernoulli'))){stop('Mispecified emission distribution')}
  
  if(emissionDist=='bernoulli'){
    res <- vapply(1:MC,function(mc){
      sum(dbeta(HSample$connectParamSample[,,mc],hyperparam$connectParam$alpha,hyperparam$connectParam$beta,log = TRUE))
    },1)
  }
  if(emissionDist=='poisson'){
    res <- vapply(1:MC,function(mc){
      sum(dgamma(HSample$connectParamSample[,,mc],hyperparam$connectParam$alpha,hyperparam$connectParam$beta,log = TRUE))
    },1)
  }
  
  return(res)
}

#-----------------------------------------------------------------------
logDirichletBlockProp <- function(HSample,MC, hyperparam,model){
 
  if( !(model %in% c('piColSBM','iidColSBM'))){stop('Mispecified model')}
  
  if (model=='piColSBM'){
    res <- vapply(1:MC,function(mc){sum(log(ddirichlet(HSample$blockPropSample[,,mc],hyperparam$blockProp)))},0)
  }
  if(model =='iidColSBM'){
    res <- vapply(1:MC,function(mc){log(ddirichlet(HSample$blockPropSample[,mc],hyperparam$blockProp))},0)
    
  }
  return(res)
}


#-----------------------------------------------------------------------
logMultinomZ <- function(HSample, M, MC, hyperparam,model){
  
  if( !(model %in% c('piColSBM','iidColSBM'))){stop('Mispecified model')}
  
 res.m  = rep(0,MC)
 for (m in 1:M){
     Zsample.m <-HSample$ZSample[[m]]
     if(is.null(hyperparam$collectTau)){
       t.m <- apply(Zsample.m,c(2,3),sum)
       if(model == 'piColSBM'){pi.m <- HSample$blockPropSample[m,, ]}
       if(model == 'iidColSBM'){pi.m <- HSample$blockPropSample}
       u.m <-  apply(t.m*log(pi.m),c(2),sum)
     }else{
       u.m <- vapply(1:MC,function(mc){sum(Zsample.m[,,mc]*log(hyperparam$collecTau[[m]]))},1)
      }
     res.m = res.m + u.m
 }
 return(res.m)
}

#-----------------------------------------------------------------------
logJointParamZ <- function(HSample,M, MC, hyperparam,emissionDist,model){
 
  if( !(model %in% c('piColSBM','iidColSBM'))){stop('Mispecified model')}
  if( !(emissionDist %in% c('poisson','bernoulli'))){stop('Mispecified emission distribution')}
  
  a <- logDistConnectParam(HSample,MC, hyperparam,emissionDist)
  b <- logDirichletBlockProp(HSample,MC, hyperparam,model)
  c <- logMultinomZ(HSample, M, MC, hyperparam,model) 
  return(a+b+c)
}

#------------------------------------------------------------------------
logPrior = logJointParamZ
logApproxPost = logJointParamZ







