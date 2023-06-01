library(parallel)

#-----------------------------------------------------------------
# Computes ESS from unnormalized log-weights
#-----------------------------------------------------------------


ComputecESS <- function(W.t.minus.1,w.increm.t) {
  cESS <- sum(W.t.minus.1  *  w.increm.t)^2 / sum(W.t.minus.1  *  w.increm.t^2)
  return(cESS)
}

#-----------------------------------------------------------------
# Dichotomie to find next alpha
#-----------------------------------------------------------------

FindAlpha.cESS <- function(W.t.minus.1, log.rho, cESS.rate, alpha.t.minus.1, tol=1e-4, op.save=FALSE) {
  
  cESS_vec <- c()
  alpha_vec <- c()
  
  threshold <- cESS.rate
  alpha.new = 1;
  cESS <- ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.new - alpha.t.minus.1) * (log.rho - mean(log.rho))))
  if ( (cESS < threshold)|is.na(cESS)) {
    alpha.left <- alpha.t.minus.1;
    cESS.left <- ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.left - alpha.t.minus.1)  *  (log.rho - mean(log.rho))))
    alpha.right <- 1;
    cESS.right <-  ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.right - alpha.t.minus.1) * (log.rho - mean(log.rho))))
    
    while (is.na(cESS.right) & (alpha.right > alpha.t.minus.1)) {
      alpha.right <- alpha.right  *  0.99;
      cESS.right <-  ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.right - alpha.t.minus.1)  *  (log.rho - mean(log.rho))))
    }
    
    alpha.new <- (alpha.left + alpha.right)/2; cESS <- ComputecESS(W.t.minus.1,w.increm.t=exp((alpha.new-alpha.t.minus.1) * (log.rho-mean(log.rho))))
    if (op.save == TRUE) {cESS_vec <- c(cESS_vec,cESS); alpha_vec <- c(alpha_vec,alpha.new)}
    diff <- 2 * tol; niter = 0
    while ((diff > tol) & (alpha.new > alpha.t.minus.1 + 1e-4) & (niter < 1e3)) {
      niter = niter + 1
      if (cESS > threshold) {alpha.left <- alpha.new; cESS.left <- cESS}else{alpha.right <- alpha.new; cESS.right <- cESS}
      alpha.new <- (alpha.left + alpha.right)/2; cESS <- ComputecESS(W.t.minus.1,w.increm.t=exp((alpha.new-alpha.t.minus.1) * (log.rho-mean(log.rho))))
      diff <- abs(cESS - threshold)
      # cat(alpha.new, cESS, diff, '\n')
      if (op.save == TRUE) {cESS_vec <- c(cESS_vec,cESS);alpha_vec<- c(alpha_vec,alpha.new)}
    }
    
  }
  
  output <-  list(alpha.new=alpha.new,cESS_vec=cESS_vec,alpha_vec=alpha_vec)
  return(output)
}

#-----------------------------------------------------------------
# Objects to sample
#-----------------------------------------------------------------

Func.resampling <- function(HSample,Resample,model) {
  
  HSample$connectParamSample = HSample$connectParamSample[,,Resample]
  for (m in 1:M){
    HSample$ZSample[[m]]$row <-   HSample$ZSample[[m]]$row[,,Resample]
    HSample$ZSample[[m]]$col <-   HSample$ZSample[[m]]$col[,,Resample]
  }
  if(model=='iidColBipartiteSBM'){
    HSample$blockPropSample$row <- HSample$blockPropSample$row[,Resample]
    HSample$blockPropSample$col <- HSample$blockPropSample$col[,Resample]
  }
  if(model=='piColBipartiteSBM'){
    HSample$blockPropSample$row <- HSample$blockPropSample$row[,,Resample]
    HSample$blockPropSample$col <- HSample$blockPropSample$col[,,Resample]
  }
  return(HSample)
}



#-----------------------------------------------------------------
# Tools functions 
#-----------------------------------------------------------------

Func.search <- function(HSample,mc,model) {
  H.mc <- list(connectParam = HSample$connectParamSample[,,mc])
  H.mc$Z <- lapply(1:M, function(m){list(row = HSample$ZSample[[m]]$row[,,mc],col = HSample$ZSample[[m]]$col[,,mc])})
  if (model == 'iidColBipartiteSBM'){
    H.mc$blockProp <- list(row = HSample$blockPropSample$row[,mc],col = HSample$blockPropSample$col[,mc])
  }
  if (model == 'piColBipartiteSBM'){
    H.mc$blockProp <- list(row = HSample$blockPropSample$row[,,mc],col = HSample$blockPropSample$col[,,mc])
  }
  return(H.mc)
}


#------------------------------------
#   Computation of the log marginal likelihood
#-----------------------------------


estim.loglikmarg.U = function(RES_SMC) {
  U = RES_SMC$U
  niter = length(U)
  l <- sum(diff(RES_SMC$alpha.vec) * (U[1:(niter-1)] + U[2:niter]))/2
  return(l)
}

#-----------------------------------------------------------------
# SMC ############################################  SMC version 2
#-----------------------------------------------------------------

SMCColBipartiteSBM<- function(data,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC){
  
  
  paramsMCMC <- estimOptionsSMC$paramsMCMC  
  MC <- estimOptionsSMC$MC
  ESS.rate <- estimOptionsSMC$ESS.rate
  cESS.rate <- estimOptionsSMC$cESS.rate
  opSave <- estimOptionsSMC$opSave #  FALSE
  op.parallel <- estimOptionsSMC$op.parallel 
  op.print <- estimOptionsSMC$op.print #TRUE
  NB.iter.max <- estimOptionsSMC$NB.iter.max # Inf
  op.SMC.classic <- estimOptionsSMC$op.SMC.classic
  os  <- op.parallel$os
  
  #-------------------------- 
  KRow <- nrow(hyperparamPrior$connectParam$alpha)
  KCol <- ncol(hyperparamPrior$connectParam$alpha)
  M <- data$M
  nbNodes <- data$nbNodes
  # ----------------- parameters of the SMC
  
  #############################################################################
  #------------------   ITERATION 0  : simulation with the first init (VBEM or Prior)
  ##############################################################################
  alpha.t <- 0; 
  t <- 1;  
  vec.log.ratio.Z <- c(0); KL <- c(0, NA);

  if (op.SMC.classic) {
    HSample <- rParamZ(MC, hyperparamPrior, emissionDist, model,nbNodes);
    }else{
    HSample <- rParamZ(MC, hyperparamApproxPost, emissionDist, model,nbNodes);
  }
  W.t <- rep(1/MC,MC)
  RES <- list(HSample.0 = HSample,W.0 = W.t);
  
  #--------------------- 
  alpha.vec <- c(alpha.t);
  vec.log.ratio.Z <- c(0); 
  MI <- array(0,c(M,2,100))
  MI[,,1] <- myMutualInformationZ(HSample$ZSample,1:10)
  vec.resampling <- c(TRUE)

  if (opSave == TRUE) {
    RES.t <- list(); 
    RES.t[[t]] <- list(HSample = HSample,W.t = W.t,alpha.vec = alpha.vec,alpha.t = alpha.t)
  }
  
  logPrior.sample <- logPrior(HSample,M, MC, hyperparamPrior,emissionDist,model)
  if (op.SMC.classic) {
    logApproxPost.sample <- logPrior.sample
    }else{
    logApproxPost.sample <- logApproxPost(HSample,M, MC, hyperparamApproxPost,emissionDist,model)
    }
  condloglik.sample <-likelihood(data, HSample,emissionDist)
  log.rho <- logPrior.sample  +  condloglik.sample - logApproxPost.sample
  
  r.Inf <- which(log.rho == -Inf)
  if (length(r.Inf) > 0) {log.rho[r.Inf] = min(log.rho[-r.Inf])/2}
  U <- c(sum(W.t * log.rho))

  
  #-------------------------------------------------------------------------------------
  #------------------   ITERATIONS suivantes
  #-------------------------------------------------------------------------------------
  while ((alpha.t < 1)  & (t < NB.iter.max)) {
    t <- t  + 1 #(first value  = 1)
    W.t.minus.1 <-  W.t
    
    #--- calculation of the new alpha
    res.find.alpha.new <- FindAlpha.cESS(W.t.minus.1,log.rho, cESS.rate,alpha.t, tol = 1e-4,op.save = FALSE)
    alpha.new <- res.find.alpha.new$alpha.new
    delta.alpha.t <- alpha.new - alpha.t;
    if (delta.alpha.t < 0) {browser()}
    alpha.t <- alpha.new
    if (alpha.t > 1) {alpha.t <- 1};
    alpha.vec <- c(alpha.vec,alpha.t);
    
    #---- Compute the weigths of the new particules. ----------------------------------------------------
    Anum = delta.alpha.t * mean(log.rho) ## nomalizing constant to prevent numerical problems
    w.increm.t <- exp(delta.alpha.t  *  log.rho - Anum)
    Ww.t <-  W.t.minus.1 * w.increm.t
    W.t <- Ww.t/sum(Ww.t)
    
    #---  compute the log ratio of the marginal likelihoods---------------------------------------
    #--- ratio.Z <- sum(W.t.minus.1 * w.increm.t) # Sum (W^(n-1)  * wtilde_n ). Cf calculs et Del Moral
    log.ratio.Z <- log(sum(W.t.minus.1 * w.increm.t)) + Anum
    vec.log.ratio.Z <- c(vec.log.ratio.Z,log.ratio.Z)
    
    #------ Resampling if necesary  -----------------------------------------------------------------
    ESS.t <- 1/sum((W.t)^2)
    if (is.nan(ESS.t)) {print("Error in ESS.t. Break"); t = t - 1; save.image(file = "fail.Rdata"); break}
    resampling.step <- (ESS.t < MC * ESS.rate)
    vec.resampling <- c(vec.resampling,resampling.step)
    
    if (resampling.step == TRUE) {
      Resample <- sample(1:MC, replace = T, prob = W.t)
      W.t <- rep(1/MC,MC)
      HSample <- Func.resampling(HSample,Resample,model)
    }
    
    
    #--------  Propagate the particles with pi_t as invariant kernel -------------
    
  
    f_MCMC <- function(mc) {
        H.mc <- Func.search(HSample,mc,model) 
        GMMB.m <- MCMCKernel(data, H.mc.init = H.mc, alpha.t, hyperparamPrior,hyperparamApproxPost, emissionDist, model,paramsMCMC, opSave=FALSE)
        return(GMMB.m)
        }
  
    if (os == "unix") {
      OUTPUT_MCMC <- mclapply(1:MC, f_MCMC, mc.preschedule = TRUE, mc.cores = op.parallel$mc.cores)
    }else{
      OUTPUT_MCMC <- lapply(1:MC,f_MCMC)
    }
    # cat('Fin MCMC \n')
    HNewSample <- list()
    HNewSample$connectParamSample <- vapply(1:MC,function(mc){OUTPUT_MCMC[[mc]]$connectParam},matrix(0,KRow,KCol))
    HNewSample$ZSample <- vector("list", length = M)
    for(m in 1:M){
      HNewSample$ZSample[[m]] <- list(row='',col='')
      HNewSample$ZSample[[m]]$row <- vapply(1:MC,function(mc){OUTPUT_MCMC[[mc]]$Z[[m]]$row},matrix(0,nbNodes[m,1],KRow))
      HNewSample$ZSample[[m]]$col <- vapply(1:MC,function(mc){OUTPUT_MCMC[[mc]]$Z[[m]]$col},matrix(0,nbNodes[m,2],KCol))
    }
    if(model=='iidColBipartiteSBM'){
      HNewSample$blockPropSample$row <- vapply(1:MC,function(mc){OUTPUT_MCMC[[mc]]$blockProp$row},rep(0,KRow))
      HNewSample$blockPropSample$col <- vapply(1:MC,function(mc){OUTPUT_MCMC[[mc]]$blockProp$col},rep(0,KCol))
    }
    if(model=='piColBipartiteSBM'){
      HNewSample$blockPropSample$row <- vapply(1:MC,function(mc){OUTPUT_MCMC[[mc]]$blockProp$row},matrix(0,M,KRow))
      HNewSample$blockPropSample$col <- vapply(1:MC,function(mc){OUTPUT_MCMC[[mc]]$blockProp$col},matrix(0,M,KCol))
    }
    HSample <- HNewSample
    if (opSave == TRUE) {
      RES.t[[t]] <- list(HSample = HSample,W.t = W.t,alpha.vec = alpha.vec,alpha.t = alpha.t);
      RES.t[[t]]$w.increm.t <- w.increm.t;
      save(RES.t,file = 'encours.Rdata');
    }
    
    MI[,,t] <- myMutualInformationZ(HSample$ZSample,1:10)
    
    #------ Computation avec the log.rho to compute the w.increment.t  and U for marginal likelihood --------------------------------
    
    logPrior.sample <- logPrior(HSample,M,MC, hyperparamPrior,emissionDist,model)
    if (op.SMC.classic) {
      logApproxPost.sample <- logPrior.sample
    }else{
      logApproxPost.sample <-  logApproxPost(HSample,M,MC, hyperparamApproxPost,emissionDist,model)
    }
    condloglik.sample <- likelihood(data, HSample,emissionDist)
    log.rho <- logPrior.sample  +  condloglik.sample - logApproxPost.sample
    log.rho[which(log.rho == -Inf)] = -1e4
    log.rho[which(logApproxPost.sample == -Inf)] = -1e4
    log.rho[which(log.rho == Inf)] = -1e4
    
    U.t = sum(W.t * log.rho)
    U <- c(U,U.t)
    
    # cat('Fin log.rho \n')
    
    
    #----- Computation of the various KL   +  ussful constants
    KLalpha <- c(alpha.t * (sum(W.t * log.rho)) - sum(vec.log.ratio.Z), -(1 - alpha.t) * sum(W.t * log.rho) - sum(vec.log.ratio.Z))
    KL <- rbind(KL, KLalpha)
    
   
    # cat('Fin moments \n')
    
    if (op.print == TRUE) {print(c(t,alpha.t,ESS.t,sum(vec.log.ratio.Z),MI[,,t]))}
    #if (op.print == TRUE) {print(c(t,alpha,ESS.t,U.t))}
  }### ----------------------------------------------------------- end  of algorithm
  
  
  
  log.Zalpha <- cumsum(vec.log.ratio.Z)
  KL[, 1] <- KL[, 1]
  log.Z <- log.Zalpha[length(log.Zalpha)]
  KL[, 2] <- KL[, 2]  +  log.Z

  
  ##########
  if (opSave == TRUE) {RES$RES.t <- RES.t}
  RES$HSample_end <- HSample
  RES$W.end <- W.t
  RES$alpha.vec <- alpha.vec
  RES$vec.log.ratio.Z <- vec.log.ratio.Z;
  #RES$vec.ratio.Z.2 <- vec.ratio.Z.2;
  RES$KL <- KL
  RES$vec.resampling <- vec.resampling
  RES$U <- U
  RES$MI <- MI[,,1:t]
  
  return(RES)
}

#------------------------------------------- 
#------------------------------------------------------- 
#--------------------------------------------------- 
transfZsampleIntoMatrix = function(ZSample){

  M <- length(ZSample)
  MC <- dim(ZSample[[1]]$row)[3]
  KRow <- dim(ZSample[[1]]$row)[2]
  KCol <- dim(ZSample[[1]]$col)[2]
  
  listMatZ <- list() 
    
  for (m in 1:M){
    
    ZSample.m <- ZSample[[m]]
    nrow.m <- nrow(ZSample[[m]]$row)
    ncol.m <- nrow(ZSample[[m]]$col)
    
    matZ.m.row <- matrix(0,nrow.m,MC)
    matZ.m.col <- matrix(0,ncol.m,MC)
    for (mc in 1:MC){
      matZ.m.row[,mc] <- c(ZSample.m$row[,,mc]%*% matrix(1:KRow,ncol=1))
      matZ.m.col[,mc] <-c(ZSample.m$col[,,mc]%*% matrix(1:KCol,ncol=1))
    }
    listMatZ[[m]] <- list(row = matZ.m.row, col = matZ.m.col)
  }
  return(listMatZ)
}


myMutualInformationZ  = function(ZSample,selecNodes = NULL){
  
  listMatZ <- transfZsampleIntoMatrix(ZSample)
  M <- length(listMatZ)
  MC <- dim(listMatZ[[1]]$row)[2]
  MI <- matrix(0,M,2)
  for (m in 1:M){
    
    #les Z 
    matZ.m.row <- listMatZ[[m]]$row 
    matZ.m.col <- listMatZ[[m]]$col
    if(!is.null(selecNodes)){
      matZ.m.row <- matZ.m.row[selecNodes,]
      matZ.m.col <- matZ.m.col[selecNodes,]
    }
    # row 
    MI[m,1]  <- mutualInformation(matZ.m.row)
    MI[m,2]  <- mutualInformation(matZ.m.col)
  }
  return(MI)
}

mutualInformation = function(matZ){
  
  ########## matrix of size n row (nb of variables) and MC columns (realisations)  
  MC <- ncol(matZ)
  n <- nrow(matZ)
  H.ind <- sum( sapply(1:n,function(i){freq_i = table(matZ[i,])/MC; return(sum(-freq_i*log(freq_i)))}))
  U <- plyr::count(t(matZ))$freq/MC
  H.conj <- -sum(U * log(U ))
  H.ind - H.conj
}


