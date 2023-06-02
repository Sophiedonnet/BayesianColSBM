initCollecTau <- function(initSimple,ref){
  
  KEstim <- initSimple[[ref]]$nbBlocks
  
  nbNodes <- t(sapply(initSimple, function(sbm){sbm$nbNodes}))

  collecTau <- lapply(1:M,function(m){
    eps <-  10^-5
    alpha_ref <- initSimple[[ref]]$connectParam$mean
    KEstim_ref <- nrow(alpha_ref)
 
    alpha_m <- initSimple[[m]]$connectParam$mean
    Perm <- permutations(KEstim_ref,nrow(alpha_m))
 
    score  = Inf
    for (i in 1:nrow(Perm)){
        score_new = sum((alpha_ref[Perm[i,],Perm[i,]] - alpha_m)^2)
        if (score_new < score){
          best_perm <- Perm[i,]
          score <- score_new
        }
    }
    tau_m <- list()
    tau_m <- matrix(eps,nbNodes[m],KEstim_ref)
    tau_m[ , best_perm] = initSimple[[m]]$probMemberships
    tau_m <- normByRow(tau_m)
    return(tau_m)}
  )
  return(collecTau)
}
