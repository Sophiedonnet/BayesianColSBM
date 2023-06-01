initCollecTau <- function(initBipartite,ref){
  
  KRowEstim <- initBipartite[[ref]]$nbBlocks[1]
  KColEstim <- initBipartite[[ref]]$nbBlocks[2]
  
  nbNodes <- t(sapply(initBipartite, function(sbm){sbm$nbNodes}))

  collecTau <- lapply(1:M,function(m){
    eps  <-  10^-5
    alpha_ref <- initBipartite[[ref]]$connectParam$mean
    KRowEstim_ref <- nrow(alpha_ref)
    KColEstim_ref <- ncol(alpha_ref)
    
    alpha_m <- initBipartite[[m]]$connectParam$mean
    Perm_row <- permutations(KRowEstim_ref,nrow(alpha_m))
    Perm_col <- permutations(KColEstim_ref,ncol(alpha_m))
    score  = Inf
    for (i in 1:nrow(Perm_row)){
      for (j in 1:nrow(Perm_col)){
        score_new = sum((alpha_ref[Perm_row[i,],Perm_col[j,]] - alpha_m)^2)
        if (score_new < score){
          best_perm_row <- Perm_row[i,]
          best_perm_col <- Perm_col[j,]
          score <- score_new
        }
      }
    }
    tau_m <- list()
    tau_m$row <- matrix(eps,nbNodes[m,1],KRowEstim_ref)
    tau_m$row[ , best_perm_row] = initBipartite[[m]]$probMemberships$row
    tau_m$row <- normByRow(tau_m$row)
    tau_m$col <- matrix(eps,nbNodes[m,2],KColEstim_ref)
    tau_m$col[ , best_perm_col] = initBipartite[[m]]$probMemberships$col
    tau_m$col <- normByRow(tau_m$col)
    return(tau_m)}
  )
  return(collecTau)
}
