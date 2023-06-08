noiseSampling <- function(collecNetworks,propNoise){
  
  M <- length(collecNetworks)
  res <- list()
  if (length(propNoise) ==1){propNoise <- rep(propNoise,M) }
  for (m in 1:M){
    n.m <- nrow(collecNetworks[[m]])
    E.m <- matrix(sample(rbinom(n.m^2,1,propNoise[m])),n.m, n.m)
    diag(E.m) <- 0
    res[[m]] <- collecNetworks[[m]]*(1-E.m) + (1-collecNetworks[[m]]) * (E.m)
  }
  return(res)
}