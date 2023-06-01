checkSample = function(HSample,MC,M,K,hyperparam,emissionDist,model){
  
  blockPropSample <- HSample$blockPropSample
  ZSample <- HSample$ZSample
  connectParamSample <- HSample$connectParamSample  
  collecTau <- hyperparam$collecTau
  #---------------  connectParamSample
  print('--------------------dim of connectParamSample--------------------')
  print(sum(dim(connectParamSample)==c(K,K,MC))==3)
  expected_mean <- hyperparam$connectParam$alpha / (hyperparam$connectParam$beta + (emissionDist == 'bernoulli')*hyperparam$connectParam$alpha)
  print(cbind(c(expected_mean),c(apply(connectParamSample,c(1,2),mean))))
  
  #--------------- blockPropSample
  print('--------------------dim of blockPropSample--------------------')
  if(model == 'iidColSBM'){
    print(dim(blockPropSample)==c(K,MC))
    
    print(cbind(hyperparam$blockProp/sum(hyperparam$blockProp),apply(blockPropSample,1,mean)))

    
  }
  if(model == 'piColSBM'){
    print(dim(blockPropSample)==c(M,K,MC))
    print(cbind(c(normByRow(hyperparam$blockProp)),c(apply(blockPropSample,c(1,2),mean))))
  }
  
  
  #---------------   Z 
  print('--------------------dim of Z-----------------')
  dimZ <- sapply(ZSample,function(Z.m){dim(Z.m)})

  print(sum(dimZ[1,]==nbNodes)==M)
  
  print(sum(dimZ[2,]==K)==M)

  print(sum(dimZ[3,]==MC)==M)

   
  ######################### 
  mrand <- sample(1:M,1)
  um <- apply(ZSample[[mrand]],c(1,2),mean)
  i <- sample(1:nbNodes[mrand],1)
  print(um[i,])
  if(is.null(collecTau)){
    if(model == 'piColSBM'){print(normByRow(hyperparam$blockProp)[mrand,])}
    if(model == 'iidColSBM'){print(hyperparam$blockProp/sum(hyperparam$blockProp))}
  }
  
  if(!is.null(collecTau)){
    print(collecTau[[mrand]][i,])
  }
  
  
  
  
  
}
