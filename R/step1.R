# Initial estimation of factor scores
step1 <-
function(model, data, sum1, pairwise, ...){
  M <- if(sum1==TRUE) apply(model$M, 2, sum1) else model$M
  blocks <- model$blocks
  if(pairwise){
    Latent <- NULL
    for(i in 1:length(model$latent)){
      latent <- names(blocks)[i]     # name of the i-th LV
      mfblock <- blocks[[i]]         # names of the MVs in i-th LVs block
      Latent <- cbind(Latent, as.matrix(data[, mfblock]) %*%
                              as.matrix(M[mfblock, latent]))
    }
    colnames(Latent) <- model$latent
    Latent <- scale(Latent)
  }
  else {Latent <- scale(as.matrix(data)%*% M)}   
  # the attributes for the scale are meaningless
  attributes(Latent)[c(3,4)] <- NULL  
  return(list(latent=Latent, outerW=M))
}
