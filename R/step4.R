# Outer estimation of the factor scores
step4 <-
function(data, outerW, blocks, pairwise){
  if(pairwise){
    Latent <- NULL           # factor scores
    latent <- names(blocks)  # just names
    for(i in 1:length(latent)){   
      mfblock <- blocks[[i]]
      Latent <- cbind(Latent, as.matrix(data[, mfblock]) %*%
                              as.matrix(outerW[mfblock, latent[i]]))
    }
    colnames(Latent) <- latent
    Latent <- scale(Latent)
  }
  else {Latent <- scale(as.matrix(data) %*% outerW)}  
  # the attributes for the scale are meaningless
  attributes(Latent)[c(3,4)] <- NULL  
  return(Latent)
}
