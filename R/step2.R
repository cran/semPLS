# Inner estimation of factor scores
step2 <-
function(Latent, innerW, blocks, pairwise){
  if(pairwise){
    fscores <- NULL
    latent <- names(blocks)
    for(i in 1:length(latent)){    
      con <- which(innerW[,i]!=0)
      fscores <- cbind(fscores, as.matrix(Latent[,con]) %*%
                                as.matrix(innerW[con,i]))
    }
    colnames(fscores) <- latent
    Latent <- scale(fscores)
  }
  else {Latent <- scale(Latent %*% innerW)}
  # the attributes for the scale are meaningless
  attributes(Latent)[c(3,4)] <- NULL  
  return(Latent)
}

