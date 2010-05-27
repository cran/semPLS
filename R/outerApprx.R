# Outer Approximation (step 3 within PLS-Algorithm)
# Calculates the new outer weights.
# uses: sum1 - function to normalize the weights to sum up to 1.
outerApprx <-
function(Latent, data, blocks, sum1=TRUE, pairwise, method){
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  nl <- ncol(Latent)                      # number of latent variables
  W <- matrix(0, ncol=nl, nrow=ncol(data))
  w <- NULL
  colnames(W) <- colnames(Latent)
  rownames(W) <- colnames(data)

  for (i in 1:nl){
    mf <- as.matrix(subset(data, select=blocks[[i]])) 
    fscores <- as.matrix(Latent[,colnames(Latent)[i]])
    ## Mode A: reflective
    if (attr(blocks[[i]], "mode") == "A") {
      W[blocks[[i]],i] <- cor(fscores, mf, use, method)
    } 
    ## Mode B: formative
    if (attr(blocks[[i]], "mode") == "B") {
      W[blocks[[i]],i] <- solve(cor(mf, mf, use, method)) %*%
                          cor(mf, fscores, use, method)
    }                              
  }

  ## Normalize weights to colwise sum up to 1?
  if(sum1==TRUE){
     W <- apply(W, 2, sum1)
  }
  return(W)
}
