# Calculates the total effects
totalEffects <- function(pathCoeff){
  ret <- pathCoeff
  step <- pathCoeff
  for (i in 2:ncol(pathCoeff)){
    step <- step %*% pathCoeff 
    ret <- step + ret
  }
  return(ret)
}
