# Converts 'plsm' object into a RAM representation for the 'sem' package.
plsm2sem <- function(model, ...){
  UseMethod("plsm2sem")
}

plsm2sem.plsm <- function(model, file=stdout(), fixedVarMV=TRUE, fixedVarLV=TRUE, fixedLoad=character(), ...){
  ## measurement model correlations (outer loadings): one way arrows "->"
  mm <- NULL
  blocks <- model$blocks
  fixed <- paste(", NA, 1\n", sep="")
  for (i in 1:length(blocks)){
    lam <- paste(", lam", i, 1:length(blocks[[i]]), ", NA\n", sep="")
    lam[blocks[[names(blocks)[i]]] %in% fixedLoad] <- fixed
    if(attr(blocks[[names(blocks)[i]]], "mode") == "A"){
      mm <- append(mm,
                   paste(names(blocks)[i], " -> ", blocks[[names(blocks)[i]]],
                         lam, sep="")
                   )
    }
    if(attr(blocks[[names(blocks)[i]]], "mode") == "B"){
      mm <- append(mm,
                   paste(blocks[[names(blocks)[i]]], " -> ", names(blocks)[i],
                         lam, sep="")
                   )
    }
  }
                     
  
  ## Not Used: correlations between manifest variables: one way arrows "->"
  
  ## structural model correlations (path coefficients): one way arrows "->"; should be free
  indx <- which(model$D!=0, arr.ind=TRUE)
  beta <- paste("beta", indx[,1], indx[,2], sep="")
  # beta <- paste("beta", 1:nrow(model$strucmod), sep="")
  sm <- paste(model$strucmod[,1], " -> ", model$strucmod[,2], ", ", beta ,", NA\n", sep="")
  
  ## variances of manifest variables: double headed arrows "<->"; should be fixed to one
  if(fixedVarMV){
    mVar <- paste(model$manifest, " <-> ", model$manifest, ", NA, 1\n", sep="" )
  }
  if(!fixedVarMV){
    the <- paste("the", 1:length(model$manifest), sep="")
    mVar <- paste(model$manifest, " <-> ", model$manifest,", ", the, ", NA\n", sep="" )
  }

  
  ## variances of latent variables: double headed arrows "<->"; should be fixed to one
  if(fixedVarLV){
    lVar <- paste(model$latent, " <-> ", model$latent, ", NA, 1\n", sep="" )
  }
  if(!fixedVarLV){
    psi <- paste("psi", 1:length(model$latent), sep="")
    lVar <- paste(model$latent, " <-> ", model$latent,", ", psi, ", NA\n", sep="" )
  }

  # print file
  if(require(sem)==FALSE) stop("Package 'XML' is required.")
  if(require(sem)==FALSE & file==1){
    cat(mm, sm, mVar, lVar, "\n", file=file)
  }
  if(require(sem)==FALSE & file!=1){
    cat(mm, sm, mVar, lVar, "\n", file=file)
    cat("Now you can run specify.model(\"", file, "\").\n",
      "See description of sem package.\n", sep="")
  }
  if(require(sem)==TRUE & file!=1){
    cat(mm, sm, mVar, lVar, "\n", file=file)
    sem_model <- specify.model(file)
    return(sem_model)
  }
  if(require(sem)==TRUE & file==1){
    tmp <- file()
    cat(mm, sm, mVar, lVar, "\n", file=tmp)
    sem_model <- specify.model(file=tmp)
    close(tmp)
    return(sem_model)
  }
}
  
