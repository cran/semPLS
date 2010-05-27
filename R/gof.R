

# Dillon-Goldstein's rho (Composite Reliability in SmartPLS)
dgrho <- function(object){
    dgr <- matrix(NA, nrow=length(object$model$latent), ncol=2)
    rownames(dgr) <- object$model$latent
    colnames(dgr) <- c("Dillon-Goldstein's rho", "reflective MVs")
    for(i in object$model$latent){
        if(attr(object$model$blocks[[i]], "mode")=="B"){
            next
        }
        x <- object$outer_loadings[, i]
        ind <- which(x!=0)
        if(length(ind)==1){
            dgr[i,2] <- 1
            next
        }
        else {
            x <- x[ind]
            dgr[i,1] <- sum(x)^2 / (sum(x)^2 + sum(1-x^2))
            dgr[i,2] <- length(ind)
        }
    }
    return(dgr)
}

comunality <- function(object){
    com <- matrix(NA, nrow=length(object$model$latent), ncol=2)
    rownames(com) <- object$model$latent
    colnames(com) <- c("comunality", "reflective MVs")
    for(i in object$model$latent){
        if(attr(object$model$blocks[[i]], "mode")=="B"){
            next
        }
        x <- object$outer_loadings[, i]
        ind <- which(x!=0)
        if(length(ind)==1){
            com[i,2] <- 1
            next
        }
        else {
            x <- x[ind]
            com[i,1] <- 1/length(x)*sum(x^2)
            com[i,2] <- length(ind)
        }
    }
    return(com)
}

print.comunality <- function(object){
    print(object, digits=3)
    aveCom <- sum(object[,2], na.rm=TRUE)^-1 * sum(object[,1] * object[,2], na.rm=TRUE)
    paste("Average comunality:", round(aveCom, digits=3))
}



# Redundancy Example:
redundancy <- function(object){
    red <- as.matrix(comunality(object)[,1] * rSquared(object)[,1])
    colnames(red) <- "redundancy"
    return(red)
}

print.redundancy <- function(object){
    print(object, digits=3)
    aveRed <- nrow(object)^-1 * sum(object[,1], na.rm=TRUE)
    paste("Average redundancy:", round(aveRed, digits=3))
}

predict <- function(object){
    Y_hat <- object$factor_scores %*% object$path_coefficients
    return(Y_hat)
}

residuals <- function(object){
    res <- object$factor_scores - predict(object)
    return(res)
}


rSquared <- function(object, na.rm=FALSE, ...){
  Y_hat <- predict(object)
  if(sum(is.na(Y_hat)) > 0 & !na.rm) stop("Use argument 'na.rm=TRUE'!")
  R_squared <- apply(Y_hat, 2, var, na.rm=na.rm) / apply(object$factor_scores, 2, var, na.rm=na.rm)
  R_squared[R_squared==0] <- 0
  R_squared <- as.matrix(R_squared)
  R_squared <- cbind(R_squared, colSums(object$model$D))
  colnames(R_squared) <- c("R-squared", "predecessors")
  return(R_squared)
}

print.rSquared <- function(object){
    print(object, digits=3)
    aveRsquared <- nrow(object)^-1 * sum(object[,1], na.rm=TRUE)
    paste("Average R-squared:", round(aveRsquared, digits=3))
}

gof <- function(object){
    rSq <- rSquared(object)
    aveRsq <- nrow(rSq)^-1 * sum(rSq[,1], na.rm=TRUE)
    com <- comunality(object)
    aveCom <- sum(com[,2], na.rm=TRUE)^-1 * sum(com[,1] * com[,2], na.rm=TRUE)
    gof <- sqrt(aveCom * aveRsq)
    return(gof)
}
