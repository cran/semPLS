#' @include resempls.R
#' roxygen()

#' Bootstrapped standard errors and confidence intervals for sempls objects.
#' The code is adapted from package 'sem' by J. Fox
#'
#'
##' @param sempls an object of class 'sempls'
##' @param nboot the number of bootstrap samples to be drawn
##' @param start a character indicating whether to restart the algorithm 'ones' or to reuse the outer weights of sempls 'old'
##' @param bootMethod a character indicating which method to use for the bootstrap, see details
##' @return object of class bootsempls inhariting from boot
##' @export
##' @callGraph
##' @seealso \code{\link{sempls}}
##' @examples
##' data(mobi)

bootsempls <- function(object, nboot=200, start=c("ones", "old"),
                method=c("ConstructLevelChanges", "IndividualSignChanges", "Standard"),
                verbose=TRUE, ...){
    method <- match.arg(method)
    if(method=="IndividualSignChanges"){
      stop("Not yet implemented.\n",
           "Try 'ConstructLevelChanges' or 'Standard'.")
    }
    refit <- function(){
        data <- data[indices,]
        refitted_model <- resempls(object, data, start, method)
        refitted_model
    }
    if (!require("boot")) stop("package boot not available")
    # the following 2 lines borrowed from boot in package boot
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    warn <- options(warn=-2)
    on.exit(options(warn)) # insure restore even in event of error
    nErrors <- 0
    data <- object$data
    N <- nrow(data)
    coefficients <- object$coefficients[,2]
    coef_names <- rownames(object$coefficients)
    coefs <- matrix(numeric(0), nrow=nboot, ncol=length(coefficients))
    attr(coefs, "path") <- object$coefficients[,1]
    colnames(coefs) <- coef_names
    clcIndex <- NULL
    tryErrorIndices <- NULL
    if (verbose) cat("Resample: ")
    for (b in 1:nboot){
      if (verbose){
        treil <- paste(rep(" ", floor(log10(nboot)) - floor(log10(b))), collapse="")
        ndel <- paste(rep("\b", floor(log10(nboot)) + 1), collapse="")
        if(b==1) cat(paste(treil, b, sep=""))
        if(b!=nboot) cat(paste(ndel, treil, b, sep=""))
        else cat(paste(ndel, b, " Done.\n", sep=""))
      }
      for (try in 1:11){
        if (try > 10) stop("more than 10 consecutive convergence failures")
        indices <- sample(N, N, replace=TRUE)
        res <- try(refit(), silent=TRUE)
        if (inherits(res, "try-error")){
          nErrors <- nErrors + 1
          tryErrorIndices <-  rbind(tryErrorIndices, indices)
        }
        else {
          # for construct level changes
          if(method=="ConstructLevelChanges"){
            clcIndex[[b]] <- res$clcIndex
          }
          # Standard
          coefs[b,] <- res$coefficients[,2]
          break()
        }
      }
    }
    options(warn)
    if (nErrors > 0) warning("There were ", nErrors,
                             " apparent convergence failures;\n",
                             "  these are discarded from the ",
                              nboot, " bootstrap replications returned.")
    res <- list(t0=coefficients, t=coefs, nboot=nboot, data=data, seed=seed,
                statistic=refit, sim="ordinary", stype="i", call=match.call(),
                tryErrorIndices=tryErrorIndices, clcIndex=clcIndex)
    class(res) <- c("bootsempls", "boot")
    res
}

##' @param x object of class bootsempls
##' @param digits number of digits to print

##' @return Prints the estimate, the bias and the standard error.
print.bootsempls <- function(x, digits = getOption("digits"), ...){
    t <- x$t
    t0 <- x$t0
    result <- data.frame("Estimate"=t0, "Bias"=colMeans(t) - t0,
        "Std.Error"=apply(t, 2, sd))
    rownames(result) <- attr(t, "path")
    cat("Call: ")
    dput(x$call)
    cat("\n")
    print(result, digits=digits)
    invisible(x)
}

##' @param object an object of class sempls
##' @param type the type of confidence interval to be computed
##' @param level the level of the confidence interval

##' @return Calculates confidence intervals for the estimates.
summary.bootsempls <- function(object,
    type=c("perc", "bca", "norm", "basic", "none"), level=0.95, ...){
    if ((!require("boot")) && (type != "none")) stop("boot package unavailable")
    type <- match.arg(type)
    t <- object$t
    t0 <- object$t0
    object$R <- object$nboot
    result <- data.frame("Estimate"=t0, "Bias"=colMeans(t) - t0,
        "Std.Error"=apply(t, 2, sd))
    if (type != "none"){
        p <- length(t0)
        lower <- upper <- rep(0, p)
        low <- if (type == "norm") 2 else 4
        up  <- if (type == "norm") 3 else 5
        noCi <- NULL
        for (i in 1:p){
          if (boot:::const(t[,i], min(1e-08, mean(t[,i], na.rm = TRUE)/1e+06))){
            lower[i] <- upper[i] <- NA
            noCi <- append(noCi, i)
          }
          else{
            ci <- as.vector(boot.ci(object, type=type, index=i,
                conf=level)[[type, exact=FALSE]])
            lower[i] <- ci[low]
            upper[i] <- ci[up]
            }
        }
        result$Lower <- lower
        result$Upper <- upper
        }
    rownames(result) <- colnames(t)
    attr(result, "path") <- attr(t, "path")
    result <- list(table=result, call=object$call, level=level, type=type)
    class(result) <- "summary.bootsempls"
    result
}

##' @param x object of summaryBootsempls
##' @param digits number ot digits to print

##' @return Prints the estimate, the bias, the standard error and confidence interval
print.summary.bootsempls <- function(x, digits = getOption("digits"), ...){
    cat("Call: ")
    dput(x$call)
    cat("\n")
    if (x$type != "none") {
        cat(paste("Lower and upper limits are for the", 100*x$level,
            "percent", x$type, "confidence interval\n\n"))
        }
    print(x$table, digits=digits)
    invisible(return(x))
}

