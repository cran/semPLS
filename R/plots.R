# lattice:::densityplot
densityplot.sempls <- function(x, data, use=c("fscores", "prediction", "residuals"),
                               main, sub, ...){
    use <- match.arg(use)
    if(use=="fscores")         val <- x$factor_scores
    else if(use=="prediction") val <- predict(x)
    else if(use=="residuals")  val <- residuals(x)

    Y <- data.frame(NULL)
    exogenous <- exogen(x$model)
    for(i in x$model$latent){
        if(i %in% exogenous & use!="fscores") next
        tmp <- data.frame(value=val[,i], name=i)
        Y <- rbind(Y, tmp)
    }
    if(missing(main)) main=paste(deparse(substitute(x)), "\n", ifelse(use=="fscores", "factor scores", use))
    if(missing(sub)){
        sub=paste("Exogenous LVs:\n", paste(exogenous, collapse=", "))
    }
    densityplot(~value|name, data=Y, main=main, sub=sub, as.table=TRUE, ...)
 }

densityplot.bootsempls <- function(x, data, pattern="beta", subset=NULL, ...){
    ind <- grep(pattern, colnames(x$t))
    ifelse(is.null(subset),
           params <- colnames(x$t)[ind],
           ifelse(is.character(subset), params <- subset, params <- colnames(x$t)[subset])
           )
    Y <- data.frame(NULL)
    for(i in params){
        if(round(var(x$t[,i], na.rm=TRUE), digits=4)==0) next
        tmp <- data.frame(value=x$t[,i], name=i)
        Y <- rbind(Y, tmp)
    }
    densityplot(~value|name, data=Y, as.table=TRUE, ...)
 }

# lattice:::parallel
parallel.bootsempls <- function(x, data, pattern="beta", subset=NULL, reflinesAt,
                                col=c("grey", "darkred", "darkred", "black"),
                                lty=c("solid", "solid", "dashed", "dotted"), ...){
    ifelse(is.null(subset), ind <- grep(pattern, colnames(x$t)), ind <- subset)
    lower <- summary(x, ...)$table$Lower
    upper <- summary(x, ...)$table$Upper
    Y <- rbind(x$t, x$t0, lower, upper, deparse.level=0)
    if(!missing(reflinesAt)){
        Y <- rbind(Y, matrix(rep(reflinesAt, each=ncol(x$t)), nrow=length(reflinesAt), byrow=TRUE))
        origin <- c(rep("1resample", x$nboot), "2sample", "3ci", "3ci",
                    rep("4reflines", times=length(reflinesAt)))
        Y <- data.frame(Y, origin)
    }
    else Y <- data.frame(Y, origin=c(rep("1resample", x$nboot), "2sample", "3ci", "3ci"))
    parallel(~Y[ind], data=Y, groups=origin, common.scale=TRUE, col=col, lty=lty, ...)
 }
