# for object of class: plsm
exogen <- function(object){
    if(!inherits(object, "plsm")) stop("Object must be of class 'plsm'!")
    ret <- names(which(colSums(object$D)==0))
    return(ret)
}

endogen <- function(object){
    if(!inherits(object, "plsm")) stop("Object must be of class 'plsm'!")
    ret <- names(which(colSums(object$D)!=0))
    return(ret)
}

formative <- function(object){
    if(!inherits(object, "plsm")) stop("Object must be of class 'plsm'!")
    ret <- names(which(lapply(object$blocks, function(x){attr(x, "mode")})=="B"))
    return(ret)
}

reflective <- function(object){
    if(class(object)!="plsm") stop("Object must be of class 'plsm'!")
    ret <- names(which(lapply(object$blocks, function(x){attr(x, "mode")})=="A"))
    return(ret)
}

indicators <- function(object, LV){
    if(!inherits(object, "plsm")) stop("Object must be of class 'plsm'!")
    if(!LV %in% object$latent) stop("The LV must be contained in the model!")
    ret <- object$blocks[[LV]]
    return(ret)
}

# used in 'pathWeighting'
predecessors <- function(object){
    if(!inherits(object, "plsm")) stop("Object must be of class 'plsm'!")
    D <- object$D
    foo <- function(x) names(which(x==1))
    pred <- apply(D, 2, foo)
    return(pred)
}
