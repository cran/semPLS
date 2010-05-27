plsm <- function(data, strucmod=NULL, measuremod=NULL, order=c("generic", "alphabetical")){
  if(is.null(strucmod)){
    cat("Choose a .csv file for the structural model!\n")
    strucmod <- as.matrix(read.csv(file.choose()))
  }
  if(is.null(measuremod)){
    cat("Choose a .csv file for the measurement model!\n")
    measuremod <- as.matrix(read.csv(file.choose()))
  }
  if(ncol(strucmod)!=2 || mode(strucmod)!="character" || class(strucmod)!="matrix")
    stop("The argument 'strucmod' must be a two column character matrix!")
  if(ncol(measuremod)!=2 || mode(measuremod)!="character" || class(measuremod)!="matrix")
    stop("The argument 'measuremod' must be a two column character matrix!")
  if(class(data)!="data.frame")
    stop("The argument 'data' must be of class 'data.frame'!")

  latent <- unique(as.vector(strucmod))
  if(any(latent %in% colnames(data)))
     stop("The latent variables are not allowed to coincide with names of observed variables!")
  manifest <- sort(setdiff(as.vector(measuremod), latent))

  if(!all(manifest %in% colnames(data)))
     stop("The manifest variables must be contained in the data.frame")

  order <- match.arg(order)
  # Adjacency matrix D for the structural model
  D <- innerW(strucmod, latent)

  # Ordering of LVs
  if (order=="generic"){
    tmp <- reorder(D)
    latent <- tmp$chain
    strucmod <- tmp$strucmod
  }
  if (order=="alphabetical"){
    latent <- sort(latent)
  }

  # Arranging the rows and columns according to the order of the LVs
  D <- D[latent, latent]

  # build blocks of manifest variables (including 'measurement mode')
  blocks <- block(latent, manifest, measuremod)

  # Ordering of MVs
  MVs <- NULL
  for(i in 1:length(blocks)) MVs <- append(MVs, blocks[[i]])

  result <- list()
  result$latent <- latent
  result$manifest <- MVs
  result$strucmod <- strucmod
  result$measuremod <- measuremod
  result$D <- D
  result$M <- initM1(model=result)
  result$blocks <- blocks
  class(result) <- "plsm"
  return(result)
}
