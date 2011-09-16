\name{plsm2sem}
\alias{plsm2sem}
\alias{plsm2sem.plsm}


\title{
  Convert Model for Use in \pkg{sem}
}
\description{
  Converts a \code{plsm} object to an object of class \code{mod} for
  usage of \code{\link[sem]{sem}} method within \pkg{sem} package.
}

\usage{
plsm2sem(model, \ldots)

\method{plsm2sem}{plsm}(model, file=stdout(), fixedVarMV=TRUE, fixedVarLV=TRUE,
          fixedLoad=character(), \ldots)
}

\arguments{
  \item{model}{An object of class \code{splsm} as returned by the method
    \code{read.splsm}.}
  \item{\ldots}{Arguments to pass down.}
  \item{file}{A \code{character} naming the file to write to. If no file argument is specified, \code{splsm2sem} writes to \code{stdout}.}
  \item{fixedVarMV}{A \code{logical} indicating whether the variances of the MVs should be fixed to one. The default is \code{TRUE}.}
  \item{fixedVarLV}{A \code{logical} indicating whether the variances of the LVs should be fixed to one. The default is \code{TRUE}.}
  \item{fixedLoad}{A \code{character} vector naming the MVs, for which the outer loading should be fixed to one.}
}


\value{
  If the \code{\link[sem]{sem}} package is available, an object of class \code{mod} (see \code{\link[sem]{specify.model}}) is returned, else a text representation of recticular action model (RAM) is written.
}

\author{
  Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}

\seealso{\code{\link[sem]{specify.model}}}

\examples{
data(ECSImobi)
if(require(sem)){
   ECSIsem <- plsm2sem(ECSImobi)
   ECSIsem
}
}
