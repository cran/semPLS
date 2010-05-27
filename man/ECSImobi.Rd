\name{ECSImobi}
\Rdversion{1.1}
\alias{ECSImobi}
\alias{ECSImm}
\alias{ECSIsm}
\docType{data}
\title{
  Measurement Instrument and the ECSI Model for the Mobile Phone Industry
}
\description{
  European customer satisfaction index (ECSI) adapted to the mobile
  phone market, see Tenenhaus et al. (2005). \code{ECSImobi} is a
  convenience function to load the \code{\link{mobi}} data togehter with the
  ECSI model, represented by the from-to-matrizes \code{ECSImm} (measurement
  model) and \code{ECSIsm} (structural model).
  as \code{plspm} object, \code{ECSI}.
}
\usage{
data(ECSImobi)
data(ECSImm)
data(ECSIsm)
}

\format{
  The from-to-matrizes represent directed edges connecting the
  variables.
  \describe{
    \item{\code{from}}{The tail of an edge}
    \item{\code{to}}{The head of an edge}
  }
}

\references{
  Tenenhaus, M., V. E. Vinzi, Y.-M. Chatelin, and C. Lauro (2005)
  PLS path modeling.
  Computational Statistics \& Data Analysis 48, 159-205.
}

\author{
Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}

\seealso{
  \code{\link{mobi}}
}
\examples{
data(ECSImobi)
ecsi <- sempls(model=ECSImobi, data=mobi, E="C")
ecsi
}

\keyword{datasets}

