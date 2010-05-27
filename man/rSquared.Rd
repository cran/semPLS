\name{rSquared}
\alias{rSquared}
\alias{rSquared.sempls}
\title{R-squared}
\description{
  Computes the R-squared values for \code{sempls} objects.
}
\usage{
rSquared(object, \ldots)

\method{rSquared}{sempls}(object, na.rm=FALSE, \ldots)
}
\arguments{
  \item{object}{An object of class \code{sempls}.
  }
  \item{\ldots}{Arguments to pass down to other methods
  }
  \item{na.rm}{logical. Should missing values be removed?
  }
}
\value{
  A matrix containing R-squared values for each endogeneous latent variable.
}

\author{
  Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}

\seealso{
  \code{\link{sempls}}, \code{\link{pathDiagram}}
}
