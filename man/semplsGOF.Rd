\name{semplsGOF}
\Rdversion{1.1}
\alias{semplsGOF}
\alias{rSquared}
\alias{rSquared.sempls}
\alias{print.rSquared}
\alias{qSquared}
\alias{qSquared.sempls}
\alias{print.qSquared}
\alias{dgrho}
\alias{dgrho.sempls}
\alias{print.dgrho}
\alias{comunality}
\alias{comunality.sempls}
\alias{print.comunality}
\alias{comunality}
\alias{comunality.sempls}
\alias{print.comunality}
\alias{redundancy}
\alias{redundancy.sempls}
\alias{print.redundancy}
\alias{rSquared2}
\alias{rSquared2.sempls}
\alias{print.rSquared2}
\alias{gof}
\alias{gof.sempls}
\alias{print.gof}

\title{
Quality Indices and Goodness of fit measures for pls path models
}
\description{
  A collection of method to validate the goodness of the model. Since
  there is no well identified global optimization criterion each part of
  the model needs to be validated.
}
\usage{
rSquared(object, \ldots)
\method{rSquared}{sempls}(object, na.rm=FALSE, \dots)
\method{print}{rSquared}(x, na.print=".", digits=2, \ldots)

qSquared(object, \ldots)
\method{qSquared}{sempls}(object, d=NULL, impfun, dlines=TRUE,
         total=FALSE, \dots)
\method{print}{qSquared}(x, na.print=".", digits=2, \ldots)

dgrho(object, \ldots)
\method{dgrho}{sempls}(object, \ldots)
\method{print}{dgrho}(x, na.print=".", digits=2, \ldots)

comunality(object, \ldots)
\method{comunality}{sempls}(object, \ldots)
\method{print}{comunality}(x, na.print=".", digits=2, \ldots)

redundancy(object, \ldots)
\method{redundancy}{sempls}(object, \ldots)
\method{print}{redundancy}(x, na.print=".", digits=2, \ldots)

rSquared2(object, \ldots)
\method{rSquared2}{sempls}(object, na.rm=FALSE, \dots)
\method{print}{rSquared2}(x, na.print=".", digits=2, \ldots)

gof(object, \ldots)
\method{gof}{sempls}(object, \ldots)
\method{print}{gof}(x, na.print=".", digits=2, \ldots)
}
\arguments{
  \item{object}{
    An object of class \code{sempls}.
  }
  \item{d}{
    A \code{numeric} value for the omission distance. Thus choosing
    \eqn{d=N-1}, where \eqn{N} is the number of complete observations,
    is leaving one out cross validation. This is done when \code{d}
    takes its default value \code{NULL}.
  }
  \item{impfun}{
    An user specified function to impute missing values.
  }
  \item{dlines}{
    If \code{TRUE} the same observations are deleted for a whole block
    of MVs, else each \code{d}th , counting from top left to bottom
    right, observation is deleted.
  }
  \item{total}{
    If \code{total=TRUE} total effects are used instead of path
    coefficients to calculate the predictions.
  }
  \item{na.rm}{
    If \code{na.rm=TRUE} observation with missing values are discarded 
    before analysis.
  }
  \item{x}{
    An object of the according class.
  }
  \item{na.print}{
    A \code{character} substituting values not to be printed.
  }
  \item{digits}{
    minimal number of _significant_ digits, see \code{\link{print.default}}.
  }
  %\item{justify}{
  %  character indicating if strings should left- or
  %  right-justified or left alone, passed to \code{format}.
  %}
  \item{\dots}{
    Arguments to be passed down.
  }
}

\value{
  Most GOF methods return a column vector with the names of the variables as
  rows and the respective measure as column.
}

\references{
  Esposito Vinzi V., Trinchera L., Amato S. (2010). PLS Path Modeling: From Foundations to Recent
  Developments and Open Issues for Model Assessment and Improvement. In
  Esposito Vinzi V., Chin W.W., Henseler J., Wang H.F. (eds.), Handbook of Partial Least Squares: Concepts,
  Methods and Applications in Marketing and Related Fields, chapter 2. Springer-Verlag
  Berlin Heidelberg.
}

\author{
  Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}

\seealso{
   \code{\link{sempls}}, \code{\link{plsLoadings}}
}
\examples{
data(ECSImobi)
ecsi <- sempls(model=ECSImobi, data=mobi, E="C")

### R-squared
rSquared(ecsi)

### Q-squared with omission distance d=4
qSquared(ecsi, d=4)

### Dillon-Goldstein's rho (aka composite reliability)
dgrho(ecsi)

### Comunalities
comunality(ecsi)

### Redundancy
redundancy(ecsi)

### R-squared (normal + corrected)
rSquared2(ecsi)

### Goodness of fit
gof(ecsi)

### check for discriminant validity using loadings
l <-plsLoadings(ecsi)
print(l, type="discriminant", cutoff=0.5, reldiff=0.2)
}