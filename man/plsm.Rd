\name{plsm}
\Rdversion{1.1}
\alias{plsm}
\alias{exogen}
\alias{endogen}
\alias{reflective}
\alias{formative}
\alias{indicators}
\alias{predecessors}

\title{
  Specification of Path Model
}
\description{
  Creates the specification of a path model used by \code{sempls}. The
  structural and the measurement must be specified in .csv files. The
  paths must be entered in the form of a from-to matrix. The variable in
  the first column represents the source and the second column
  represents the target of each path.
}
\usage{
plsm(data, strucmod=NULL, measuremod=NULL, order=c("generic",
  "alphabetical"))

exogen(object)
endogen(object)
reflective(object)
formative(object)
indicators(object, LV)
predecessors(object)
}

\arguments{
  \item{data}{
    A \code{data.frame} intendend to use for the fitting method, \code{sempls}
  }
  \item{strucmod}{A from-to-matrix representing the structural model. If
    \code{NULL} the path to an .csv file representing the structural
    model can be entered interactively.}
  \item{measuremod}{A from-to-matrix representing the measurement models. If
    \code{NULL} the path to an .csv file representing the the measurement
    models can be entered interactively.}
  \item{order}{
     A \code{character} describing how to order the latent variables
     (LVs). If
    \itemize{
      \item \code{"generic"} the LVs will be ordered according to their
      appearance in the causal chain of the structural model (default).
      \item \code{"alphabetical"} the LVs will be ordered alphabetically.
    }
  }

  \item{object}{An object of class \code{plsm}.}
  \item{LV}{A character value naming a latent variable contained in the
    model}
}

\value{
  The object returned is of class \code{plsm} with the elements:
\item{latent}{A \code{character} vector naming the latent variables.}
\item{manifest}{A \code{character} vector naming the manifest variables.}
\item{strucmod}{Contains only the subset of \code{path} representing the structural model.}
\item{measuremod}{Contains only the subset of \code{path} representing the measurement model.}
\item{D}{The adjacency matrix \eqn{\bm{D}} for the structural model.}
\item{M}{The adjacency matrix \eqn{\bm{M}} for the measurement model.}
\item{blocks}{A \code{list} naming the MVs belonging to each LV's block
  and telling their measurement mode.}

  \code{exogen} returns the exogenous latent variables contained in a path
  model.

  \code{endogen} returns the endogenous latent variables contained in a path
  model.

  \code{reflective} returns the reflectively measured latent variables
  contained in a path model.

  \code{formative} returns the formatively measured latent variables
  contained in a path model.

  \code{indicators} returns the index or observed variables related to
  the given latent variable.

  \code{predecessors} returns a list with the names of predecessors for
  each latent variable.
}

\author{
  Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}

\seealso{
\code{\link{sempls}}, \code{\link{read.splsm}},
}

\examples{
# getting the path to the .csv file representing the inner Model
ptf_Struc <- system.file("ECSIstrucmod.csv", package="semPLS")

# getting the path to the .csv file representing the outer Models
ptf_Meas <- system.file("ECSImeasuremod.csv", package="semPLS")

sm <- as.matrix(read.csv(ptf_Struc))
mm <- as.matrix(read.csv(ptf_Meas))
data(mobi)

ECSI <- plsm(data=mobi, strucmod=sm, measuremod=mm)
ECSI

### Interactive mode
\dontrun{
ECSI <- plsm(data=mobi)
ECSI}

exogen(ECSI)
endogen(ECSI)
reflective(ECSI)
formative(ECSI)
indicators(ECSI, "Image")
predecessors(ECSI)
}
