\name{plsm}
\Rdversion{1.1}
\alias{plsm}
\alias{mvplot}
\alias{mvpairs}
\alias{mvplot.plsm}
\alias{mvpairs.plsm}

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
plsm(data, strucmod, measuremod, order=c("generic",
     "alphabetical"), interactive=FALSE)

mvplot(model, \ldots)
\method{mvplot}{plsm}(model, data, LVs, ask=TRUE, \ldots)

mvpairs(model, \ldots)
\method{mvpairs}{plsm}(model, data, LVs, ask=TRUE, \ldots)
}

\arguments{
  \item{data}{
    A \code{data.frame} intendend to use for the fitting method, \code{sempls}.
  }
  \item{strucmod}{
    Either a from-to-matrix representing the inner/structural model or
    the path to an .csv file representing the inner/structural
    model.
  }
  \item{measuremod}{
    Either from-to-matrix representing the outer/measurement models or
    the path to an .csv file representing the outer/measurement
    models.
  }
  \item{order}{
     A \code{character} describing how to order the latent variables
     (LVs). If
    \itemize{
      \item \code{"generic"} the LVs will be ordered according to their
      appearance in the causal chain of the structural model (default).
      \item \code{"alphabetical"} the LVs will be ordered alphabetically.
    }
  }
  \item{interactive}{
    Logical indicating whether to specify the model
    interactively using \code{edit}. The default is \code{FALSE}
  }
  \item{model}{
    An object of class \code{plsm}.
  }
  \item{LVs}{
    A \code{character} vector naming the blocks of LVs for which to
  create the plots.
  }
  \item{ask}{
    See \code{?par}
  }
  \item{...}{
    Arguments to pass down to other methods, e.g., \code{par}.
  }
}

\value{
  The object returned is of class \code{plsm} with the elements:
  \item{latent}{
    A \code{character} vector naming the latent variables.
  }
  \item{manifest}{
    A \code{character} vector naming the manifest variables.
  }
  \item{strucmod}{
    Contains only the subset of \code{path} representing the structural
    model.
  }
  \item{measuremod}{
    Contains only the subset of \code{path} representing the measurement
    model.
  }
  \item{D}{
    The adjacency matrix \eqn{\bm{D}} for the structural model.
  }
  \item{M}{
    The adjacency matrix \eqn{\bm{M}} for the measurement model.
  }
  \item{blocks}{
    A \code{list} naming the MVs belonging to each LV's block
    and telling their measurement mode.
  }
  \item{order}{
    See armguments section.
  }
}

\author{
  Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}

\seealso{
  \code{\link{sempls}}, \code{\link{read.splsm}}
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

# Adjacency matrix of the structural model
ECSI$D

# Adjacency matrix of the measurement model
ECSI$M

# return all elements
ECSI

### Interactive mode
\dontrun{
# specify model in a spreadsheets
ECSI <- plsm(data=mobi, interactive=TRUE)
ECSI
}

### explore blocks of MVs
mvplot(model=ECSI, data=mobi, LVs="Expectation")
mvpairs(model=ECSI, data=mobi, LVs="Expectation")
}
