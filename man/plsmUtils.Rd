\name{plsmUtils}
\Rdversion{1.1}
\alias{plsmUtils}
\alias{plsmEdit}
\alias{plsmEdit.plsm}
\alias{addLV}
\alias{addLV.plsm}
\alias{invertLVs}
\alias{invertLVs.plsm}
\alias{removeLVs}
\alias{removeLVs.plsm}
\alias{addMVs}
\alias{addMVs.plsm}
\alias{removeMVs}
\alias{removeMVs.plsm}
\alias{addPath}
\alias{addPath.plsm}
\alias{removePath}
\alias{removePath.plsm}
\alias{exogenous}
\alias{endogenous}
\alias{reflective}
\alias{formative}
\alias{indicators}
\alias{predecessors}
\alias{successors}
\alias{connected}
\alias{acyclic}

\title{
Utility methods for \code{plsm} objects.
}
\description{
  By means of the utility methods pls model descrpitions inheriting from
  class \code{plsm} can easily be altered or investigated.% All altering methods perform sanity
  %checks for the proposed model.
}
\usage{
plsmEdit(model, \ldots)
\method{plsmEdit}{plsm}(model, data, \ldots)

addLV(model, \ldots)
\method{addLV}{plsm}(model, data, LV, mode, MVs, pred, succ, \ldots)

invertLVs(model, \ldots)
\method{invertLVs}{plsm}(model, LVs, \ldots)

removeLVs(model, \ldots)
\method{removeLVs}{plsm}(model, which, \ldots)

addMVs(model, \ldots)
\method{addMVs}{plsm}(model, data, LV, MVs, \ldots)

removeMVs(model, \ldots)
\method{removeMVs}{plsm}(model, MVs, \ldots)

addPath(model, \ldots)
\method{addPath}{plsm}(model, from, to, \ldots)

removePath(model, \ldots)
\method{removePath}{plsm}(model, from, to, \ldots)

exogenous(model)
endogenous(model)
reflective(model)
formative(model)
indicators(model, LV)
predecessors(model)
successors(model)
connected(model)
acyclic(model)
}
\arguments{
  \item{model}{
    An object inheriting from class \code{plsm} as returned from
    \code{\link{plsm}} or \code{\link{read.splsm}}.
  }
  \item{data}{
    A \code{data.frame} containing the observed variables
    (MVs). The storage mode for all the MVs included in the model must
    be \code{numeric}.
  }
  \item{LV}{
    A \code{character} value naming an LV to add or the LV to associate
    with newly created MVs.
  }
  \item{mode}{
    A \code{character} value, specifying the measurement mode of the LV
    to add:
    \itemize{
      \item \code{"A"} for a reflective LV,
      \item \code{"B"} for a formative LV.
    }
  }
  \item{MVs}{
    A \code{character} vector naming MVs, which are either to be added,
    deleted or associated with a newly created LV. Note: the names of the MVs
    must exist as \code{numeric} variables in \code{data}.
  }
  \item{pred}{
    A \code{character} vector naming the predecessors of an LV. If
    \code{pred} is \code{NULL} the LV is exogenous.
  }
  \item{succ}{
    A \code{character} vector naming the successors of an LV.
  }
  \item{LVs}{
    A \code{character} vector naming LVs to invert. Invert means to
    invert the direction of pathes of associated measurement
    models. This means changing reflective measurement to formative and
    vice versa.
  }
  \item{which}{
    A \code{character} vector naming LVs to delete.
  }
  \item{from}{
    A \code{character} vector naming predecessors of LVs
    specified in \code{to}.
  }
  \item{to}{
    A \code{character} vector naming successors of LVs
    specified in \code{from}.
  }
  \item{\dots}{
    Currently unused.
    %Arguments to be passed down.
  }
}

\value{
  All utility methods return an object of class \code{\link{plsm}}.
}

\details{
  \code{exogenous} returns the exogenous latent variables contained in a path
  model.

  \code{endogenous} returns the endogenous latent variables contained in a path
  model.

  \code{reflective} returns the reflectively measured latent variables
  contained in a path model.

  \code{formative} returns the formatively measured latent variables
  contained in a path model.

  \code{indicators} returns the index or observed variables related to
  the given latent variable.

  \code{predecessors} returns a list with the names of predecessors for
  each latent variable.

  \code{successors} returns a list with the names of successors for
  each latent variable.

  \code{connected} returns \code{TRUE} if the inner model is a connected
  graph.

  \code{acyclic} returns \code{TRUE} if the inner model is recursive,
  thus represented by an acyclic graph.
}

\seealso{
   \code{\link{plsm}}, \code{\link{read.splsm}}
}
 
\examples{
data(mobi)
data(ECSIsm)
data(ECSImm)

ECSI <- plsm(data=mobi, strucmod=ECSIsm, measuremod=ECSImm)


### Some sense free examples

# Print the block of MVs associated with "Expectation".
ECSI[["blocks"]]["Expectation"]

# Change measurement model for "Expectation" from reflective to
# formative and print its block MVs
invertLVs(model=ECSI, LVs=c("Expectation"))[["blocks"]]["Expectation"]

# Print the adjacency matrix for the inner model.
ECSI[["D"]]

# Add an additional path from "Quality" to "Loyalty" and
# print the resulting adjacency matrix.
addPath(model=ECSI, from="Quality", to="Loyalty")[["D"]]

# Remove the path previously added.
removePath(model=ECSI, from="Image", to=c("Satisfaction", "Loyalty"))[["D"]]

# Print all MVs used in the model
ECSI$manifest

# Remove some MVs and print the MVs used in the model
removeMVs(model=ECSI, MVs=c("IMAG3", "CUEX1", "PERQ7"))$manifest


### some handy functions

exogenous(ECSI)
endogenous(ECSI)
reflective(ECSI)
formative(ECSI)
indicators(ECSI, "Image")
predecessors(ECSI)
}
