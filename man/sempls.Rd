\name{sempls}
\Rdversion{1.1}
\alias{sempls}
\alias{sempls.plsm}
\alias{print.sempls}
\alias{plot.sempls}
\alias{densityplot.sempls}
\alias{pathCoeff}
\alias{pathCoeff.sempls}
\alias{print.pathCoeff}
\alias{totalEffects}
\alias{totalEffects.sempls}
\alias{print.totalEffects}
\alias{plsWeights}
\alias{plsWeights.sempls}
\alias{print.plsWeights}
\alias{plsLoadings}
\alias{plsLoadings.sempls}
\alias{print.plsLoadings}


\title{
Structural Equation Modeling by Partial Least Squares
}
\description{
  \code{sempls} fits structural equation models by the patial least
  squares (PLS) method. The estimation is based on the raw data and
  requires no distributional assumptions.
}
\usage{
sempls(model, \ldots)

\method{sempls}{plsm}(model, data, maxit=20, tol=1e-7,
         scaled=TRUE, sum1=FALSE, wscheme="centroid", pairwise=FALSE,
         method=c("pearson", "kendall", "spearman"),
         convCrit=c("relative", "square"),
         verbose=TRUE, \ldots)

\method{print}{sempls}(x, digits=2, \ldots)
\method{plot}{sempls}(x, \ldots)
\method{densityplot}{sempls}(x, data, use=c("fscores", "prediction",
            "residuals"), \ldots)

pathCoeff(object, \ldots)
\method{pathCoeff}{sempls}(object, \ldots)
\method{print}{pathCoeff}(x, na.print=".", digits=2, abbreviate=FALSE, \ldots)

totalEffects(object)
\method{totalEffects}{sempls}(object)
\method{print}{totalEffects}(x, na.print=".", digits=2, abbreviate=FALSE,\ldots)

plsWeights(object)
\method{plsWeights}{sempls}(object)
\method{print}{plsWeights}(x, na.print=".", digits=2, abbreviate=FALSE, \ldots)

plsLoadings(object)
\method{plsLoadings}{sempls}(object)
\method{print}{plsLoadings}(x, type=c("discriminant", "outer", "cross"),
     cutoff=NULL, reldiff=0.2, na.print=".", digits=2, abbreviate=FALSE, \ldots)
}
\arguments{
  \item{model}{
    An object inheriting from class \code{plsm} as returned from
    \code{\link{plsm}} or \code{\link{read.splsm}}.
  }
  \item{\dots}{
    Arguments to be passed down.
  }
  \item{data}{
    A \code{data.frame} containing the observed variables
    (MVs). The storage mode for all the MVs included in the model must
    be \code{numeric}.
  }
  \item{maxit}{
    A \code{numeric} value, which determines the maximum number of
    iterations performed by the PLS algorithm. The default is \eqn{20}
    iterations.
  }
  \item{tol}{
    A \code{numeric} value, specifying the tolerance for the maximum relative
    differences in the outer weights. The default value is
    \eqn{10^{-7}}.
  }
  \item{scaled}{
    A \code{logical} value indicating, whether the observed
    data shall be scaled to zero mean and unit variance. The default is
    \code{TRUE}.
  }
  \item{sum1}{
    A \code{logical} value indicating, whether the outer
    weights foreach latent variable (LV) shall be standardized to sum up
    to one. The default is \code{FALSE}. Since the factor scores are
    scaled in each step of the PLS algorithm, changing this value to
    \code{TRUE} does not affect the results.
  }

  \item{wscheme}{
    A \code{character} naming the weighting scheme to
    use. Possible values are:
    \itemize{
      \item \code{"A"} or \code{"centroid"} for the centroid scheme, the default,
      \item \code{"B"} or \code{"factorial"}for the factorial scheme and
      \item \code{"C"}, \code{"pw"} or \code{"pathWeighting"} for the path weighting scheme.
    }
  }

  \item{pairwise}{A \code{logical} value indicating, whether
    correlations shall be calculated pairwise. If the observed data
    does not contain missing values, the results are not affected.
    The default is \code{FALSE}. For more details the R help,
    \code{?cor}, can be consulted.}
  \item{method}{A \code{character} naming the method to calculate
    the correlations. Possible values are:
     \itemize{
       \item \code{"pearson"} , the default,
       \item \code{"kendall"},
       \item \code{"spearman"}.
     }
    For more details on the method, the R help, \code{?cor}, can be
    consulted. Note, that despite of the \code{method} argument, pearson
    correlations are always used for the inner approximation (step 2).
  }
  \item{convCrit}{
    The convergence criteria to use:
    \itemize{
      \item \code{"relative"}, the default,
      \item \code{"square"}.
    }
  }
  \item{verbose}{
    Logical: If \code{FALSE} no status messages are printed.
  }
  \item{object}{
    An object of class \code{sempls}.
    }
  \item{x}{
    An object of the according class.
  }
  \item{type}{
    If the argument \code{what="loadings"}, \code{type} describes the
    loadings to be extracted -- those are:
    \itemize{
      \item \code{"discriminant"}, the default, contrasts outer against cross
      loadings to check for discrimant validity of the measurement model,
      \item \code{"outer"} for the outer loadings and
      \item \code{"cross"} for the cross loadings.
    }
  }
  \item{cutoff}{
    A numerical value at which to cutoff the loadings -- this means
    loadings smaller than the cutoff value will not be printed.
  }
  \item{reldiff}{
    The argument is only effectiv when \code{type="discriminant"}. It is
    a \code{numeric} value, specifying the relative difference between
    outer and cross loadings at which cross loadings will still be
    printed.
  }
  \item{na.print}{
    A \code{character} substituting values not to be printed.
  }
  \item{digits}{
    minimal number of _significant_ digits, see \code{\link{print.default}}.
  }
  \item{use}{
    The values for which the density plots are created. If
    \itemize{
      \item \code{"fscores"}: the factor scores are used,
      \item \code{"prediction"}: the estimated factor scores are used,
      \item \code{"residuals"}: the residuals are used.
    }
  }
  \item{abbreviate}{
    A logical indicating whether dimnames should be abbreviated. For
    Details see \code{\link{abbreviate}}. The default is \code{FALSE}.
  }
}

\value{
  \code{sempls} returns an object of class \code{sempls}, with the following elements:

  \item{coefficients}{
    A \code{data.frame} containing the estimates for
    all the arcs in the path model, those are the outer loadings for
    mode \sQuote{A} type LVs and outer weights for mode \sQuote{B} type LVs and path
    coefficients for those belonging to the structural model.
  }
  \item{path_coefficient}{
    The \code{matrix} of path coefficients.
  }
  \item{outer_loadings}{
    The \code{matrix} of outer loadings.
  }
  \item{cross_loadings}{
    The \code{matrix} of cross loadings.
  }
  \item{total_effects}{
    The \code{matrix} of total effects.
  }
  \item{inner_weights}{
    The \code{matrix} of inner weights.
  }
  \item{outer_weights}{
    The \code{matrix} of outer weights.
  }
  \item{factor_scores}{
    A \code{data.frame} containing the estimated factor scores for the
    LVs.
  }
  \item{data}{
    A \code{data.frame} containing the preprocessed obseravtions of the
    MVs.
  }
  \item{incomplete}{
    The index of the incomplete observations.
  }
  \item{...}{
    All the other values are just storing information used in
    the \code{call}.
  }
}

\author{
  Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}

\seealso{
   \code{\link{plsm}}, \code{\link{read.splsm}},
   \code{\link{rSquared}}, \code{\link{pathDiagram}},
   \code{\link{bootsempls}}, \code{\link{plsm2sem}},
   \code{\link[sem]{sem}}
}
\examples{
data(ECSImobi)
ecsi <- sempls(model=ECSImobi, data=mobi, wscheme="pathWeighting")
ecsi

## create plots
densityplot(ecsi)
densityplot(ecsi, use="prediction")
densityplot(ecsi, use="residuals")

## Values of 'sempls' objects
names(ecsi)
ecsi$outer_weights
ecsi$outer_loadings
ecsi$path_coefficients
ecsi$total_effects


### using convenience methods to sempls results
## path coefficients
pathCoeff(ecsi)

## total effects
totalEffects(ecsi)

## get loadings and check for discriminant validity
(l <- plsLoadings(ecsi))
# outer loadings
print(l, type="outer", digits=2)
# outer loadings greater than 0.5
print(l,type="outer", cutoff=0.5, digits=2)
# cross loadings greater than 0.5
print(ecsi, type="cross", cutoff=0.5, digits=2)


### R-squared
rSquared(ecsi)


### Create .dot representation of the path diagram and
### create .pdf file if graphviz is available.
\dontrun{
pathDiagram(ecsi, file="ecsiPLS1", edge.labels="both",
            output.type="graphics", digits=3, graphics.fmt = "pdf")

# include R-squared values
pathDiagram(ecsi, file="ecsiPLS2", edge.labels="both",
            output.type="graphics", digits=3, graphics.fmt = "pdf",
            rSquared=rSquared(ecsi))

# only the structural model
pathDiagram(ecsi, file="ecsiPLS3", edge.labels="both",
            output.type="graphics", digits=3, graphics.fmt = "pdf",
            rSquared=rSquared(ecsi), full=FALSE)
}
}
