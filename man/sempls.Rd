\name{sempls}
\Rdversion{1.1}
\alias{sempls}
\alias{sempls.plsm}
\alias{print.sempls}
\alias{densityplot.sempls}

\title{
Structural Equation Modeling by Partial Least Squares
}
\description{
  \code{sempls} fits structural equation models by the patial least
  squares (PLS) method. The estimation is based on the raw data and
  requires no distributional assumptions.
}
\usage{
sempls(model, ...)

\method{sempls}{plsm}(model, data, maxit=20, tol=1e-7, scaled=TRUE,sum1=TRUE, E="A",
    pairwise=FALSE, method=c("pearson", "kendall", "spearman"),
    convCrit=c("relative", "square"), \ldots)

\method{print}{sempls}(x, \ldots)

\method{densityplot}{sempls}(x, data, use=c("fscores", "prediction",
            "residuals"), main, sub, \ldots)
}
\arguments{
  \item{model}{An object of class \code{splsm} as returned from \code{read.splsm}.
  }
  \item{\dots}{Arguments to be passed down.
  }
  \item{data}{A \code{data.frame} containing the the observed variables
    (MVs). The storage mode for all the MVs included in the model must
    be \code{numeric}.
}
  \item{maxit}{A \code{numeric} value, which determines the maximum number of
    iterations performed by the PLS algorithm. The default is \eqn{20}
    iterations.}
  \item{tol}{A \code{numeric} value, specifying the tolerance for the maximum relative
    differences in the outer weights. The default value is \eqn{10^{-7}}.}
  \item{scaled}{A \code{logical} value indicating, whether the observed
    data shall be scaled to zero mean and unit variance. The default is
    \code{TRUE}.}
  \item{sum1}{A \code{logical} value indicating, whether the outer
    weights foreach latent variable (LV) shall be standardized to sum up
    to one. The default is \code{TRUE}. Since the factor scores are
    scaled in each step of the PLS algorithm, changing this value to
    \code{FALSE} does not affect the results.}

  \item{E}{A \code{character} naming the weighting scheme to
    use. Possible values are:
    \itemize{
      \item \code{"A"} for the centroid scheme, the default,
      \item \code{"B"} for the factorial scheme and
      \item \code{"C"} for the path weighting scheme.
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
  \item{convCrit}{The convergence criteria to use:
    \itemize{
       \item \code{"relative"} , the default,
       \item \code{"square"}.
     }
  }
  \item{x}{An object of class \code{sempls} returned by the
    \code{sempls} method.
  }
  \item{use}{The values for which the density plots are created. If
    \itemize{
      \item \code{"fscores"}: the factor scores are used,
      \item \code{"prediction"}: the estimated factor scores are used,
      \item \code{"residuals"}: the residuals are used.
    }
  }
  \item{main, sub}{Graphical parameters.}
}

\value{
  \code{sempls} returns an object of class \code{sempls}, with the following elements:

  \item{coefficients}{A \code{data.frame} containing the estimates for
    all the arcs in the path model, those are the outer loadings for
    mode \sQuote{A} type LVs and outer weights for mode \sQuote{B} type LVs and path
    coefficients for those belonging to the structural model.}
  \item{path_coefficient}{The \code{matrix} of path coefficients.}
  \item{outer_loadings}{The \code{matrix} of outer loadings.}
  \item{cross_loadings}{The \code{matrix} of cross loadings.}
  \item{total_effects}{The \code{matrix} of total effects.}
  \item{inner_weights}{The \code{matrix} of inner weights.}
  \item{outer_weights}{The \code{matrix} of outer weights.}
  \item{factor_scores}{A \code{data.frame} containing the estimated factor scores for the LVs.}
  \item{data}{A \code{data.frame} containing the preprocessed obseravtions of the MVs.}
  \item{incomplete}{The index of the incomplete observations.}
  \item{...}{All the other values are just storing information used in
    the \code{call}.}
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
ecsi <- sempls(model=ECSImobi, data=mobi, E="C")
ecsi

densityplot(ecsi)
densityplot(ecsi, use="prediction")
densityplot(ecsi, use="residuals")

names(ecsi)
ecsi$outer_weights
ecsi$outer_loadings
ecsi$path_coefficients
ecsi$total_effects

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
            rSquared=Rsquared)

# only the structural model
pathDiagram(ecsi, file="ecsiPLS3", edge.labels="both",
            output.type="graphics", digits=3, graphics.fmt = "pdf",
            rSquared=Rsquared, full=FALSE)
}
}
