\name{bootsempls}
\alias{bootsempls}
\alias{print.bootsempls}
\alias{summary.bootsempls}
\alias{print.summary.bootsempls}
\alias{densityplot.bootsempls}
\alias{parallelplot.bootsempls}
\alias{parallel.bootsempls}
\title{Bootstrap a PLS path model}
\description{
Bootstraps a PLS path model in a \code{sempls} object (as returned by
the \code{\link{sempls}} method).
}

\usage{
bootsempls(object, nboot=200, start=c("ones", "old"), method=
            c("ConstructLevelChanges", "IndividualSignChanges",
            "Standard"), verbose=TRUE, strata, \ldots)

\method{print}{bootsempls}(x, digits=3, \ldots)

\method{summary}{bootsempls}(object, type=c("perc", "bca", "norm", "basic", "none"),
        level=0.95, \ldots)

\method{print}{summary.bootsempls}(x, na.print, digits = 3, \ldots)

\method{densityplot}{bootsempls}(x, data, pattern="beta", subset=NULL, \ldots)

\method{parallelplot}{bootsempls}(x, data, pattern="beta", subset=NULL, reflinesAt,
         col=c("grey", "darkred", "darkred", "black"), lty=c("solid",
         "solid","dashed", "dotted"), \ldots)
}

\arguments{
  \item{object}{
    An object of class \code{sempls} as returned by the method
    \code{sempls}.
  }
  \item{nboot}{
    The number of bootstrap replications; the default is \code{200}.
  }
  \item{start}{
    A \code{character} value defining the initialisation of
    outer weights. If
    \itemize{
     \item \code{start="ones"}, then the outer weights for each block are initialised by ones.
     \item \code{start="old"}, the final outer weights taken from the
     \code{sempls} object are used for initialisation.
    }
  }
  \item{method}{
    A \code{character} value, which can take the values:
    \itemize{
      \item \code{"ConstructLevelChanges"} (default): The vector of
      loadings for each LV in each resample is compared to the corresponding
      vector of loadings in the original sample. The signs of the weights,
      and consequently the signs of the loadings, are reversed if the
      absolute value of the sum differences between estimated loadings
      obtained from the original sample and the estimated loadings obtained
      from the resample is bigger than the absolute value of the sum sums of
      the latter, see Tenenhaus et al (2005).
      \item \code{"IndividualSignChanges"}: not implemented yet.
      \item \code{"Standard"}: No compensation for sign changes of
      resampled statistics.
    }
  }
  \item{verbose}{
    A \code{logical} indicating, whether progress of bootstrap shall
    be displayed.
  }
  \item{x}{
    An object of class \code{bootsempls} and
    \code{summary.bootsempls} respectively.
  }
  \item{na.print}{
    A \code{character} substituting values not to be printed.
  }
  \item{digits}{
    Controls the number of digits to print.
  }
  \item{type}{
    Type of bootstrapped confidence intervals to compute; the
    default is \code{"perc"} (percentile); see
    \code{\link[boot]{boot.ci}} for details.
  }
  \item{strata}{
    An integer vector or factor specifying the strata for
    multi-sample problems. If the argument is not provided, all data is
    assumed to come from the same sample. For details, see
    \code{\link[boot]{boot.ci}}.
  }
  \item{level}{
    Level for confidence intervals; default is \code{0.95}.
  }
  \item{\ldots}{
    Arguments to be passed down to other methods.
  }
  \item{data}{
    The \code{data} is not used because the \code{bootsempls}
    object already contains the required data.
  }
  \item{pattern}{
    A regular expression passed on to
    \code{\link[base]{grep}}. It is helpful to easily identify a set of
    coefficients, e.g. \code{pattern="beta"} plots the path coefficients
    only.
  }
  \item{subset}{
    Index or character vector of coefficients to
    include. Note, that \code{subset} overrides \code{pattern}.
  }
  \item{reflinesAt}{
    A vector of values at which to plot reference lines
    into the parallel cooordinates.
  }
  \item{col}{
    Colors for bootstrap statistics, sample statistic,
    lower and upper bootstrap confidence levels and reference lines.
  }
  \item{lty}{
    Line type for  bootstrap statistics, sample statistic,
    lower and upper bootstrap confidence levels and reference lines.
  }
}

\details{
  \code{boot.sempls} implements the nonparametric bootstrap, assuming an
  independent random sample. Convergence failures in the bootstrap resamples
  are discarded (and a warning printed); \eqn{10} consecutive convergence failures
  result in an error. You can use the \code{\link[boot]{boot}} function
  in the boot package for more complex sampling schemes and additional options.
}

\value{
  \code{boot.sempls} returns an object of class \code{bootsempls}, which inherits
    from class \code{boot}, supported by the \code{boot} package. The returned
    object contains the following components:
  \item{t0}{
      The estimated parameters in the model fit to the original data
      set.
  }
  \item{t}{
    a matrix containing the bootstrapped estimates, one bootstrap
    replication per row.
  }
  \item{data}{
    The data frame containing the data to which the model was fit.
  }
  \item{seed}{
    The value of \code{.Random.seed} when \code{boot.sempls} was called.
  }
  \item{statistic}{
    The function used to produce the bootstrap replications;
    this is always the local function \code{refit} from
    \code{boot.sempls}.
  }
  \item{sim}{
    Always set to \code{"ordinary"}; see the documentation for the
    \code{\link[boot]{boot}} function.
  }
  \item{stype}{
    Always set to \code{"i"}; see the documentation for the
    \code{\link[boot]{boot}} function.
  }
  \item{call}{
    The call of the \code{boot.sem} function.
  }
  \item{tryErrorIndices}{
    Contains the indices for each resample, which
    returned \code{try-error}.
  }
  \item{clcIndices}{
    When the method \code{ConstructLevelChanges} is used, it stores the
    indices of the blocks of MVs for which a sign change occured.
  }
  \item{bootIndices}{
    A matrix containing the indices of the converged
    bootstrap samples as rows.
  }
  \item{outer_weights}{
    A martrix containing, as rows, the outer weights
    for each bootsrap sample.
  }
  \item{fitted_model}{
    The fitted sempls model returned from \code{\link{sempls}}.
  }
  \item{strata}{
    The strata used. This is the vector passed to boot, if
    it was supplied or a vector of ones if there were no strata.
  }
}

\seealso{
  \code{\link[boot]{boot}}, \code{\link[sem]{boot.sem}}
}

\references{
  Tenenhaus, M., V. E. Vinzi, Y.-M. Chatelin, and C. Lauro (2005)
  PLS path modeling. Computational Statistics & Data Analysis 48,
  159-205.
}

\examples{
\dontrun{
data(ECSImobi)
ecsi <- sempls(model=ECSImobi, data=mobi)

### Bootstrapping
set.seed(123)
ecsiBoot <- bootsempls(ecsi, nboot=200, start="ones", verbose=TRUE)
summary(ecsiBoot, type="perc", level=0.95)



## inspectation of bootstrap samples
parallelplot(ecsiBoot, subset=1:ncol(ecsiBoot$t), reflinesAt=0)

# only inspecting the path coefficients
parallelplot(ecsiBoot, pattern="beta", reflinesAt=c(0,1))
densityplot(ecsiBoot, pattern="beta")

# only inspecting the outer loadings
parallelplot(ecsiBoot, pattern="lam")

# only inspecting the outer loadings for Loyalty
parallelplot(ecsiBoot, pattern="lam7", type="perc", level=0.90,
         main="Loyalty\n 200 bootstrapped outer loadings")
}
}
