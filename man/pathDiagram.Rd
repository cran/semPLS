\name{pathDiagram}
\Rdversion{1.1}
\alias{pathDiagram}
\alias{pathDiagram.sempls}

\title{
  Draw Path Diagram
}
\description{
 \code{pathDiagram} creates a description of the path diagram
  for a PLS path model object of class \code{sempls} to be processed by the
  graph-drawing program \emph{dot}; see Gansner, Koutsofios, and North (2006)
  and \url{http://www.graphviz.org/}.
}
\usage{
pathDiagram(object, \ldots)


\method{pathDiagram}{sempls}(object, file, min.rank=NULL, max.rank=NULL, same.rank=NULL,
            edge.labels=c("names", "values", "both"), size=c(8,8),
            node.font=c("Helvetica", 14), edge.font=c("Helvetica", 10),
            rank.direction=c("LR", "TB"), digits=2, output.type =
            c("graphics", "dot"), graphics.fmt = "pdf", dot.options=NULL,
            rSquared=NULL, full=TRUE, \ldots)
}

\arguments{
  \item{object}{A \code{sempls} object}
  \item{\ldots}{Arguments to be passed down to \code{pathDiagram.sempls}.}
  \item{file}{
    A file in which to write the \emph{dot} description of the
    path diagram; if not specified, the description is written to standard
    output (normally the R console).}
  \item{min.rank}{A \code{character} string listing names of variables to be assigned
    minimum rank (order) in the graph; the names should be separated by commas.}
  \item{max.rank}{A \code{character} string listing names of variables to be assigned
    maximum rank in the graph; the names should be separated by commas.}
  \item{same.rank}{A \code{character} string or vector of \code{character} strings of variables
    to be assigned equivalent rank in the graph; names in each string should be
    separated by commas.}
  \item{edge.labels}{\code{"names"} to label arrows with parameter names;
    \code{"values"} to label arrows with parameter estimates, or \code{both}}
  \item{size}{The size of the graph, in inches.}
  \item{node.font}{Font name and point-size for printing variable names.}
  \item{edge.font}{Font name and point-size for printing arrow names or values.}
  \item{rank.direction}{Draw graph left-to-right, \code{"LR"}, the default,
    or top-to-bottom, \code{"TB"}.}
  \item{digits}{Number of digits after the decimal point (default, 2) to which
    to round parameter estimates.}
  \item{output.type}{If \code{"graphics"}, the default, both a \dQuote{.dot} file and
    a graphics file will be created.}
  \item{graphics.fmt}{A graphics format recognized by the dot program;
    the default is \dQuote{pdf}; \code{graphics.fmt} is also used for
    the extension of the graphics file that is created.}
  \item{dot.options}{Options to be passed to the \emph{dot} program,
    given as a \code{character} string.}
  \item{rSquared}{A matrix returned from method \code{rSquared}. If
    supplied, the R-squared values are printed in the diagram.}
  \item{full}{If \code{FALSE} the .dot file is created only for the structural model.}
}

\value{
  \code{pathDiagram} is used for its side-effect, producing a
  graph description of the model.
}

\references{
  Gansner, E., E. Koutsofios, and S. North (2006)
  Drawing graphs with \emph{dot}.
  \url{http://www.graphviz.org/Documentation/dotguide.pdf}.
}

\author{
  Armin Monecke \email{armin.monecke@stat.uni-muenchen.de}
}


\seealso{
  \code{\link{sempls}}, \code{\link{rSquared}}
}

\examples{
### create .pdf file for the path diagram
### Note that graphviz (www.graphviz.org) must be available.
\dontrun{
pathDiagram(ecsi, file="ecsiPLS1", edge.labels="both",
            output.type="graphics", digits=3, graphics.fmt = "pdf")

# include R-squared values
rSquared <- rSquared(ecsi)
pathDiagram(ecsi, file="ecsiPLS2", edge.labels="both",
            output.type="graphics", digits=3, graphics.fmt = "pdf",
            rSquared=Rsquared)

# only the structural model
pathDiagram(ecsi, file="ecsiPLS3", edge.labels="both",
            output.type="graphics", digits=3, graphics.fmt = "pdf",
            rSquared=Rsquared, full=FALSE)
}
}
