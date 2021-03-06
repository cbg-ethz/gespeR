% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/gespeR-methods.R
\docType{methods}
\name{loadValues}
\alias{loadValues}
\alias{loadValues,TargetRelations-method}
\alias{loadValues,gespeR-method}
\alias{unloadValues}
\alias{unloadValues,TargetRelations-method}
\alias{unloadValues,gespeR-method}
\alias{writeValues}
\alias{writeValues,TargetRelations-method}
\title{Methods for values of \code{\linkS4class{TargetRelations}} objects}
\usage{
loadValues(object)

\S4method{loadValues}{TargetRelations}(object)

\S4method{loadValues}{gespeR}(object)

unloadValues(object, ...)

\S4method{unloadValues}{TargetRelations}(object, writeValues = TRUE,
  overwrite = FALSE, path = NULL)

\S4method{unloadValues}{gespeR}(object, writeValues = TRUE,
  overwrite = FALSE, path = NULL)

writeValues(object, ...)

\S4method{writeValues}{TargetRelations}(object, overwrite = FALSE)
}
\arguments{
\item{object}{A \code{\linkS4class{TargetRelations}} object or \code{\linkS4class{gespeR}} object}

\item{...}{Additional arguments}

\item{writeValues}{Indicator, whether to write values}

\item{overwrite}{Indicator, wheter to overwrite values if file exists at path}

\item{path}{The path to write out values}
}
\value{
A \code{\linkS4class{TargetRelations}} object or \code{\linkS4class{gespeR}} object
}
\description{
Load, unload or write to file the values of a \code{\linkS4class{TargetRelations}} object
}
\examples{
data(stabilityfits)
\dontrun{
loadValues(stabilityfits$A)
}
}
\author{
Fabian Schmich
}

