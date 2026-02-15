\name{print.structest}
\alias{print.structest}
\alias{summary.structest}
\title{Print and Summary Methods for structest Objects}
\usage{
\method{print}{structest}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{summary}{structest}(object, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
  \item{x, object}{an object of class \code{"structest"}.}
  \item{digits}{number of significant digits for printing.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  Print and summary methods for objects returned by \code{\link{test_t0}}
  and \code{\link{test_t1}}.
}
