\name{test_t0}
\alias{test_t0}
\title{T0 Test: Reliability-Dependent Test of Structural Interpretation}
\usage{
test_t0(X, z, na.rm = TRUE, verbose = FALSE)
}
\arguments{
  \item{X}{numeric matrix (n x d) of indicator variables, d >= 3.}
  \item{z}{numeric vector of length n encoding the auxiliary variable. Must have at least 2 distinct levels.}
  \item{na.rm}{logical; if \code{TRUE}, rows with any \code{NA} are removed.}
  \item{verbose}{logical; if \code{TRUE}, print progress information.}
}
\value{
  An object of class \code{c("structest_t0", "structest", "htest")}.
}
\description{
  Tests whether the structural interpretation of a univariate latent factor
  model can be rejected, using a GMM test that accounts for uncertainty in
  the estimated reliability coefficients (Section 3.2 of VanderWeele and
  Vansteelandt, 2022).
}
\references{
  VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
  reject the structural interpretation of a latent factor model.
  \emph{Journal of the Royal Statistical Society Series B}, 84, 1063--1089.
}
