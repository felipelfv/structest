\name{test_t1}
\alias{test_t1}
\title{T1 Test: Reliability-Independent Test of Structural Interpretation}
\usage{
test_t1(X, z, na.rm = TRUE, max_iter = 1000, tol = 1e-25, verbose = FALSE)
}
\arguments{
  \item{X}{numeric matrix (n x d) of indicator variables, d >= 2.}
  \item{z}{numeric vector of length n encoding the auxiliary variable. Must have at least 3 distinct levels.}
  \item{na.rm}{logical; if \code{TRUE}, rows with any \code{NA} are removed.}
  \item{max_iter}{integer; maximum iterations for the alternating LS procedure.}
  \item{tol}{numeric; convergence tolerance for the alternating LS procedure.}
  \item{verbose}{logical; if \code{TRUE}, print progress information.}
}
\value{
  An object of class \code{c("structest_t1", "structest", "htest")}.
}
\description{
  Tests whether the structural interpretation of a univariate latent factor
  model can be rejected, without requiring estimation of reliability
  coefficients. Uses centered indicators and a two-step efficient GMM
  procedure (Section 3.3 of VanderWeele and Vansteelandt, 2022).
}
\references{
  VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
  reject the structural interpretation of a latent factor model.
  \emph{Journal of the Royal Statistical Society Series B}, 84, 1063--1089.
}
