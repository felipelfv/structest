\name{estimate_reliability}
\alias{estimate_reliability}
\title{Estimate Reliability Coefficients (Section 3.1)}
\usage{
estimate_reliability(X, na.rm = TRUE)
}
\arguments{
  \item{X}{numeric matrix (n x d) of indicator variables, d >= 3.}
  \item{na.rm}{logical; if \code{TRUE}, rows with any \code{NA} are removed.}
}
\value{
  A list with components:
  \describe{
    \item{lambda}{numeric vector of estimated reliabilities.}
    \item{glm_fit}{the fitted \code{glm} object.}
    \item{V_k}{matrix of per-subject estimating function contributions.}
  }
}
\description{
  Estimates the reliability coefficients for each indicator in a univariate
  latent factor model, using pairwise covariances and a quasi-Poisson GLM
  with log link (Section 3.1 of VanderWeele and Vansteelandt, 2022).
}
\references{
  VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
  reject the structural interpretation of a latent factor model.
  \emph{Journal of the Royal Statistical Society Series B}, 84, 1063--1089.
}
