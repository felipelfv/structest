#' structest: Testing the Structural Interpretation of a Latent Factor Model
#'
#' Implements two statistical tests from VanderWeele and Vansteelandt (2022)
#' to reject the structural interpretation of a univariate latent factor model.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{test_t0}}}{Reliability-dependent test (Section 3.2).
#'     Requires d >= 3 indicators and p >= 2 Z-levels.}
#'   \item{\code{\link{test_t1}}}{Reliability-independent test (Section 3.3).
#'     Requires d >= 2 indicators and p >= 3 Z-levels.}
#'   \item{\code{\link{estimate_reliability}}}{Estimate reliability coefficients
#'     via quasi-Poisson GLM (Section 3.1).}
#' }
#'
#' @references
#' VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
#' reject the structural interpretation of a latent factor model.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 84, 2032--2054.
#'
#' @docType package
#' @name structest-package
"_PACKAGE"
