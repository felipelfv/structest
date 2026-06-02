#' Equivalence Test for the Structural Interpretation
#'
#' Inverts the burden of proof of \code{\link{test_t0}} / \code{\link{test_t1}}.
#' Those tests can only \emph{reject} the structural interpretation; a large
#' p-value never confirms it (absence of evidence is not evidence of absence).
#' This equivalence test asks the complementary question: are the deviations
#' from the structural model small enough to be considered negligible, relative
#' to a pre-specified margin?
#'
#' @param object a fitted \code{"structest"} object, as returned by
#'   \code{\link{test_t0}} or \code{\link{test_t1}}.
#' @param margin numeric scalar, the equivalence margin on the standardised
#'   misfit \eqn{\varepsilon} (an RMSEA-style index, see Details). There is no
#'   default: the margin is a substantive judgement about how much departure
#'   from the structural interpretation is tolerable, and must be supplied
#'   explicitly.
#' @param alpha numeric significance level for the equivalence test
#'   (default 0.05).
#'
#' @return An object of class \code{c("structest_equiv", "structest", "htest")}
#'   containing:
#' \describe{
#'   \item{statistic}{the underlying T0/T1 distance-metric statistic \eqn{J}.}
#'   \item{parameter}{degrees of freedom \eqn{df}.}
#'   \item{p.value}{equivalence p-value (lower-tail probability of a non-central
#'     \eqn{\chi^2}; small values reject poor fit).}
#'   \item{margin}{the supplied equivalence margin \eqn{\varepsilon_0}.}
#'   \item{epsilon}{point estimate of the standardised misfit \eqn{\hat\varepsilon}.}
#'   \item{epsilon_upper}{one-sided \eqn{(1-\alpha)} upper confidence limit on
#'     \eqn{\varepsilon}. Equivalence holds iff this is below \code{margin}.}
#'   \item{equivalence}{logical; \code{TRUE} if poor fit is rejected at level
#'     \code{alpha}.}
#'   \item{method, data.name, n_obs, d, p}{carried over from \code{object}.}
#' }
#'
#' @details
#' The T0/T1 distance-metric statistic \eqn{J} is asymptotically
#' \eqn{\chi^2_{df}} under the null that the structural interpretation holds, and
#' non-central \eqn{\chi^2_{df}(\lambda)} under local departures from it
#' (Newey, 1985). The non-centrality \eqn{\lambda} indexes population misfit. We
#' summarise it on an RMSEA-style scale,
#' \deqn{\varepsilon = \sqrt{\lambda / (n\,df)},}
#' which is roughly invariant to sample size. The margin \eqn{\varepsilon_0}
#' corresponds to a boundary non-centrality \eqn{\lambda_0 = n\,df\,\varepsilon_0^2}.
#'
#' The test is a one-sided ("close-fit") equivalence test:
#' \deqn{H_0:\ \varepsilon \ge \varepsilon_0 \quad\text{vs}\quad
#'       H_1:\ \varepsilon < \varepsilon_0.}
#' Rejecting \eqn{H_0} licenses the positive claim that departures from the
#' structural interpretation are smaller than the margin (so the
#' \code{\link{fit_structural}} estimates may be read structurally). The
#' equivalence p-value is
#' \eqn{P(\chi^2_{df}(\lambda_0) \le J_{\text{obs}})}, and the decision is
#' equivalent to checking whether \code{epsilon_upper} falls below \code{margin}.
#'
#' Note: unlike SEM, there is no canonical \eqn{\varepsilon} threshold for
#' "approximate structural interpretation". The margin is yours to justify, and
#' the non-central \eqn{\chi^2} approximation is asymptotic; calibration by
#' simulation under known local misfit is advisable.
#'
#' @references
#' VanderWeele, T. J., & Vansteelandt, S. (2022). A statistical test to
#' reject the structural interpretation of a latent factor model.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology)}, \emph{84}(5), 2032--2054. \doi{10.1111/rssb.12555}
#'
#' Newey, W. K. (1985). Generalized method of moments specification testing.
#' \emph{Journal of Econometrics}, \emph{29}(3), 229--256.
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and
#' determination of sample size for covariance structure modeling.
#' \emph{Psychological Methods}, \emph{1}(2), 130--149.
#'
#' @seealso \code{\link{test_t0}}, \code{\link{test_t1}} for the rejection tests
#'   this inverts; \code{\link{fit_structural}} for the structural estimates.
#'
#' @examples
#' set.seed(12345)
#' n <- 1000
#' z <- rbinom(n, 1, 0.5)
#' eta <- 1 + 0.5 * z + rnorm(n)
#' lambda <- c(1.0, 0.8, 0.6)
#' X <- cbind(
#'   2 + lambda[1] * eta + rnorm(n, sd = 0.5),
#'   3 + lambda[2] * eta + rnorm(n, sd = 0.5),
#'   1 + lambda[3] * eta + rnorm(n, sd = 0.5)
#' )
#'
#' fit <- test_t0(X, z)
#' equiv_test(fit, margin = 0.05)
#'
#' @importFrom stats pchisq uniroot
#' @export
equiv_test <- function(object, margin, alpha = 0.05) {
  if (!inherits(object, "structest")) {
    stop("'object' must be a fitted structest object (from test_t0() or test_t1()).",
         call. = FALSE)
  }
  if (missing(margin)) {
    stop("'margin' must be supplied: choose a substantive equivalence margin ",
         "on the standardised misfit epsilon (there is no canonical default).",
         call. = FALSE)
  }
  if (!is.numeric(margin) || length(margin) != 1L || is.na(margin) || margin <= 0) {
    stop("'margin' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }

  J  <- unname(object$statistic)
  df <- unname(object$parameter)
  n  <- object$n_obs

  if (is.na(df) || df < 1) {
    stop("the fitted model is just-identified (df < 1); there is no ",
         "over-identifying misfit to test for equivalence.", call. = FALSE)
  }

  # Boundary non-centrality implied by the margin, and the close-fit p-value.
  ncp0    <- n * df * margin^2
  p.value <- pchisq(J, df = df, ncp = ncp0)

  # Point estimate of standardised misfit (clamped at 0).
  epsilon <- sqrt(max(0, (J - df) / (n * df)))

  # One-sided (1 - alpha) upper confidence limit on epsilon: the ncp at which
  # pchisq(J, df, ncp) == alpha. pchisq is decreasing in ncp, so if even ncp = 0
  # already puts J in the lower alpha tail, the limit is 0.
  if (pchisq(J, df = df, ncp = 0) <= alpha) {
    epsilon_upper <- 0
  } else {
    hi <- max(J, df) + 10
    while (pchisq(J, df = df, ncp = hi) > alpha) hi <- hi * 2
    ncp_u <- uniroot(function(ncp) pchisq(J, df = df, ncp = ncp) - alpha,
                     interval = c(0, hi))$root
    epsilon_upper <- sqrt(ncp_u / (n * df))
  }

  equivalence <- p.value < alpha

  test_label <- names(object$statistic)
  if (is.null(test_label)) test_label <- "J"

  result <- list(
    statistic   = object$statistic,
    parameter   = object$parameter,
    p.value     = p.value,
    margin      = margin,
    alpha       = alpha,
    epsilon     = epsilon,
    epsilon_upper = epsilon_upper,
    equivalence = equivalence,
    method      = paste0("Equivalence test for the structural interpretation ",
                         "(close-fit, based on ", test_label, ")"),
    data.name   = object$data.name,
    n_obs       = n,
    d           = object$d,
    p           = object$p,
    base_test   = test_label
  )
  class(result) <- c("structest_equiv", "structest", "htest")
  result
}


#' Print method for structest equivalence tests
#'
#' @param x an object of class \code{"structest_equiv"}.
#' @param digits number of significant digits for printing.
#' @param ... further arguments passed to or from other methods.
#' @return Invisibly returns \code{x}.
#' @export
print.structest_equiv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n\t", x$method, "\n\n")
  cat("data:  ", x$data.name, "\n")
  cat("H0: misfit epsilon >= ", format(x$margin, digits = digits),
      "  vs  H1: epsilon < ", format(x$margin, digits = digits), "\n", sep = "")
  cat(x$base_test, " = ", format(x$statistic, digits = digits),
      ", df = ", x$parameter,
      ", equivalence p-value = ", format.pval(x$p.value, digits = digits),
      "\n", sep = "")
  cat("epsilon (estimate) = ", format(x$epsilon, digits = digits),
      ", ", round(100 * (1 - x$alpha)), "% upper limit = ",
      format(x$epsilon_upper, digits = digits), "\n", sep = "")
  cat("\n")
  if (x$equivalence) {
    cat("Equivalence supported at alpha = ", format(x$alpha, digits = digits),
        ": departures from the structural interpretation are below the margin.\n",
        sep = "")
  } else {
    cat("Equivalence NOT supported at alpha = ", format(x$alpha, digits = digits),
        ": cannot rule out misfit at or beyond the margin.\n", sep = "")
  }
  cat("\n")
  invisible(x)
}
