#' Print method for structest objects
#'
#' @param x an object of class \code{"structest"}.
#' @param digits number of significant digits for printing.
#' @param ... further arguments passed to or from other methods.
#' @return Invisibly returns \code{x}.
#' @export
print.structest <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n\t", x$method, "\n\n")
  cat("data:  ", x$data.name, "\n")
  cat("statistic = ", format(x$statistic, digits = digits),
      ", df = ", x$parameter,
      ", p-value = ", format.pval(x$p.value, digits = digits), "\n", sep = "")
  cat("n = ", x$n_obs, ", d = ", x$d, " indicators, p = ", x$p, " Z-levels\n",
      sep = "")
  cat("\n")
  invisible(x)
}


#' Summary method for structest objects
#'
#' @param object an object of class \code{"structest"}.
#' @param digits number of significant digits for printing.
#' @param ... further arguments passed to or from other methods.
#' @return Invisibly returns \code{object}.
#' @export
summary.structest <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n\t", object$method, "\n\n")
  cat("data:  ", object$data.name, "\n\n")
  cat("Test statistic: ", format(object$statistic, digits = digits), "\n")
  cat("Degrees of freedom: ", object$parameter, "\n")
  cat("P-value: ", format.pval(object$p.value, digits = digits), "\n\n")
  cat("Sample size: ", object$n_obs, "\n")
  cat("Indicators (d): ", object$d, "\n")
  cat("Z-levels (p): ", object$p, "\n")
  cat("Convergence code: ", object$convergence, "\n\n")

  cat("Parameter estimates:\n")
  if (!is.null(object$estimates$gamma)) {
    cat("  gamma (intercepts):\n")
    print(round(object$estimates$gamma, digits))
  }
  if (!is.null(object$estimates$alpha)) {
    cat("  alpha (centered loadings):\n")
    print(round(object$estimates$alpha, digits))
  }
  if (!is.null(object$estimates$beta)) {
    cat("  beta (Z effects):\n")
    print(round(object$estimates$beta, digits))
  }
  if (!is.null(object$estimates$lambda)) {
    cat("  lambda (reliabilities):\n")
    print(round(object$estimates$lambda, digits))
  }
  cat("\n")
  invisible(object)
}
