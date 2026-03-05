#' Print method for structest_fit objects
#'
#' @param x an object of class \code{"structest_fit"}.
#' @param digits number of significant digits for printing.
#' @param ... further arguments passed to or from other methods.
#' @return Invisibly returns \code{x}.
#' @export
print.structest_fit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  type_label <- if (x$type == "t1") {
    "T1: Reliability-independent"
  } else {
    "T0: Reliability-dependent"
  }
  cat("\nStructural Model Fit (", type_label, ")\n\n", sep = "")
  cat("data:  ", x$data.name, "\n")
  cat("n = ", x$n_obs, ", d = ", x$d, " indicators, p = ", x$p,
      " Z-levels\n\n", sep = "")

  # Print coefficient table
  ct <- x$coefficients
  out <- data.frame(
    Estimate = formatC(ct$Estimate, digits = digits, format = "f"),
    Std.Err = formatC(ct$Std.Err, digits = digits, format = "f"),
    `z value` = formatC(ct$z_value, digits = digits, format = "f"),
    `Pr(>|z|)` = format.pval(ct$p_value, digits = digits),
    check.names = FALSE
  )
  rownames(out) <- rownames(ct)
  cat("Coefficients:\n")
  print(out, right = TRUE)

  cat("\nOveridentification test: J = ",
      format(x$test$statistic, digits = digits),
      ", df = ", x$test$df,
      ", p = ", format.pval(x$test$p.value, digits = digits), "\n", sep = "")
  cat("\n")
  invisible(x)
}


#' Summary method for structest_fit objects
#'
#' @param object an object of class \code{"structest_fit"}.
#' @param digits number of significant digits for printing.
#' @param ... further arguments passed to or from other methods.
#' @return Invisibly returns \code{object}.
#' @export
summary.structest_fit <- function(object,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  type_label <- if (object$type == "t1") {
    "T1: Reliability-independent"
  } else {
    "T0: Reliability-dependent"
  }
  cat("\nStructural Model Fit (", type_label, ")\n\n", sep = "")
  cat("data:  ", object$data.name, "\n")
  cat("n = ", object$n_obs, ", d = ", object$d, " indicators, p = ", object$p,
      " Z-levels\n\n", sep = "")

  # Full coefficient table with CIs
  ct <- object$coefficients
  ci_pct <- paste0(format(100 * object$conf_level, digits = 3), "%")
  ci_label <- paste0(ci_pct, " CI")
  out <- data.frame(
    Estimate = formatC(ct$Estimate, digits = digits, format = "f"),
    Std.Err = formatC(ct$Std.Err, digits = digits, format = "f"),
    `z value` = formatC(ct$z_value, digits = digits, format = "f"),
    `Pr(>|z|)` = format.pval(ct$p_value, digits = digits),
    CI = paste0("[", formatC(ct$CI_lower, digits = digits, format = "f"),
                ", ", formatC(ct$CI_upper, digits = digits, format = "f"), "]"),
    check.names = FALSE
  )
  names(out)[5] <- ci_label
  rownames(out) <- rownames(ct)
  cat("Coefficients:\n")
  print(out, right = TRUE)

  cat("\nOveridentification test: J = ",
      format(object$test$statistic, digits = digits),
      ", df = ", object$test$df,
      ", p = ", format.pval(object$test$p.value, digits = digits), "\n",
      sep = "")
  if (object$test$p.value > 0.05) {
    cat("(Fail to reject structural interpretation)\n")
  } else {
    cat("(Reject structural interpretation at 5% level)\n")
  }

  if (object$type == "t0" && !is.null(object$estimates$lambda)) {
    cat("\nEstimated reliabilities (lambda):\n")
    print(round(object$estimates$lambda, digits))
  }

  cat("\nConvergence code: ", object$convergence, "\n\n", sep = "")
  invisible(object)
}


#' Extract coefficients from a structest_fit object
#'
#' @param object an object of class \code{"structest_fit"}.
#' @param ... further arguments (ignored).
#' @return Named numeric vector of parameter estimates.
#' @export
coef.structest_fit <- function(object, ...) {
  setNames(object$coefficients$Estimate, rownames(object$coefficients))
}


#' Confidence intervals for structest_fit parameters
#'
#' @param object an object of class \code{"structest_fit"}.
#' @param parm character vector of parameter names, or numeric indices.
#'   If missing, all parameters are included.
#' @param level confidence level (default uses the level from the fit).
#' @param ... further arguments (ignored).
#' @return A matrix with columns for lower and upper confidence bounds.
#' @importFrom stats qnorm
#' @export
confint.structest_fit <- function(object, parm, level = object$conf_level,
                                  ...) {
  est <- object$coefficients$Estimate
  se <- object$coefficients$Std.Err
  pnames <- rownames(object$coefficients)

  if (missing(parm)) {
    parm <- seq_along(est)
  } else if (is.character(parm)) {
    parm <- match(parm, pnames)
  }

  crit <- qnorm(1 - (1 - level) / 2)
  ci <- cbind(est[parm] - crit * se[parm],
              est[parm] + crit * se[parm])
  pct <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                        digits = 3), " %")
  colnames(ci) <- pct
  rownames(ci) <- pnames[parm]
  ci
}


#' Variance-covariance matrix for structest_fit parameters
#'
#' @param object an object of class \code{"structest_fit"}.
#' @param ... further arguments (ignored).
#' @return The variance-covariance matrix of the parameter estimates.
#' @importFrom stats vcov
#' @export
vcov.structest_fit <- function(object, ...) {
  object$vcov
}
