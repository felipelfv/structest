#' Fit a Structural Model with Standard Errors
#'
#' Fits the structural model from VanderWeele & Vansteelandt (2022) and returns
#' a coefficient table with point estimates, GMM sandwich standard errors,
#' z-values, p-values, and confidence intervals. The overidentification test
#' statistic (T0 or T1) is included as a model diagnostic.
#'
#' @param X numeric matrix (n x d) of indicator variables.
#'   For \code{type = "t1"}, d >= 2; for \code{type = "t0"}, d >= 3.
#' @param z numeric vector or factor of length n encoding a discrete auxiliary
#'   variable. For \code{type = "t1"}, p >= 3 levels; for \code{type = "t0"},
#'   p >= 2 levels.
#' @param type character; \code{"t1"} for the reliability-independent model or
#'   \code{"t0"} for the reliability-dependent model.
#' @param conf_level numeric; confidence level for confidence intervals
#'   (default 0.95).
#' @param na.rm logical; if \code{TRUE}, rows with any \code{NA} are removed.
#' @param max_iter integer; maximum iterations for \code{nlm}.
#' @param verbose logical; if \code{TRUE}, print progress information.
#'
#' @return An object of class \code{c("structest_fit_t1", "structest_fit")} or
#'   \code{c("structest_fit_t0", "structest_fit")} containing:
#' \describe{
#'   \item{coefficients}{data.frame with columns \code{Estimate}, \code{Std.Err},
#'     \code{z_value}, \code{p_value}, \code{CI_lower}, \code{CI_upper}.}
#'   \item{vcov}{variance-covariance matrix of the parameter estimates.}
#'   \item{estimates}{list with named vectors of all parameters (including
#'     constrained ones): \code{gamma}, \code{alpha}/\code{lambda}, \code{beta}.}
#'   \item{test}{list with \code{statistic}, \code{df}, \code{p.value} for the
#'     overidentification test.}
#'   \item{type}{character; \code{"t1"} or \code{"t0"}.}
#'   \item{n_obs}{number of observations.}
#'   \item{d}{number of indicators.}
#'   \item{p}{number of Z-levels.}
#'   \item{data.name}{character string describing the data.}
#'   \item{conf_level}{confidence level used.}
#'   \item{convergence}{convergence code from \code{nlm}.}
#' }
#'
#' @details
#' Standard errors are computed using the GMM sandwich variance estimator:
#' \deqn{\mathrm{Var}(\hat\theta) = \frac{1}{n} (G' \Omega^{-1} G)^{-1}}
#' where \eqn{G = E[\partial g / \partial \theta']} is the analytically derived
#' Jacobian of moment conditions and \eqn{\Omega = \mathrm{Var}(g_k)} is the
#' variance of per-observation moment contributions evaluated at the final
#' estimates.
#'
#' For \code{type = "t0"}, the variance \eqn{\Omega} is additionally adjusted
#' for uncertainty in the reliability estimates, as described in the paper
#' (p. 2041).
#'
#' @references
#' VanderWeele, T. J., & Vansteelandt, S. (2022). A statistical test to
#' reject the structural interpretation of a latent factor model.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology)}, \emph{84}(5), 2032--2054.
#' \doi{10.1111/rssb.12555}
#'
#' @examples
#' set.seed(12345)
#' n <- 1000
#' z <- sample(0:2, n, replace = TRUE)
#' eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
#' lambda <- c(1.0, 0.8, 0.6)
#' X <- cbind(
#'   2 + lambda[1] * eta + rnorm(n, sd = 0.5),
#'   3 + lambda[2] * eta + rnorm(n, sd = 0.5),
#'   1 + lambda[3] * eta + rnorm(n, sd = 0.5)
#' )
#'
#' # Reliability-independent fit
#' fit1 <- fit_structural(X, z, type = "t1")
#' fit1
#' summary(fit1)
#' coef(fit1)
#' confint(fit1)
#'
#' # Reliability-dependent fit
#' fit0 <- fit_structural(X, z, type = "t0")
#' summary(fit0)
#'
#' @importFrom stats qnorm pnorm
#' @export
fit_structural <- function(X, z, type = c("t1", "t0"), conf_level = 0.95,
                           na.rm = TRUE, max_iter = 1000L, verbose = FALSE) {
  type <- match.arg(type)
  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(z)))

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number between 0 and 1.", call. = FALSE)
  }

  # Run the appropriate test to get estimates and overidentification statistic
  test_res <- switch(type,
    t1 = test_t1(X, z, na.rm = na.rm, max_iter = max_iter, verbose = verbose),
    t0 = test_t0(X, z, na.rm = na.rm, max_iter = max_iter, verbose = verbose)
  )

  # Prepare clean data (same validation as test functions)
  inp <- validate_inputs(X, z, na.rm = na.rm,
                         min_d = if (type == "t1") 2L else 3L,
                         min_p = if (type == "t1") 3L else 2L)
  X <- inp$X; z <- inp$z; n <- inp$n; d <- inp$d; p <- inp$p
  z_levels <- inp$z_levels
  Z_mat <- z_indicator_matrix(z, z_levels)
  p_z <- colMeans(Z_mat)

  crit <- qnorm(1 - (1 - conf_level) / 2)

  if (type == "t1") {
    se_result <- .sandwich_se_t1(test_res, X, z, z_levels, Z_mat, p_z, n, d, p)
  } else {
    se_result <- .sandwich_se_t0(test_res, inp, X, z, z_levels, Z_mat, p_z, n, d, p)
  }

  estimates_vec <- se_result$estimates_vec
  vcov_mat <- se_result$vcov_mat
  param_names <- se_result$param_names
  se <- sqrt(pmax(diag(vcov_mat), 0))

  # Build coefficient table
  z_val <- estimates_vec / se
  p_val <- 2 * pnorm(-abs(z_val))
  ci_lower <- estimates_vec - crit * se
  ci_upper <- estimates_vec + crit * se

  coef_table <- data.frame(
    Estimate = estimates_vec,
    Std.Err = se,
    z_value = z_val,
    p_value = p_val,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    row.names = param_names
  )

  subclass <- paste0("structest_fit_", type)

  structure(list(
    coefficients = coef_table,
    vcov = vcov_mat,
    estimates = test_res$estimates,
    test = list(
      statistic = unname(test_res$statistic),
      df = unname(test_res$parameter),
      p.value = test_res$p.value
    ),
    type = type,
    n_obs = n,
    d = d,
    p = p,
    data.name = dname,
    conf_level = conf_level,
    convergence = test_res$convergence
  ), class = c(subclass, "structest_fit"))
}


# Internal helpers for sandwich SE computation

#' Sandwich SEs for T1 (reliability-independent) model
#' @noRd
.sandwich_se_t1 <- function(test_res, X, z, z_levels, Z_mat, p_z, n, d, p) {
  gamma_hat <- test_res$estimates$gamma
  alpha_hat <- test_res$estimates$alpha
  beta_hat <- test_res$estimates$beta

  # Reconstruct moment conditions at final estimates
  beta_k <- beta_hat[match(z, z_levels)]
  resid <- X - matrix(gamma_hat, n, d, byrow = TRUE) - outer(beta_k, alpha_hat)

  # Column ordering: for j in 1:p, for i in 1:d
  g <- matrix(0, nrow = n, ncol = d * p)
  col <- 0L
  for (j in seq_len(p)) {
    for (i in seq_len(d)) {
      col <- col + 1L
      g[, col] <- Z_mat[, j] * resid[, i]
    }
  }

  Omega <- var(g) * ((n - 1) / n)
  Omega_inv <- tryCatch(solve(Omega), error = function(e) {
    solve(Omega + 1e-6 * diag(ncol(Omega)))
  })

  # Analytical Jacobian G: (d*p) x (2d+p-2)
  # Parameter order: gamma_1,...,gamma_d, alpha_2,...,alpha_d, beta_2,...,beta_p
  n_params <- 2L * d + p - 2L
  n_moments <- d * p
  G <- matrix(0, nrow = n_moments, ncol = n_params)

  col_m <- 0L
  for (j in seq_len(p)) {
    for (i in seq_len(d)) {
      col_m <- col_m + 1L
      # dg/d(gamma_i)
      G[col_m, i] <- -p_z[j]
      # dg/d(alpha_i) for i > 1
      if (i > 1L) {
        G[col_m, d + (i - 1L)] <- -p_z[j] * beta_hat[j]
      }
      # dg/d(beta_j) for j > 1
      if (j > 1L) {
        G[col_m, 2L * d - 1L + (j - 1L)] <- -p_z[j] * alpha_hat[i]
      }
    }
  }

  # Sandwich variance: (1/n) * (G' Omega^{-1} G)^{-1}
  bread <- crossprod(G, Omega_inv %*% G)
  vcov_mat <- solve(bread) / n

  # Parameter names
  ind_names <- names(gamma_hat)
  param_names <- c(
    paste0("gamma[", ind_names, "]"),
    paste0("alpha[", ind_names[-1], "]"),
    paste0("beta[", names(beta_hat)[-1], "]")
  )

  estimates_vec <- c(gamma_hat, alpha_hat[-1], beta_hat[-1])
  names(estimates_vec) <- param_names
  rownames(vcov_mat) <- colnames(vcov_mat) <- param_names

  list(estimates_vec = estimates_vec, vcov_mat = vcov_mat,
       param_names = param_names)
}


#' Sandwich SEs for T0 (reliability-dependent) model
#' @noRd
.sandwich_se_t0 <- function(test_res, inp, X, z, z_levels, Z_mat, p_z,
                            n, d, p) {
  gamma_hat <- test_res$estimates$gamma
  beta_hat <- test_res$estimates$beta
  lambda <- test_res$estimates$lambda
  ratio <- lambda / lambda[1]

  # Re-estimate reliability to get V_k, U_pairs for variance adjustment
  rel <- estimate_reliability(inp$X, na.rm = FALSE)
  V_k <- rel$V_k
  U_pairs <- rel$U_pairs
  pairs <- rel$pairs
  mean_U_pairs <- colMeans(U_pairs)

  # Reconstruct moment conditions (column ordering: for i in 1:d, for j in 1:p)
  u <- matrix(0, nrow = n, ncol = d * p)
  col <- 0L
  for (i in seq_len(d)) {
    for (j in seq_len(p)) {
      col <- col + 1L
      pred <- gamma_hat[i] + ratio[i] * beta_hat[j]
      u[, col] <- Z_mat[, j] * (X[, i] - pred)
    }
  }

  # Reliability variance adjustment (same as test_t0)
  diag_pz <- numeric(d * p)
  col <- 0L
  for (i in seq_len(d)) {
    for (j in seq_len(p)) {
      col <- col + 1L
      diag_pz[col] <- p_z[j]
    }
  }

  dpred_dlambda <- matrix(0, nrow = d * p, ncol = d)
  col <- 0L
  for (i in seq_len(d)) {
    for (j in seq_len(p)) {
      col <- col + 1L
      if (i > 1L && j > 1L) {
        dpred_dlambda[col, 1] <- -ratio[i] * beta_hat[j] / lambda[1]
        dpred_dlambda[col, i] <- beta_hat[j] / lambda[1]
      }
    }
  }

  dudlambda <- -diag(diag_pz) %*% dpred_dlambda

  dvdlambda <- matrix(0, nrow = d, ncol = d)
  for (i in seq_len(d)) {
    dvdlambda[i, i] <- -sum(lambda[-i]^2)
    for (k in seq_len(d)) {
      if (k != i) {
        pair_idx <- which((pairs[, 1] == min(i, k)) & (pairs[, 2] == max(i, k)))
        dvdlambda[i, k] <- mean_U_pairs[pair_idx] - lambda[i] * lambda[k]
      }
    }
  }

  dvdlambda_inv <- tryCatch(solve(dvdlambda), error = function(e) {
    solve(dvdlambda + 1e-6 * diag(ncol(dvdlambda)))
  })
  uadjust <- u - t(dudlambda %*% dvdlambda_inv %*% t(V_k))

  Omega <- var(uadjust) * ((n - 1) / n)
  Omega_inv <- tryCatch(solve(Omega), error = function(e) {
    solve(Omega + 1e-6 * diag(ncol(Omega)))
  })

  # Analytical Jacobian G: (d*p) x (d+p-1)
  # Parameter order: gamma_1,...,gamma_d, beta_2,...,beta_p
  n_params <- d + p - 1L
  n_moments <- d * p
  G <- matrix(0, nrow = n_moments, ncol = n_params)

  col_m <- 0L
  for (i in seq_len(d)) {
    for (j in seq_len(p)) {
      col_m <- col_m + 1L
      # dg/d(gamma_i)
      G[col_m, i] <- -p_z[j]
      # dg/d(beta_t) for j > 1
      if (j > 1L) {
        G[col_m, d + (j - 1L)] <- -p_z[j] * ratio[i]
      }
    }
  }

  # Sandwich variance
  bread <- crossprod(G, Omega_inv %*% G)
  vcov_mat <- solve(bread) / n

  # Parameter names
  ind_names <- names(gamma_hat)
  param_names <- c(
    paste0("gamma[", ind_names, "]"),
    paste0("beta[", names(beta_hat)[-1], "]")
  )

  estimates_vec <- c(gamma_hat, beta_hat[-1])
  names(estimates_vec) <- param_names
  rownames(vcov_mat) <- colnames(vcov_mat) <- param_names

  list(estimates_vec = estimates_vec, vcov_mat = vcov_mat,
       param_names = param_names)
}
