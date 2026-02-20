#' T1 Test: Reliability-Independent Test of Structural Interpretation (Section 3.3)
#'
#' Tests whether the structural interpretation of a univariate latent factor
#' model can be rejected, without requiring estimation of factor
#' loadings.
#'
#' @param X numeric matrix (n x d) of indicator variables, d >= 2.
#' @param z numeric vector of length n encoding the auxiliary variable.
#'   Must have at least 3 distinct levels.
#' @param na.rm logical; if \code{TRUE}, rows with any \code{NA} are removed.
#' @param max_iter integer; maximum iterations for \code{nlm}.
#' @param tol numeric; convergence tolerance for the alternating LS initializer.
#' @param verbose logical; if \code{TRUE}, print progress information.
#'
#' @return An object of class \code{c("structest_t1", "structest", "htest")}
#'   containing:
#' \describe{
#'   \item{statistic}{the T1 test statistic.}
#'   \item{parameter}{degrees of freedom, \eqn{(d-1)(p-2)}.}
#'   \item{p.value}{p-value from chi-squared distribution.}
#'   \item{method}{description of the test.}
#'   \item{data.name}{name of the data objects.}
#'   \item{estimates}{list with \code{gamma} (\eqn{\gamma_i}), \code{alpha}
#'     (\eqn{\alpha_i}, with \eqn{\alpha_1 = 1}), and \code{beta}
#'     (\eqn{\beta_j}, with \eqn{\beta_1 = 0}).}
#'   \item{n_obs}{number of observations used.}
#'   \item{d}{number of indicators.}
#'   \item{p}{number of Z-levels.}
#'   \item{convergence}{convergence code from \code{nlm}.}
#'   \item{optim_details}{full output from \code{nlm}.}
#' }
#'
#' @details
#' Under the structural model \eqn{X_i = \mu_i + \lambda_i \eta + \varepsilon_i},
#' the conditional expectations satisfy
#' \deqn{E(X_i \mid Z = z_j) = \gamma_i + \alpha_i \beta_j}
#' where \eqn{\gamma_i} are intercepts absorbing the reference-level means,
#' \eqn{\alpha_i} are parameters (with \eqn{\alpha_1 = 1} for
#' identification), and \eqn{\beta_j} are parameters (with
#' \eqn{\beta_1 = 0} for the reference level).
#'
#' This gives \eqn{d \times p} moment conditions:
#' \deqn{E[I(Z = z_j)(X_i - \gamma_i - \alpha_i \beta_j)] = 0}
#' with \eqn{2d + p - 2} free parameters (\eqn{d} intercepts, \eqn{d - 1}
#' alphas, \eqn{p - 1} betas), yielding
#' \eqn{dp - (2d + p - 2) = (d-1)(p-2)} degrees of freedom.
#'
#' The test checks whether the mean-difference matrix
#' \eqn{\Delta_{ij} = E(X_i \mid Z = z_j) - E(X_i \mid Z = z_1)} has rank
#' \eqn{\le 1}, which is the testable implication of the structural model
#' (Theorem 2 in the paper).
#'
#' The procedure first uses two-step GMM to obtain stable initial estimates,
#' then minimizes the GMM criterion with the weight matrix recomputed at each
#' parameter value, matching the authors' reference implementation. The test
#' statistic is asymptotically \eqn{\chi^2_{(d-1)(p-2)}} under the null.
#'
#' @references
#' VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
#' reject the structural interpretation of a latent factor model.
#' \emph{Journal of the Royal Statistical Society Series B}, 84, 2032--2054.
#'
#' @export
test_t1 <- function(X, z, na.rm = TRUE, max_iter = 1000L, tol = 1e-25,
                    verbose = FALSE) {
  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(z)))

  # Validate inputs
  inp <- validate_inputs(X, z, na.rm = na.rm, min_d = 2L, min_p = 3L)
  X <- inp$X; z <- inp$z; n <- inp$n; d <- inp$d; p <- inp$p
  z_levels <- inp$z_levels

  # Build Z indicator matrix (n x p)
  Z_mat <- z_indicator_matrix(z, z_levels)

  # Paper's formulation (Equation 3, p. 2042):
  #   E(X_i | Z = z_j) = gamma_i + alpha_i * beta_j
  #   Identification: alpha_1 = 1 (scale), beta_1 = 0 (reference level)
  #   Moment conditions: d * p
  #   Free parameters: d + (d-1) + (p-1) = 2d + p - 2
  #   df = d*p - (2d + p - 2) = (d-1)(p-2)

  # Initialization from cell means
  # Cell means: mean(X_i | Z = z_j)
  cell_means <- matrix(0, d, p)
  for (j in seq_len(p)) {
    mask <- z == z_levels[j]
    for (i in seq_len(d)) {
      cell_means[i, j] <- mean(X[mask, i])
    }
  }

  # gamma_i = E(X_i | Z = z_1), the reference-level mean
  gamma_init <- cell_means[, 1]

  # Mean differences: D[i,j] = E(X_i|Z=z_j) - gamma_i = alpha_i * beta_j
  D <- cell_means - gamma_init  # d x p, first column is 0

  # beta_j = D[1, j] since alpha_1 = 1
  beta_init <- D[1, ]  # length p, beta_init[1] = 0

  # alpha_i from LS fit of D[i, -1] on beta[-1], for i >= 2
  alpha_init <- numeric(d)
  alpha_init[1] <- 1
  beta_sq_sum <- sum(beta_init[-1]^2)
  if (beta_sq_sum > 1e-10) {
    for (i in 2:d) {
      alpha_init[i] <- sum(D[i, -1] * beta_init[-1]) / beta_sq_sum
    }
  } else {
    alpha_init[2:d] <- 1  # fallback if beta is near zero
  }

  # Parameter vector theta = (gamma_1,...,gamma_d, alpha_2,...,alpha_d, beta_2,...,beta_p)
  # Length: d + (d-1) + (p-1) = 2d + p - 2
  theta_init <- c(gamma_init, alpha_init[-1], beta_init[-1])

  # GMM moment matrix builder
  # Returns n x (d*p) matrix of per-observation moment contributions
  build_gmm <- function(theta_val) {
    gam <- theta_val[1:d]
    alp <- c(1, theta_val[(d + 1):(2 * d - 1)])       # alpha_1 = 1 fixed
    bet <- c(0, theta_val[(2 * d):(2 * d + p - 2)])    # beta_1 = 0 fixed

    # Map each observation to its beta value
    beta_k <- numeric(n)
    for (j in seq_len(p)) {
      beta_k[z == z_levels[j]] <- bet[j]
    }

    # Residual matrix: X[k,i] - gamma_i - alpha_i * beta(z_k)
    resid <- X - matrix(gam, n, d, byrow = TRUE) - outer(beta_k, alp)

    # Moment conditions: g[k, col] = I(z_k == z_j) * resid[k, i]
    g <- matrix(0, nrow = n, ncol = d * p)
    col <- 0L
    for (j in seq_len(p)) {
      for (i in seq_len(d)) {
        col <- col + 1L
        g[, col] <- Z_mat[, j] * resid[, i]
      }
    }
    g
  }

  # Two-step GMM for initial estimates
  # Step 1: weight matrix from initial estimates (fixed during optimization)
  g0 <- build_gmm(theta_init)
  Omega0 <- var(g0) * ((n - 1) / n)
  W1 <- tryCatch(solve(Omega0), error = function(e) {
    solve(Omega0 + 1e-6 * diag(ncol(Omega0)))
  })

  q_step1 <- function(theta_val) {
    g <- build_gmm(theta_val)
    gm <- colMeans(g)
    n * drop(crossprod(gm, W1 %*% gm))
  }

  res1 <- nlm(q_step1, p = theta_init, iterlim = max_iter)
  theta1 <- res1$estimate

  if (verbose) {
    message("Step 1 converged with code: ", res1$code,
            ", J = ", format(res1$minimum, digits = 4))
  }

  # Step 2: efficient weight matrix from step-1 residuals (FIXED)
  g1 <- build_gmm(theta1)
  Omega1 <- var(g1) * ((n - 1) / n)
  W2 <- tryCatch(solve(Omega1), error = function(e) {
    solve(Omega1 + 1e-6 * diag(ncol(Omega1)))
  })

  q_step2 <- function(theta_val) {
    g <- build_gmm(theta_val)
    gm <- colMeans(g)
    n * drop(crossprod(gm, W2 %*% gm))
  }

  res2 <- nlm(q_step2, p = theta1, iterlim = max_iter)

  if (verbose) {
    message("Step 2 converged with code: ", res2$code,
            ", J = ", format(res2$minimum, digits = 4))
  }

  # Final test statistic: weight matrix recomputed at each theta
  q_final <- function(theta_val) {
    g <- build_gmm(theta_val)
    gm <- colMeans(g)
    Sigma <- var(g) * ((n - 1) / n)
    Sigma_inv <- tryCatch(solve(Sigma), error = function(e) {
      solve(Sigma + 1e-6 * diag(ncol(Sigma)))
    })
    n * drop(crossprod(gm, Sigma_inv %*% gm))
  }

  res <- nlm(q_final, p = res2$estimate, iterlim = max_iter)

  if (verbose) {
    message("Final step converged with code: ", res$code,
            ", J = ", format(res$minimum, digits = 4))
  }

  # Degrees of freedom: (d-1)(p-2)
  df <- (d - 1L) * (p - 2L)

  # Extract estimates
  theta_hat <- res$estimate
  gamma_hat <- theta_hat[1:d]
  alpha_hat <- c(1, theta_hat[(d + 1):(2 * d - 1)])
  beta_hat <- c(0, theta_hat[(2 * d):(2 * d + p - 2)])

  names(gamma_hat) <- if (!is.null(colnames(inp$X))) colnames(inp$X) else paste0("X", seq_len(d))
  names(alpha_hat) <- names(gamma_hat)
  names(beta_hat) <- paste0("z=", z_levels)

  result <- list(
    statistic = c(T1 = res$minimum),
    parameter = c(df = df),
    p.value = 1 - pchisq(res$minimum, df = df),
    method = "T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022)",
    data.name = dname,
    estimates = list(gamma = gamma_hat, alpha = alpha_hat, beta = beta_hat),
    n_obs = n,
    d = d,
    p = p,
    convergence = res$code,
    optim_details = res
  )
  class(result) <- c("structest_t1", "structest", "htest")
  result
}
