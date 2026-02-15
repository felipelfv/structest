#' T1 Test: Reliability-Independent Test of Structural Interpretation (Section 3.3)
#'
#' Tests whether the structural interpretation of a univariate latent factor
#' model can be rejected, without requiring estimation of reliability
#' coefficients. Uses centered indicators and a two-step efficient GMM procedure.
#'
#' @param X numeric matrix (n x d) of indicator variables, d >= 3.
#' @param z numeric vector of length n encoding the auxiliary variable.
#'   Must have at least 3 distinct levels.
#' @param na.rm logical; if \code{TRUE}, rows with any \code{NA} are removed.
#' @param max_iter integer; maximum iterations for the alternating LS procedure.
#' @param tol numeric; convergence tolerance for the alternating LS procedure.
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
#'   \item{estimates}{list with \code{alpha} and \code{beta}.}
#'   \item{n_obs}{number of observations used.}
#'   \item{d}{number of indicators.}
#'   \item{p}{number of Z-levels.}
#'   \item{convergence}{convergence code from \code{nlm}.}
#'   \item{optim_details}{full output from \code{nlm}.}
#' }
#'
#' @details
#' The indicators are centered by their row mean:
#' \eqn{X^c_{ik} = X_{ik} - \bar{X}_{i\cdot}}, which eliminates the latent
#' factor \eqn{Y}. Under the structural interpretation,
#' \eqn{E[X^c_i | Z = z_j] = \alpha_i \beta_j} where \eqn{\alpha_i} and
#' \eqn{\beta_j} are identifiable nuisance parameters. The test has
#' \eqn{(d-1)(p-2)} degrees of freedom. The full model (Equation 3 in the
#' paper) has \eqn{d \times p} moment conditions and \eqn{2d + p - 2} free
#' parameters (\eqn{d} intercepts, \eqn{d} alphas, \eqn{p-1} betas, minus 1
#' for the alpha-beta scale indeterminacy), giving
#' \eqn{dp - (2d+p-2) = (d-1)(p-2)} degrees of freedom.
#'
#' The procedure uses a two-step efficient GMM: first with a fixed weight
#' matrix from initial estimates, then with the efficient weight matrix.
#'
#' @references
#' VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
#' reject the structural interpretation of a latent factor model.
#' \emph{Journal of the Royal Statistical Society Series B}, 84, 1063--1089.
#'
#' @export
test_t1 <- function(X, z, na.rm = TRUE, max_iter = 1000L, tol = 1e-25,
                    verbose = FALSE) {
  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(z)))

  # Validate inputs
  inp <- validate_inputs(X, z, na.rm = na.rm, min_d = 3L, min_p = 3L)
  X <- inp$X; z <- inp$z; n <- inp$n; d <- inp$d; p <- inp$p
  z_levels <- inp$z_levels

  # Step 1: Center indicators by row mean
  row_means <- rowMeans(X)
  Xc <- X - row_means  # n x d, centered

  # Build Z indicator matrix
  Z_mat <- z_indicator_matrix(z, z_levels)

  # Step 2: Initial estimates via alternating LS (Appendix of paper)
  alpha <- numeric(d)
  beta <- numeric(p)

  # Initialize alpha from the first Z-level
  for (i in seq_len(d)) {
    alpha[i] <- mean(Xc[z == z_levels[1], i], na.rm = TRUE)
  }

  # Alternating LS iterations
  # Model: E[Xc_i | z == z_j] = alpha_i * beta_j
  xalpha <- Xc  # will hold Xc * alpha (broadcasted)
  valpha <- Xc  # will hold alpha values (broadcasted)

  # Compute alpha-weighted quantities
  for (i in seq_len(d)) {
    xalpha[, i] <- Xc[, i] * alpha[i]
    valpha[, i] <- alpha[i]
  }

  su <- 10
  iter <- 0L

  repeat {
    iter <- iter + 1L

    # Update beta: for each Z-level j,
    # beta_j = sum_k sum_i Xc_ik * alpha_i * I(z_k == z_j) /
    #          sum_k sum_i alpha_i^2 * I(z_k == z_j)
    for (j in seq_len(p)) {
      mask <- z == z_levels[j]
      num <- sum(xalpha[mask, ], na.rm = TRUE)
      den <- sum(valpha[mask, ]^2, na.rm = TRUE)
      beta[j] <- num / den
    }

    # Compute beta-weighted quantities
    xbeta <- Xc
    vbeta <- Xc
    for (j in seq_len(p)) {
      mask <- z == z_levels[j]
      for (i in seq_len(d)) {
        xbeta[mask, i] <- Xc[mask, i] * beta[j]
        vbeta[mask, i] <- beta[j]
      }
    }

    # Update alpha: for each indicator i,
    # alpha_i = sum_k Xc_ik * beta(z_k) / sum_k beta(z_k)^2
    for (i in seq_len(d)) {
      num <- sum(xbeta[, i], na.rm = TRUE)
      den <- sum(vbeta[, i]^2, na.rm = TRUE)
      alpha[i] <- num / den
      xalpha[, i] <- Xc[, i] * alpha[i]
      valpha[, i] <- alpha[i]
    }

    # Check convergence using estimating equations
    u_check <- numeric(d + p)
    for (i in seq_len(d)) {
      u_check[i] <- sum(xbeta[, i] - valpha[, i] * vbeta[, i]^2, na.rm = TRUE)
    }
    for (j in seq_len(p)) {
      mask <- z == z_levels[j]
      u_check[d + j] <- sum(xalpha[mask, ] - valpha[mask, ]^2 * vbeta[mask, ],
                            na.rm = TRUE)
    }

    su_new <- sum(u_check^2)
    if (verbose) message("Iteration ", iter, ": criterion = ", su_new)
    if (abs(su_new - su) < tol || iter >= max_iter) break
    su <- su_new
  }

  # Save initial theta
  theta <- c(alpha, beta)

  # Build GMM moment matrix for given theta.
  # Moment conditions: E[I(z==z_j)*(Xc_i - alpha_i*beta_j)] = 0
  # for i=2,...,d (indicator 1 dropped) and j=1,...,p, giving (d-1)*p conditions.
  # df = (d-1)(p-2) as per Eq. 3 in the paper.

  build_gmm_matrix <- function(theta_val, Xc_mat, z_vec, z_lvls, d_val, p_val, n_val) {
    alpha_val <- theta_val[seq_len(d_val)]
    beta_val <- theta_val[(d_val + 1):(d_val + p_val)]

    # Map each observation to its beta value
    beta_k <- numeric(n_val)
    for (j in seq_len(p_val)) {
      beta_k[z_vec == z_lvls[j]] <- beta_val[j]
    }

    resid <- Xc_mat - outer(beta_k, alpha_val)  # n x d

    # Drop first indicator
    resid <- resid[, -1, drop = FALSE]  # n x (d-1)

    # Moment conditions: for each z-level j, multiply residuals by I(z==z_j)
    Z_ind <- z_indicator_matrix(z_vec, z_lvls)
    g <- matrix(0, nrow = n_val, ncol = (d_val - 1L) * p_val)
    col <- 0L
    for (j in seq_len(p_val)) {
      for (i in seq_len(d_val - 1L)) {
        col <- col + 1L
        g[, col] <- Z_ind[, j] * resid[, i]
      }
    }
    g
  }

  # Step 3: First-step GMM with fixed weight matrix
  g_init <- build_gmm_matrix(theta, Xc, z, z_levels, d, p, n)
  W <- solve(var(g_init, na.rm = TRUE))

  q_step1 <- function(theta_val) {
    g <- build_gmm_matrix(theta_val, Xc, z, z_levels, d, p, n)
    gm <- colMeans(g, na.rm = TRUE)
    n * drop(t(gm) %*% W %*% gm)
  }

  res1 <- nlm(q_step1, p = theta)
  theta <- res1$estimate

  # Step 4: Second-step efficient GMM with updated weight matrix
  q_efficient <- function(theta_val) {
    g <- build_gmm_matrix(theta_val, Xc, z, z_levels, d, p, n)
    gm <- colMeans(g, na.rm = TRUE)
    W_eff <- solve(var(g, na.rm = TRUE))
    n * drop(t(gm) %*% W_eff %*% gm)
  }

  res <- nlm(q_efficient, p = theta)

  if (verbose) {
    message("T1 optimization converged with code: ", res$code)
    message("Minimum: ", res$minimum)
  }

  # Degrees of freedom: (d-1)(p-2), as per the paper (p. 2042).
  # The full model has d*p moment conditions and 2d+p-2 free parameters,
  # giving df = d*p - (2d+p-2) = (d-1)(p-2).
  df <- (d - 1L) * (p - 2L)

  # Extract estimates
  alpha_hat <- res$estimate[seq_len(d)]
  beta_hat <- res$estimate[(d + 1):(d + p)]
  names(alpha_hat) <- if (!is.null(colnames(inp$X))) colnames(inp$X) else paste0("X", seq_len(d))
  names(beta_hat) <- paste0("z=", z_levels)

  result <- list(
    statistic = c(T1 = res$minimum),
    parameter = c(df = df),
    p.value = 1 - pchisq(res$minimum, df = df),
    method = "T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022)",
    data.name = dname,
    estimates = list(alpha = alpha_hat, beta = beta_hat),
    n_obs = n,
    d = d,
    p = p,
    convergence = res$code,
    optim_details = res
  )
  class(result) <- c("structest_t1", "structest", "htest")
  result
}
