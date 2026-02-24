#' T0 Test: Reliability-Dependent Test of Structural Interpretation (Section 3.2)
#'
#' Tests whether the structural interpretation of a univariate latent factor
#' model can be rejected, using estimated reliability coefficients. The
#' variance of the estimating equations is adjusted for uncertainty in the
#' reliability estimates (see paper, p. 2041).
#'
#' @param X numeric matrix (n x d) of indicator variables, d >= 3.
#' @param z numeric vector of length n encoding the auxiliary variable.
#'   Must have at least 2 distinct levels.
#' @param na.rm logical; if \code{TRUE}, rows with any \code{NA} are removed.
#' @param verbose logical; if \code{TRUE}, print progress information.
#'
#' @return An object of class \code{c("structest_t0", "structest", "htest")}
#'   containing:
#' \describe{
#'   \item{statistic}{the T0 test statistic.}
#'   \item{parameter}{degrees of freedom, \eqn{(d-1)(p-1)}.}
#'   \item{p.value}{p-value from chi-squared distribution.}
#'   \item{method}{description of the test.}
#'   \item{data.name}{name of the data objects.}
#'   \item{estimates}{list with \code{gamma} (intercepts), \code{beta} (Z effects),
#'     and \code{lambda} (reliabilities).}
#'   \item{n_obs}{number of observations used.}
#'   \item{d}{number of indicators.}
#'   \item{p}{number of Z-levels.}
#'   \item{convergence}{convergence code from \code{nlm}.}
#'   \item{optim_details}{full output from \code{nlm}.}
#' }
#'
#' @details
#' Under the structural model, \eqn{E[X_i | Z = z_j] = \gamma_i + (\lambda_i / \lambda_1) \beta_j}.
#' This gives \eqn{d \times p} moment conditions with \eqn{d + p - 1} free parameters
#' (\eqn{\gamma_1, \ldots, \gamma_d} and \eqn{\beta_2, \ldots, \beta_p} with \eqn{\beta_1 = 0}),
#' yielding \eqn{(d-1)(p-1)} degrees of freedom.
#'
#' The variance of the estimating equations is adjusted for the
#' uncertainty in the reliability estimates \eqn{\lambda_i}
#' (see paper, p. 2041). Consistent generalised methods of moments
#' estimators (Newey & McFadden, 1994) are obtained by minimising
#' the resulting distance metric statistic.
#'
#' @references
#' VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
#' reject the structural interpretation of a latent factor model.
#' \emph{Journal of the Royal Statistical Society Series B}, 84, 2032--2054.
#'
#' @export
test_t0 <- function(X, z, na.rm = TRUE, verbose = FALSE) {
  dname <- paste(deparse(substitute(X)), "and", deparse(substitute(z)))

  # Validate inputs
  inp <- validate_inputs(X, z, na.rm = na.rm, min_d = 3L, min_p = 2L)
  X <- inp$X; z <- inp$z; n <- inp$n; d <- inp$d; p <- inp$p
  z_levels <- inp$z_levels

  # Step 1: Estimate reliabilities (Section 3.1)
  rel <- estimate_reliability(X, na.rm = FALSE)
  lambda <- rel$lambda
  V_k <- rel$V_k         # n x d: per-subject reliability scores
  U_pairs <- rel$U_pairs  # n x choose(d,2): pairwise estimating functions
  pairs <- rel$pairs
  mean_U_pairs <- colMeans(U_pairs)

  # Step 2: Fit constrained model for initial parameter estimates
  # theta = (gamma_1, ..., gamma_d, beta_2, ..., beta_p)
  # with beta_1 = 0 (reference level = first z-level)
  # Model: E[X_i | z=z_j] = gamma_i + (lambda_i/lambda_1) * beta_j

  # Build design matrix for constrained lm
  x_long <- as.vector(X)
  y_long <- rep(seq_len(d), each = n)
  z_long <- rep(z, times = d)

  # For each z-level j>1, create an indicator I(z=z_j) * lambda_i/lambda_1
  # This gives p-1 regressors
  ratio <- lambda / lambda[1]
  Zy_list <- list()
  for (j in 2:p) {
    Zy_list[[j - 1]] <- as.numeric(z_long == z_levels[j]) * rep(ratio, each = n)
  }
  Zy_mat <- do.call(cbind, Zy_list)

  # lm: x_long ~ -1 + factor(y_long) + Zy_mat
  design <- cbind(model.matrix(~ -1 + factor(y_long)), Zy_mat)
  mod <- lm.fit(design, x_long)
  theta_init <- mod$coefficients  # d + (p-1) values

  # Build Z indicator matrix
  Z_mat <- z_indicator_matrix(z, z_levels)
  p_z <- colMeans(Z_mat)  # proportion in each z-level

  # Step 3: Distance metric statistic with variance adjustment for reliability uncertainty
  # u_{ij} = I(z==z_j) * (X_i - gamma_i - (lambda_i/lambda_1)*beta_j), beta_1 = 0

  q_t0 <- function(theta) {
    gamma <- theta[seq_len(d)]
    beta <- c(0, theta[(d + 1):(d + p - 1)])  # beta_1=0, beta_2,...,beta_p

    # Build n x (d*p) matrix of moment conditions
    # Column ordering: (i=1,j=1), (i=1,j=2), ..., (i=1,j=p), (i=2,j=1), ...
    u <- matrix(0, nrow = n, ncol = d * p)
    col <- 0L
    for (i in seq_len(d)) {
      for (j in seq_len(p)) {
        col <- col + 1L
        pred <- gamma[i] + ratio[i] * beta[j]
        u[, col] <- Z_mat[, j] * (X[, i] - pred)
      }
    }
    g <- colMeans(u)

    # Variance adjustment for estimated reliabilities (paper p. 2041):
    # Sigma = Var(U*_k) where U*_k = U_k - E[dU/dlambda] (E[dV/dlambda])^{-1} V_k
    #
    # All matrices below use the paper's sign conventions directly:
    #   dudlambda = E[dU/dlambda]           (d*p) x d
    #   dvdlambda = E[dV/dlambda]           d x d
    #   V_k       = paper's V_k             n x d

    # E[dU_{(i,j)}/dlambda_s]:
    #   U_{(i,j),k} = I(z_k=z_j)(X_ik - gamma_i - ratio_i * beta_j)
    #   dU/dlambda_s = -I(z_k=z_j) * d(ratio_i * beta_j)/dlambda_s
    #   E[dU/dlambda_s] = -p_z[j] * d(ratio_i * beta_j)/dlambda_s
    # Only nonzero for i > 1 AND j > 1:
    #   d(lambda_i/lambda_1 * beta_j)/dlambda_1 = -ratio_i * beta_j / lambda_1
    #   d(lambda_i/lambda_1 * beta_j)/dlambda_i = beta_j / lambda_1
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
        if (i > 1 && j > 1) {
          dpred_dlambda[col, 1] <- -ratio[i] * beta[j] / lambda[1]
          dpred_dlambda[col, i] <- beta[j] / lambda[1]
        }
      }
    }

    dudlambda <- -diag(diag_pz) %*% dpred_dlambda  # E[dU/dlambda], (d*p) x d

    # E[dV/dlambda] where V_{ik} = sum_{j!=i} lambda_j [(X_i-Xbar)(X_j-Xbar) - lambda_i lambda_j]
    #   dV_i/dlambda_i = sum_{j!=i} lambda_j * (-lambda_j) = -sum_{j!=i} lambda_j^2
    #   dV_i/dlambda_m = (X_i-Xbar)(X_m-Xbar) - 2*lambda_i*lambda_m  (m != i)
    #   E[dV_i/dlambda_m] = mean_obs[...] = mean(u_im) - lambda_i*lambda_m
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

    # Adjusted moment conditions: U*_k = U_k - E[dU/dlambda] (E[dV/dlambda])^{-1} V_k
    uadjust <- u - t(dudlambda %*% solve(dvdlambda) %*% t(V_k))

    # Distance metric statistic with adjusted variance (1/N normalization, matching paper)
    Sigma <- var(uadjust) * ((n - 1) / n)
    n * drop(t(g) %*% solve(Sigma) %*% g)
  }

  # Run optimization
  res <- nlm(q_t0, p = theta_init)

  if (verbose) {
    message("T0 optimization converged with code: ", res$code)
    message("Minimum: ", res$minimum)
  }

  # Degrees of freedom
  df <- (d - 1L) * (p - 1L)

  # Extract estimates
  gamma_hat <- res$estimate[seq_len(d)]
  beta_hat <- c(0, res$estimate[(d + 1):(d + p - 1)])
  names(gamma_hat) <- if (!is.null(colnames(inp$X))) colnames(inp$X) else paste0("X", seq_len(d))
  names(beta_hat) <- paste0("z=", z_levels)

  result <- list(
    statistic = c(T0 = res$minimum),
    parameter = c(df = df),
    p.value = 1 - pchisq(res$minimum, df = df),
    method = "T0: Reliability-dependent test of structural interpretation (VanderWeele & Vansteelandt, 2022)",
    data.name = dname,
    estimates = list(gamma = gamma_hat, beta = beta_hat, lambda = lambda),
    n_obs = n,
    d = d,
    p = p,
    convergence = res$code,
    optim_details = res
  )
  class(result) <- c("structest_t0", "structest", "htest")
  result
}
