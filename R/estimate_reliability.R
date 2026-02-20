#' Estimate Factor Loadings (Section 3.1)
#'
#' Estimates the factor loadings \eqn{\lambda_i} for each indicator in a
#' univariate latent factor model, using pairwise covariances and a
#' quasi-Poisson GLM with log link.
#'
#' @param X numeric matrix (n x d) of indicator variables, d >= 3.
#' @param na.rm logical; if \code{TRUE}, rows with any \code{NA} are removed.
#'
#' @return A list with components:
#' \describe{
#'   \item{lambda}{numeric vector of length d with estimated factor loadings.}
#'   \item{glm_fit}{the fitted \code{glm} object.}
#'   \item{V_k}{matrix (n x d) of per-subject estimating function contributions
#'     for the loading parameters (used internally by \code{test_t0}).}
#' }
#'
#' @details
#' Under the structural model \eqn{X_i = \alpha_i + \lambda_i Y + \varepsilon_i},
#' we have \eqn{\mathrm{Cov}(X_i, X_j) = \lambda_i \lambda_j}. The log of
#' each pairwise covariance is modelled as a linear combination of
#' \eqn{\log(\lambda_i)}, yielding a quasi-Poisson GLM.
#'
#' @export
estimate_reliability <- function(X, na.rm = TRUE) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("'X' must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (na.rm) {
    cc <- complete.cases(X)
    if (sum(cc) == 0L) stop("No complete cases remain.", call. = FALSE)
    X <- X[cc, , drop = FALSE]
  }

  n <- nrow(X)
  d <- ncol(X)
  if (d < 3L) {
    stop("Need at least 3 indicators for loading estimation.", call. = FALSE)
  }

  # Compute pairwise sample covariances
  pairs <- pairwise_indices(d)
  n_pairs <- nrow(pairs)
  A <- numeric(n_pairs)
  for (k in seq_len(n_pairs)) {
    A[k] <- cov(X[, pairs[k, 1]], X[, pairs[k, 2]])
  }

  # Check for negative pairwise covariances (paper p. 2040)
  if (any(A <= 0)) {
    neg_pairs <- which(A <= 0)
    pair_labels <- paste0("(X", pairs[neg_pairs, 1], ", X", pairs[neg_pairs, 2], ")")
    warning(
      "Non-positive pairwise covariance(s) detected for pairs: ",
      paste(pair_labels, collapse = ", "), ". ",
      "Loading estimation assumes all pairwise covariances are positive.",
      call. = FALSE
    )
  }

  # Design matrix for quasi-Poisson GLM
  B <- build_design_matrix(d)
  colnames(B) <- paste0("B", seq_len(d))

  # Fit quasi-Poisson GLM with log link (no intercept)
  df <- as.data.frame(B)
  df$A <- A
  fmla <- stats::as.formula(paste("A ~ -1 +", paste(colnames(B), collapse = " + ")))
  fit <- glm(fmla, data = df, family = quasi(link = "log"))
  lambda <- as.numeric(exp(coef(fit)))

  # Per-subject estimating function contributions for lambda
  # u_ij[k] = (X_i[k] - mean_i)(X_j[k] - mean_j) - lambda_i * lambda_j
  # u_i = -sum_{j != i} lambda_j * u_ij  (score contribution for lambda_i)
  means <- colMeans(X)
  Xc <- sweep(X, 2, means)  # centered indicators

  # Matrix of pairwise estimating functions: n x n_pairs
  U_pairs <- matrix(0, nrow = n, ncol = n_pairs)
  for (k in seq_len(n_pairs)) {
    i <- pairs[k, 1]
    j <- pairs[k, 2]
    U_pairs[, k] <- Xc[, i] * Xc[, j] - lambda[i] * lambda[j]
  }

  # Per-subject estimating function contributions for factor loadings (paper p. 2041):
  # V_{ik} = sum_{j != i} lambda_j * {(X_i - Xbar_i)(X_j - Xbar_j) - lambda_i lambda_j}
  V_k <- matrix(0, nrow = n, ncol = d)
  for (i in seq_len(d)) {
    for (k in seq_len(n_pairs)) {
      i1 <- pairs[k, 1]
      i2 <- pairs[k, 2]
      if (i1 == i) {
        V_k[, i] <- V_k[, i] + lambda[i2] * U_pairs[, k]
      } else if (i2 == i) {
        V_k[, i] <- V_k[, i] + lambda[i1] * U_pairs[, k]
      }
    }
  }

  list(lambda = lambda, glm_fit = fit, V_k = V_k,
       U_pairs = U_pairs, pairs = pairs, Xc = Xc)
}
