# Internal helper functions for structest

#' Validate and prepare inputs for structest functions
#' @param X matrix of indicator variables (n x d)
#' @param z auxiliary variable (length n)
#' @param na.rm logical; remove rows with NA?
#' @param min_d minimum number of indicators required
#' @param min_p minimum number of Z-levels required
#' @return list with cleaned X (matrix), z (vector), n, d, p, z_levels
#' @noRd
validate_inputs <- function(X, z, na.rm = TRUE, min_d = 2L, min_p = 2L) {
  # Coerce X to matrix
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) {
    stop("'X' must be a matrix or data.frame of indicator variables.", call. = FALSE)
  }
  if (!is.numeric(X)) {
    stop("'X' must be numeric.", call. = FALSE)
  }

  if (is.factor(z)) {
    z <- as.integer(z) - 1L
  } else {
    z <- as.vector(z)
    if (!is.numeric(z) && !is.integer(z)) {
      stop("'z' must be numeric, integer, or factor.", call. = FALSE)
    }
  }

  if (nrow(X) != length(z)) {
    stop("Number of rows in 'X' must equal length of 'z'.", call. = FALSE)
  }

  d <- ncol(X)
  if (d < min_d) {
    stop(sprintf("Need at least %d indicators (columns in X), got %d.", min_d, d),
         call. = FALSE)
  }

  # Handle missing data
  if (na.rm) {
    cc <- complete.cases(cbind(X, z))
    if (sum(cc) == 0L) stop("No complete cases remain after removing NAs.", call. = FALSE)
    X <- X[cc, , drop = FALSE]
    z <- z[cc]
  }

  n <- nrow(X)
  z_levels <- sort(unique(z))
  p <- length(z_levels)

  if (p < min_p) {
    stop(sprintf("Need at least %d distinct Z-levels, got %d.", min_p, p),
         call. = FALSE)
  }

  list(X = X, z = z, n = n, d = d, p = p, z_levels = z_levels)
}


#' Generate all choose(d, 2) pairwise index combinations
#' @param d number of indicators
#' @return matrix with 2 columns, each row a pair (i, j) with i < j
#' @noRd
pairwise_indices <- function(d) {
  idx <- combn(d, 2)
  t(idx)
}


#' Build the design matrix B for the quasi-Poisson GLM (loading estimation)
#'
#' For d indicators, there are choose(d,2) pairwise covariances.
#' B is choose(d,2) x d, where B[k, i] = 1 if indicator i is in pair k.
#'
#' @param d number of indicators
#' @return matrix B of dimension choose(d,2) x d
#' @noRd
build_design_matrix <- function(d) {
  pairs <- pairwise_indices(d)
  n_pairs <- nrow(pairs)
  B <- matrix(0L, nrow = n_pairs, ncol = d)
  for (k in seq_len(n_pairs)) {
    B[k, pairs[k, 1]] <- 1L
    B[k, pairs[k, 2]] <- 1L
  }
  B
}


#' Build indicator matrix for Z-levels
#'
#' Returns an n x p matrix where column j is 1 if z == z_levels[j], else 0.
#'
#' @param z numeric vector
#' @param z_levels sorted unique values of z
#' @return matrix of dimension n x p
#' @noRd
z_indicator_matrix <- function(z, z_levels) {
  n <- length(z)
  p <- length(z_levels)
  Z_mat <- matrix(0, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    Z_mat[, j] <- as.numeric(z == z_levels[j])
  }
  Z_mat
}
