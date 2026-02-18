test_that("estimate_reliability returns correct structure", {
  set.seed(123)
  n <- 500
  Y <- rnorm(n)
  lambda_true <- c(1.0, 0.8, 0.6)
  X <- cbind(
    lambda_true[1] * Y + rnorm(n, sd = 0.5),
    lambda_true[2] * Y + rnorm(n, sd = 0.5),
    lambda_true[3] * Y + rnorm(n, sd = 0.5)
  )

  result <- estimate_reliability(X)
  expect_type(result, "list")
  expect_true("lambda" %in% names(result))
  expect_true("glm_fit" %in% names(result))
  expect_true("V_k" %in% names(result))
  expect_length(result$lambda, 3)
  expect_equal(nrow(result$V_k), n)
  expect_equal(ncol(result$V_k), 3)
})

test_that("estimate_reliability recovers approximate loadings", {
  set.seed(456)
  n <- 5000
  Y <- rnorm(n)
  lambda_true <- c(1.0, 0.8, 0.6, 0.4)
  X <- cbind(
    lambda_true[1] * Y + rnorm(n, sd = 0.3),
    lambda_true[2] * Y + rnorm(n, sd = 0.3),
    lambda_true[3] * Y + rnorm(n, sd = 0.3),
    lambda_true[4] * Y + rnorm(n, sd = 0.3)
  )

  result <- estimate_reliability(X)
  # Loadings are identified up to sign/scale from covariances
  # The ratios should be approximately correct
  ratios_true <- lambda_true / lambda_true[1]
  ratios_est <- result$lambda / result$lambda[1]
  expect_equal(ratios_est, ratios_true, tolerance = 0.15)
})

test_that("estimate_reliability requires d >= 3", {
  X <- matrix(rnorm(20), ncol = 2)
  expect_error(estimate_reliability(X), "at least 3 indicators")
})

test_that("estimate_reliability accepts a data.frame", {
  set.seed(100)
  n <- 200
  Y <- rnorm(n)
  X <- data.frame(
    a = 1.0 * Y + rnorm(n, sd = 0.3),
    b = 0.8 * Y + rnorm(n, sd = 0.3),
    c = 0.6 * Y + rnorm(n, sd = 0.3)
  )
  result <- estimate_reliability(X)
  expect_length(result$lambda, 3)
})

test_that("estimate_reliability rejects non-numeric input", {
  X <- matrix(letters[1:9], ncol = 3)
  expect_error(estimate_reliability(X), "numeric matrix or data.frame")
  expect_error(estimate_reliability("not a matrix"), "numeric matrix or data.frame")
})

test_that("estimate_reliability drops rows with NAs", {
  set.seed(200)
  n <- 300
  Y <- rnorm(n)
  X <- cbind(
    1.0 * Y + rnorm(n, sd = 0.3),
    0.8 * Y + rnorm(n, sd = 0.3),
    0.6 * Y + rnorm(n, sd = 0.3)
  )
  X[1, 1] <- NA
  X[5, 3] <- NA
  result <- estimate_reliability(X, na.rm = TRUE)
  expect_equal(nrow(result$V_k), n - 2)
  expect_length(result$lambda, 3)
})

test_that("estimate_reliability errors when all rows are NA", {
  X <- matrix(NA_real_, nrow = 5, ncol = 3)
  expect_error(estimate_reliability(X), "No complete cases")
})

test_that("estimate_reliability warns on non-positive covariances", {
  set.seed(300)
  n <- 200
  X <- cbind(
    rnorm(n),
    rnorm(n),
    rnorm(n)
  )
  # Force a negative pairwise covariance by making col 2 anti-correlated with col 1
  X[, 2] <- -X[, 1] + rnorm(n, sd = 0.01)
  # The warning fires before the GLM (which may also error), so just check the warning
  expect_warning(
    tryCatch(estimate_reliability(X), error = function(e) NULL),
    "Non-positive pairwise covariance"
  )
})
