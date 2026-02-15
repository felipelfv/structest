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
