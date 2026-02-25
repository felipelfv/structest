test_that("test_t0 does not reject under structural model", {
  set.seed(789)
  n <- 2000
  z <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * z + rnorm(n)
  lambda <- c(1.0, 0.8, 0.6)
  X <- cbind(
    2 + lambda[1] * Y + rnorm(n, sd = 0.5),
    3 + lambda[2] * Y + rnorm(n, sd = 0.5),
    1 + lambda[3] * Y + rnorm(n, sd = 0.5)
  )

  result <- test_t0(X, z)
  expect_s3_class(result, "htest")
  expect_s3_class(result, "structest")
  expect_s3_class(result, "structest_t0")
  # Should not reject at alpha = 0.05 (structural model holds)
  expect_gt(result$p.value, 0.01)
  expect_equal(result$parameter, c(df = 2))  # (3-1)*(2-1)=2
})

test_that("test_t0 rejects when structural model is violated", {
  set.seed(101)
  n <- 2000
  z <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  X <- cbind(
    2 + 1.0 * Y + 0.5 * z + rnorm(n, sd = 0.3),
    3 + 0.8 * Y + 2.0 * z + rnorm(n, sd = 0.3),  # disproportionate z effect
    1 + 0.6 * Y - 0.3 * z + rnorm(n, sd = 0.3)   # opposite z effect
  )

  result <- test_t0(X, z)
  # Should reject (structural model violated)
  expect_lt(result$p.value, 0.05)
})

test_that("test_t0 returns correct components", {
  set.seed(202)
  n <- 500
  z <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  X <- cbind(
    Y + rnorm(n, sd = 0.5),
    0.8 * Y + rnorm(n, sd = 0.5),
    0.6 * Y + rnorm(n, sd = 0.5)
  )

  result <- test_t0(X, z)
  expect_true("statistic" %in% names(result))
  expect_true("parameter" %in% names(result))
  expect_true("p.value" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("estimates" %in% names(result))
  expect_true("lambda" %in% names(result$estimates))
  expect_true("gamma" %in% names(result$estimates))
  expect_true("beta" %in% names(result$estimates))
  expect_equal(result$d, 3)
  expect_equal(result$p, 2)
})

test_that("test_t0 recovers known parameters under structural model", {
  set.seed(888)
  n <- 10000
  z <- sample(0:2, n, replace = TRUE)  # p = 3 levels

  # True parameters (Equation 2):
  #   E(X_i | Z = z_j) = gamma_i + (lambda_i/lambda_1) * beta_j
  #   beta_1 = 0 (reference level)
  lambda_true <- c(1.0, 0.8, 0.6)
  ratio_true  <- lambda_true / lambda_true[1]  # c(1.0, 0.8, 0.6)
  gamma_true  <- c(2.0, 5.0, -1.0)
  beta_true   <- c(0.0, 0.5, 1.2)  # beta_1 = 0

  # Generate data: X_ik = mu_i + lambda_i * eta_k + eps_ik
  # where eta_k ~ N(0, 1), eps independent
  eta <- rnorm(n)
  # Shift eta mean by Z-level to produce the right conditional expectations
  # E(eta | Z = z_j) = beta_j / lambda_1 = beta_j (since lambda_1 = 1)
  eta <- eta + beta_true[match(z, sort(unique(z)))]
  mu <- c(2.0, 5.0 - 0.8 * 0, -1.0 - 0.6 * 0)  # adjust so gamma works out
  X <- cbind(
    gamma_true[1] + ratio_true[1] * beta_true[match(z, sort(unique(z)))] + rnorm(n, sd = 0.5),
    gamma_true[2] + ratio_true[2] * beta_true[match(z, sort(unique(z)))] + rnorm(n, sd = 0.5),
    gamma_true[3] + ratio_true[3] * beta_true[match(z, sort(unique(z)))] + rnorm(n, sd = 0.5)
  )

  result <- test_t0(X, z)

  # Should not reject (structural model holds)
  expect_gt(result$p.value, 0.05)

  # Check parameter recovery (n = 10000, errors should be < 0.03)
  expect_equal(unname(result$estimates$gamma), gamma_true, tolerance = 0.03)
  expect_equal(unname(result$estimates$beta),  beta_true,  tolerance = 0.03)
})

test_that("test_t0 works with p > 2 Z-levels", {
  set.seed(303)
  n <- 1500
  z <- sample(0:2, n, replace = TRUE)
  Y <- 1 + 0.3 * z + rnorm(n)
  X <- cbind(
    Y + rnorm(n, sd = 0.4),
    0.7 * Y + rnorm(n, sd = 0.4),
    0.5 * Y + rnorm(n, sd = 0.4)
  )

  result <- test_t0(X, z)
  expect_equal(result$parameter, c(df = 4))  # (3-1)*(3-1)=4
  expect_gt(result$p.value, 0.01)
})
