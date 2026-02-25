test_that("test_t1 does not reject under structural model", {
  set.seed(111)
  n <- 3000
  z <- sample(0:2, n, replace = TRUE)
  Y <- 1 + 0.5 * (z == 1) + 1.0 * (z == 2) + rnorm(n)
  lambda <- c(1.0, 0.8, 0.6)
  # Use equal intercepts so that centering by row mean yields rank-1 structure
  X <- cbind(
    lambda[1] * Y + rnorm(n, sd = 0.3),
    lambda[2] * Y + rnorm(n, sd = 0.3),
    lambda[3] * Y + rnorm(n, sd = 0.3)
  )

  result <- test_t1(X, z)
  expect_s3_class(result, "htest")
  expect_s3_class(result, "structest")
  expect_s3_class(result, "structest_t1")
  # Should not reject (structural model holds)
  expect_gt(result$p.value, 0.01)
  expect_equal(result$parameter, c(df = 2))  # (3-1)*(3-2)=2
})

test_that("test_t1 rejects when structural model is violated", {
  set.seed(222)
  n <- 2000
  z <- sample(0:2, n, replace = TRUE)
  Y <- rnorm(n)
  # Disproportionate effects of z on different indicators
  X <- cbind(
    2 + 1.0 * Y + 0.5 * (z == 1) + 1.0 * (z == 2) + rnorm(n, sd = 0.3),
    3 + 0.8 * Y + 2.0 * (z == 1) - 0.5 * (z == 2) + rnorm(n, sd = 0.3),
    1 + 0.6 * Y - 0.3 * (z == 1) + 0.1 * (z == 2) + rnorm(n, sd = 0.3)
  )

  result <- test_t1(X, z)
  # Should reject (structural model violated)
  expect_lt(result$p.value, 0.05)
})

test_that("test_t1 returns correct components", {
  set.seed(333)
  n <- 500
  z <- sample(0:2, n, replace = TRUE)
  Y <- rnorm(n)
  X <- cbind(
    Y + rnorm(n, sd = 0.5),
    0.8 * Y + rnorm(n, sd = 0.5),
    0.6 * Y + rnorm(n, sd = 0.5)
  )

  result <- test_t1(X, z)
  expect_true("statistic" %in% names(result))
  expect_true("parameter" %in% names(result))
  expect_true("p.value" %in% names(result))
  expect_true("estimates" %in% names(result))
  expect_true("gamma" %in% names(result$estimates))
  expect_true("alpha" %in% names(result$estimates))
  expect_true("beta" %in% names(result$estimates))
  expect_equal(result$d, 3)
  expect_equal(result$p, 3)
  expect_equal(result$parameter, c(df = 2))  # (3-1)*(3-2)=2
})

test_that("test_t1 works with d = 2 indicators (paper Section 3.3)", {
  set.seed(444)
  n <- 3000
  z <- sample(0:2, n, replace = TRUE)
  Y <- 1 + 0.5 * (z == 1) + 1.0 * (z == 2) + rnorm(n)
  X <- cbind(
    1.0 * Y + rnorm(n, sd = 0.3),
    0.7 * Y + rnorm(n, sd = 0.3)
  )

  result <- test_t1(X, z)
  expect_equal(result$d, 2)
  expect_equal(result$parameter, c(df = 1))  # (2-1)*(3-2)=1
  expect_gt(result$p.value, 0.01)  # structural model holds
})

test_that("test_t1 recovers known parameters under structural model", {
  set.seed(999)
  n <- 10000
  z <- sample(0:3, n, replace = TRUE)  # p = 4 levels

  # True parameters (Equation 3):
  #   E(X_i | Z = z_j) = gamma_i + alpha_i * beta_j
  #   alpha_1 = 1 (identification), beta_1 = 0 (reference)
  gamma_true <- c(2.0, 5.0, -1.0)    # d = 3
  alpha_true <- c(1.0, 1.5, 0.7)      # alpha_1 = 1
  beta_true  <- c(0.0, 0.8, -0.5, 1.2) # beta_1 = 0

  # Generate data: X_ik = gamma_i + alpha_i * beta_{z_k} + noise
  X <- matrix(NA, n, 3)
  beta_k <- beta_true[match(z, sort(unique(z)))]
  for (i in 1:3) {
    X[, i] <- gamma_true[i] + alpha_true[i] * beta_k + rnorm(n, sd = 0.5)
  }

  result <- test_t1(X, z)

  # Should not reject (structural model holds)
  expect_gt(result$p.value, 0.05)

  # Check parameter recovery (n = 10000, errors should be < 0.03)
  expect_equal(unname(result$estimates$gamma), gamma_true, tolerance = 0.03)
  expect_equal(unname(result$estimates$alpha), alpha_true, tolerance = 0.03)
  expect_equal(unname(result$estimates$beta),  beta_true,  tolerance = 0.03)
})

test_that("test_t1 requires d >= 2", {
  set.seed(445)
  n <- 500
  z <- sample(0:2, n, replace = TRUE)
  X <- matrix(rnorm(n), ncol = 1)

  expect_error(test_t1(X, z), "at least 2")
})

test_that("test_t1 df is correct for d=5, p=3", {
  set.seed(555)
  n <- 2000
  z <- sample(0:2, n, replace = TRUE)
  Y <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
  lambda <- c(1.0, 0.85, 1.1, 0.95, 0.75)
  X <- matrix(NA, n, 5)
  for (i in 1:5) X[, i] <- lambda[i] * Y + rnorm(n, sd = 0.5)

  result <- test_t1(X, z)
  expect_equal(result$parameter, c(df = 4))  # (5-1)*(3-2)=4
  expect_gt(result$p.value, 0.01)
})
