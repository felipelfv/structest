test_that("validate_inputs rejects non-matrix X", {
  expect_error(test_t0(1:10, rep(0:1, 5)), "'X' must be a matrix")
})

test_that("validate_inputs rejects mismatched dimensions", {
  X <- matrix(rnorm(30), ncol = 3)
  z <- rep(0:1, 3)  # length 6, but X has 10 rows
  expect_error(test_t0(X, z), "must equal length")
})

test_that("test_t0 requires d >= 3", {
  X <- matrix(rnorm(20), ncol = 2)
  z <- rep(0:1, 5)
  expect_error(test_t0(X, z), "at least 3 indicators")
})

test_that("test_t0 requires p >= 2", {
  X <- matrix(rnorm(30), ncol = 3)
  z <- rep(0, 10)
  expect_error(test_t0(X, z), "at least 2 distinct Z-levels")
})

test_that("test_t1 requires d >= 2", {
  X <- matrix(rnorm(10), ncol = 1)
  z <- rep(0:2, length.out = 10)
  expect_error(test_t1(X, z), "at least 2 indicators")
})

test_that("test_t1 requires p >= 3", {
  X <- matrix(rnorm(30), ncol = 3)
  z <- rep(0:1, 5)
  expect_error(test_t1(X, z), "at least 3 distinct Z-levels")
})

test_that("NA handling works", {
  set.seed(12345)
  n <- 100
  z <- rep(0:1, length.out = n)
  Y <- 1 + 0.5 * z + rnorm(n)
  X <- cbind(
    Y + rnorm(n, sd = 0.3),
    0.8 * Y + rnorm(n, sd = 0.3),
    0.6 * Y + rnorm(n, sd = 0.3)
  )
  X[1, 1] <- NA
  X[5, 3] <- NA
  # Should run without error when na.rm = TRUE
  result <- test_t0(X, z)
  expect_s3_class(result, "structest")
  expect_equal(result$n_obs, 98L)
})

test_that("factor z is accepted", {
  set.seed(12345)
  n <- 500
  z_fac <- factor(rep(c("A", "B"), length.out = n))
  z_num <- as.integer(z_fac) - 1L
  Y <- 1 + 0.5 * z_num + rnorm(n)
  X <- cbind(
    Y + rnorm(n, sd = 0.3),
    0.8 * Y + rnorm(n, sd = 0.3),
    0.6 * Y + rnorm(n, sd = 0.3)
  )
  result <- test_t0(X, z_fac)
  expect_s3_class(result, "structest")
})
