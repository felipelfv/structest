# Structural-holding data reused across tests: T1 statistic ~ chi^2_df, so the
# fitted object is a realistic input for equiv_test().
make_fit <- function(seed = 111, n = 2000) {
  set.seed(seed)
  z <- sample(0:2, n, replace = TRUE)
  Y <- 1 + 0.5 * (z == 1) + 1.0 * (z == 2) + rnorm(n)
  lambda <- c(1.0, 0.8, 0.6)
  X <- cbind(
    lambda[1] * Y + rnorm(n, sd = 0.3),
    lambda[2] * Y + rnorm(n, sd = 0.3),
    lambda[3] * Y + rnorm(n, sd = 0.3)
  )
  test_t1(X, z)
}

test_that("equiv_test rejects non-structest input", {
  expect_error(equiv_test(list(statistic = 1), margin = 0.05),
               "structest")
})

test_that("equiv_test requires a margin (no default)", {
  fit <- make_fit()
  expect_error(equiv_test(fit), "margin")
})

test_that("equiv_test validates margin and alpha", {
  fit <- make_fit()
  expect_error(equiv_test(fit, margin = -0.1), "positive")
  expect_error(equiv_test(fit, margin = 0),    "positive")
  expect_error(equiv_test(fit, margin = c(0.05, 0.1)), "single")
  expect_error(equiv_test(fit, margin = 0.05, alpha = 0),  "\\(0, 1\\)")
  expect_error(equiv_test(fit, margin = 0.05, alpha = 1),  "\\(0, 1\\)")
})

test_that("equiv_test returns the documented structure", {
  fit <- make_fit()
  res <- equiv_test(fit, margin = 0.05)
  expect_s3_class(res, "structest_equiv")
  expect_s3_class(res, "structest")
  expect_s3_class(res, "htest")
  for (nm in c("statistic", "parameter", "p.value", "margin", "alpha",
               "epsilon", "epsilon_upper", "equivalence")) {
    expect_true(nm %in% names(res), info = nm)
  }
  expect_type(res$equivalence, "logical")
  expect_identical(unname(res$statistic), unname(fit$statistic))
  expect_identical(unname(res$parameter), unname(fit$parameter))
})

test_that("epsilon point estimate matches sqrt(max(0,(J-df)/(n*df)))", {
  fit <- make_fit()
  res <- equiv_test(fit, margin = 0.05)
  J  <- unname(fit$statistic); df <- unname(fit$parameter); n <- fit$n_obs
  expect_equal(res$epsilon, sqrt(max(0, (J - df) / (n * df))))
})

test_that("decision duality: equivalence <=> epsilon_upper < margin <=> p < alpha", {
  fit <- make_fit()
  # check across several margins so we hit both decisions
  for (m in c(0.01, 0.05, 0.1, 0.3)) {
    res <- equiv_test(fit, margin = m, alpha = 0.05)
    expect_equal(isTRUE(res$equivalence), res$epsilon_upper < res$margin)
    expect_equal(isTRUE(res$equivalence), res$p.value < res$alpha)
  }
})

test_that("a generous margin supports equivalence, a tiny one does not", {
  fit <- make_fit()
  expect_true(equiv_test(fit,  margin = 0.5)$equivalence)
  expect_false(equiv_test(fit, margin = 1e-4)$equivalence)
})

test_that("epsilon_upper is a valid (1-alpha) upper bound on epsilon", {
  fit <- make_fit()
  res <- equiv_test(fit, margin = 0.05, alpha = 0.05)
  expect_gte(res$epsilon_upper, res$epsilon)
  expect_gte(res$epsilon_upper, 0)
})

test_that("equiv_test errors on a just-identified model (df < 1)", {
  fake <- structure(
    list(statistic = c(T1 = 0.5), parameter = c(df = 0L), p.value = 0.5,
         data.name = "X and z", n_obs = 500L, d = 2L, p = 2L),
    class = c("structest_t1", "structest", "htest")
  )
  expect_error(equiv_test(fake, margin = 0.05), "just-identified")
})

test_that("print.structest_equiv returns its argument invisibly", {
  fit <- make_fit()
  res <- equiv_test(fit, margin = 0.05)
  expect_output(print(res), "Equivalence")
  expect_invisible(print(res))
})
