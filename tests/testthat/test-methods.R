# Helper to build a mock structest object
make_structest <- function() {
  obj <- list(
    statistic  = c(T0 = 3.456),
    parameter  = c(df = 2),
    p.value    = 0.1778,
    method     = "T0: Reliability-dependent test of structural interpretation (VanderWeele & Vansteelandt, 2022)",
    data.name  = "X, Y, Z",
    estimates  = list(
      gamma  = c(X1 = 0.5, X2 = 0.3),
      alpha  = c(a1 = 0.1),
      beta   = c(`z=0` = 0, `z=1` = 0.7),
      lambda = c(l1 = 0.8, l2 = 0.9)
    ),
    n_obs       = 200,
    d           = 2,
    p           = 2,
    convergence = 0
  )
  class(obj) <- c("structest_t0", "structest", "htest")
  obj
}

# print.structest

test_that("print.structest produces expected output", {
  obj <- make_structest()
  out <- capture.output(print(obj))
  combined <- paste(out, collapse = "\n")

  expect_match(combined, "T0: Reliability-dependent")
  expect_match(combined, "X, Y, Z")
  expect_match(combined, "statistic")
  expect_match(combined, "df = 2")
  expect_match(combined, "p-value")
  expect_match(combined, "n = 200")
  expect_match(combined, "d = 2 indicators")
  expect_match(combined, "p = 2 Z-levels")
})

test_that("print.structest returns object invisibly", {
  obj <- make_structest()
  out <- capture.output(ret <- print(obj))
  expect_identical(ret, obj)
})

test_that("print respects the digits argument", {
  obj <- make_structest()
  out <- capture.output(print(obj, digits = 2))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "3.5")
})

# summary.structest

test_that("summary.structest produces expected output", {
  obj <- make_structest()
  out <- capture.output(summary(obj))
  combined <- paste(out, collapse = "\n")

  expect_match(combined, "T0: Reliability-dependent")
  expect_match(combined, "X, Y, Z")
  expect_match(combined, "Test statistic")
  expect_match(combined, "Degrees of freedom:.*2")
  expect_match(combined, "P-value")
  expect_match(combined, "Sample size:.*200")
  expect_match(combined, "Indicators \\(d\\):.*2")
  expect_match(combined, "Z-levels \\(p\\):.*2")
  expect_match(combined, "Convergence code:.*0")
  expect_match(combined, "gamma \\(intercepts\\)")
  expect_match(combined, "alpha")
  expect_match(combined, "beta \\(Z effects\\)")
  expect_match(combined, "lambda \\(factor loadings\\)")
})

test_that("summary.structest returns object invisibly", {
  obj <- make_structest()
  out <- capture.output(ret <- summary(obj))
  expect_identical(ret, obj)
})

test_that("summary.structest handles missing estimate components", {
  obj <- make_structest()
  obj$estimates$alpha  <- NULL
  obj$estimates$lambda <- NULL

  out <- capture.output(summary(obj))
  combined <- paste(out, collapse = "\n")

  expect_match(combined, "gamma \\(intercepts\\)")
  expect_match(combined, "beta \\(Z effects\\)")
  expect_false(grepl("alpha", combined))
  expect_false(grepl("lambda", combined))
})
