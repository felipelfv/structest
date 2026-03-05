# Tests for fit_structural()

# Shared simulation helper
sim_structural <- function(n = 1000, seed = 12345) {
  set.seed(seed)
  z <- sample(0:2, n, replace = TRUE)
  eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
  lambda <- c(1.0, 0.8, 0.6)
  X <- cbind(
    2 + lambda[1] * eta + rnorm(n, sd = 0.5),
    3 + lambda[2] * eta + rnorm(n, sd = 0.5),
    1 + lambda[3] * eta + rnorm(n, sd = 0.5)
  )
  list(X = X, z = z, lambda = lambda)
}


test_that("fit_structural T1 returns correct structure", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  expect_s3_class(fit, "structest_fit")
  expect_s3_class(fit, "structest_fit_t1")
  expect_true(is.data.frame(fit$coefficients))
  expect_true(is.matrix(fit$vcov))
  expect_equal(fit$type, "t1")
  expect_equal(fit$d, 3L)
  expect_equal(fit$p, 3L)
  expect_equal(fit$n_obs, 1000L)

  # Coefficient table columns
  expect_named(fit$coefficients,
               c("Estimate", "Std.Err", "z_value", "p_value",
                 "CI_lower", "CI_upper"))

  # Number of free parameters: 2d + p - 2 = 2*3 + 3 - 2 = 7
  expect_equal(nrow(fit$coefficients), 7L)
  expect_equal(nrow(fit$vcov), 7L)
  expect_equal(ncol(fit$vcov), 7L)
})

test_that("fit_structural T1 point estimates match test_t1", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")
  tst <- test_t1(dat$X, dat$z)

  expect_equal(fit$estimates$gamma, tst$estimates$gamma, tolerance = 1e-6)
  expect_equal(fit$estimates$alpha, tst$estimates$alpha, tolerance = 1e-6)
  expect_equal(fit$estimates$beta, tst$estimates$beta, tolerance = 1e-6)
  expect_equal(fit$test$statistic, unname(tst$statistic), tolerance = 1e-6)
})

test_that("fit_structural T1 SEs are positive and finite", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  se <- fit$coefficients$Std.Err
  expect_true(all(se > 0))
  expect_true(all(is.finite(se)))
})

test_that("fit_structural T1 CIs contain estimates", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  expect_true(all(fit$coefficients$CI_lower < fit$coefficients$Estimate))
  expect_true(all(fit$coefficients$CI_upper > fit$coefficients$Estimate))
})

test_that("fit_structural T1 vcov is symmetric positive definite", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  expect_equal(fit$vcov, t(fit$vcov), tolerance = 1e-10)
  evals <- eigen(fit$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))
})


test_that("fit_structural T0 returns correct structure", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t0")

  expect_s3_class(fit, "structest_fit")
  expect_s3_class(fit, "structest_fit_t0")
  expect_equal(fit$type, "t0")

  # Number of free parameters: d + p - 1 = 3 + 3 - 1 = 5
  expect_equal(nrow(fit$coefficients), 5L)
  expect_equal(nrow(fit$vcov), 5L)
})

test_that("fit_structural T0 point estimates match test_t0", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t0")
  tst <- test_t0(dat$X, dat$z)

  expect_equal(fit$estimates$gamma, tst$estimates$gamma, tolerance = 1e-6)
  expect_equal(fit$estimates$beta, tst$estimates$beta, tolerance = 1e-6)
  expect_equal(fit$test$statistic, unname(tst$statistic), tolerance = 1e-6)
})

test_that("fit_structural T0 SEs are positive and finite", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t0")

  se <- fit$coefficients$Std.Err
  expect_true(all(se > 0))
  expect_true(all(is.finite(se)))
})


test_that("coef.structest_fit returns named numeric vector", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  co <- coef(fit)
  expect_true(is.numeric(co))
  expect_equal(length(co), 7L)
})

test_that("confint.structest_fit returns matrix with correct dimensions", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 7L)
  expect_equal(ncol(ci), 2L)

  # Custom level
  ci90 <- confint(fit, level = 0.90)
  # 90% CIs should be narrower than 95% CIs
  width90 <- ci90[, 2] - ci90[, 1]
  width95 <- ci[, 2] - ci[, 1]
  expect_true(all(width90 < width95))
})

test_that("confint.structest_fit works with parm subset", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  # By index
  ci_sub <- confint(fit, parm = 1:3)
  expect_equal(nrow(ci_sub), 3L)

  # By name
  ci_name <- confint(fit, parm = "gamma[X1]")
  expect_equal(nrow(ci_name), 1L)
})

test_that("vcov.structest_fit returns the variance-covariance matrix", {
  dat <- sim_structural()
  fit <- fit_structural(dat$X, dat$z, type = "t1")

  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(V, fit$vcov)
})

test_that("print and summary methods work without error", {
  dat <- sim_structural()
  fit1 <- fit_structural(dat$X, dat$z, type = "t1")
  fit0 <- fit_structural(dat$X, dat$z, type = "t0")

  expect_output(print(fit1), "Structural Model Fit")
  expect_output(print(fit1), "Coefficients:")
  expect_output(print(fit1), "Overidentification test")

  expect_output(summary(fit1), "Structural Model Fit")
  expect_output(summary(fit1), "95%")
  expect_output(summary(fit0), "T0: Reliability-dependent")
  expect_output(summary(fit0), "lambda")
})

test_that("fit_structural rejects invalid conf_level", {
  dat <- sim_structural()
  expect_error(fit_structural(dat$X, dat$z, conf_level = 0),
               "conf_level")
  expect_error(fit_structural(dat$X, dat$z, conf_level = 1),
               "conf_level")
  expect_error(fit_structural(dat$X, dat$z, conf_level = -0.5),
               "conf_level")
})
