# ── CFA Fit vs Structural Interpretation ──────────────────────────────
#
# Demonstrates that good CFA fit does NOT guarantee the structural
# interpretation holds. The T1 test (VanderWeele & Vansteelandt, 2022)
# can detect violations that CFA misses entirely.
#
# Key insight: CFA tests the covariance structure of X (no role for Z).
# T1 tests whether E(X_i | Z=j) has rank-1 structure, which is a
# necessary condition for the structural interpretation.
#
# For T1 to reject, Z must affect X through TWO non-proportional
# channels (e.g., Z -> eta -> X AND Z -> X_i directly, with different
# Z-level patterns). A single channel always gives rank 1.
#
# We use d=5 indicators so CFA has df=5 (with d=3, CFA is saturated
# and fit indices are trivially perfect).
# ─────────────────────────────────────────────────────────────────────

library(structest)
library(lavaan)

set.seed(42)
n <- 2000
d <- 5
mu     <- c(2, 3, 1, 4, 2.5)
lambda <- c(1.0, 0.9, 0.8, 0.7, 0.6)
sd_eps <- 0.5

# ── Helpers ──────────────────────────────────────────────────────────

make_X <- function(eta) {
  X <- sapply(seq_len(d), function(i) mu[i] + lambda[i] * eta + rnorm(n, sd = sd_eps))
  colnames(X) <- paste0("X", seq_len(d))
  X
}

cfa_fit <- function(X) {
  mod <- paste0("F =~ ", paste(colnames(X), collapse = " + "))
  fit <- cfa(mod, data = as.data.frame(X))
  fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
}

run <- function(X, z, label) {
  fi <- cfa_fit(X)
  t1 <- test_t1(X, z)
  tag <- if (t1$p.value < 0.05) "** REJECT **" else "do not reject"
  cat(sprintf(
    paste0("\n%s\n",
           "  CFA:  chisq = %6.2f  df = %g  p = %.4f\n",
           "        CFI = %.3f  TLI = %.3f  RMSEA = %.3f  SRMR = %.3f\n",
           "  T1:   stat = %6.2f  df = %d  p = %.4f  %s\n"),
    label,
    fi["chisq"], fi["df"], fi["pvalue"],
    fi["cfi"], fi["tli"], fi["rmsea"], fi["srmr"],
    t1$statistic, t1$parameter, t1$p.value, tag
  ))
}

cat("===================================================================\n")
cat("  CFA Fit vs Structural Interpretation (T1 Test)\n")
cat("  n =", n, " d =", d, " indicators\n")
cat("===================================================================\n")

# ── 1. Baseline: Z -> eta -> X (structural holds) ───────────────────
# The standard structural model. Z affects X only through eta.
# CFA fits well. T1 does not reject.

z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.4 * (z == 1) + 0.9 * (z == 2) + rnorm(n)
X <- make_X(eta)
run(X, z, "1. Z -> eta -> X  (structural holds)")

# ── 2. Each X_i affected by Z with different effect sizes ───────────
# Z affects eta (channel 1), AND Z affects each X_i directly with
# non-proportional loadings (channel 2). The two channels have
# different Z-level ratios, so the mean structure is rank > 1.
# CFA still fits well (direct effects are small). T1 rejects.

z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
X <- make_X(eta)
# Direct effects: ratio across Z-levels = 0.6/0.2 = 3.0
# (vs eta channel: 0.3/0.7 = 0.43 -- very different)
delta <- c(0.3, 0.5, 0.1, 0.4, 0.2)  # not proportional to lambda
X <- X + outer(0.6 * (z == 1) + 0.2 * (z == 2), delta)
run(X, z, "2. Z -> each X_i with different effects  (reject)")

# ── 3. Only X_1 directly affected by Z ──────────────────────────────
# Z affects eta (all indicators shift proportionally), but X_1 also
# gets an extra direct kick from Z with a different pattern.
# CFA fits well. T1 rejects.

z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
X <- make_X(eta)
X[, 1] <- X[, 1] + 0.6 * (z == 1) + 0.2 * (z == 2)
run(X, z, "3. Z -> X1 only  (reject)")

# ── 4. Z -> X_1, X_2 directly, bypassing eta ────────────────────────
# Z operates through eta AND has direct effects on X_1 and X_2.
# The direct effects have a different Z-level pattern than the eta
# channel, breaking the rank-1 structure.
# CFA fits well. T1 rejects.

z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
X <- make_X(eta)
X[, 1] <- X[, 1] + 0.5 * (z == 1) + 0.1 * (z == 2)
X[, 2] <- X[, 2] + 0.1 * (z == 1) + 0.5 * (z == 2)
run(X, z, "4. Z -> X1, X2 directly, bypassing eta  (reject)")

# ── 5. Confounded but structural interpretation holds ────────────────
# A confounder C causes both Z and eta. But X_i still depend on Z
# only through eta (no direct Z -> X_i paths).
# E(X_i | Z) = gamma_i + lambda_i * E(eta | Z) -- still rank 1.
# CFA fits well. T1 does not reject.

C <- rnorm(n)
z_latent <- 0.8 * C + rnorm(n)
z <- as.integer(cut(z_latent,
  breaks = quantile(z_latent, c(0, 1/3, 2/3, 1)),
  include.lowest = TRUE
)) - 1L
eta <- 1 + 0.6 * C + rnorm(n)
X <- make_X(eta)
run(X, z, "5. Confounded (C -> Z, C -> eta)  (structural holds)")

# ── 6. Z -> X_1 directly, but NOT through eta (T1 blind spot) ───────
# Z affects only X_1 directly. There is NO Z -> eta path.
# The structural interpretation is clearly wrong (Z doesn't act
# through a latent factor), yet the mean structure is rank 1:
#   E(X_i | Z) = gamma_i + delta_i * f(Z),  with delta = (c, 0, 0, 0, 0)
# One channel = one outer product = rank 1. T1 cannot detect this.

z <- sample(0:2, n, replace = TRUE)
eta <- 1 + rnorm(n)  # eta independent of Z
X <- make_X(eta)
X[, 1] <- X[, 1] + 0.6 * (z == 1) + 1.2 * (z == 2)  # only direct path
run(X, z, "6. Z -> X1 only, NO Z -> eta  (T1 blind spot, does NOT reject)")

# ── 7. Correlated residuals: CFA fails, structural holds ────────────
# Z -> eta -> X is the true model (structural interpretation valid).
# But X_1 and X_2 share method variance (correlated residuals),
# e.g., similar item wording or same response format.
# CFA with 1 factor can't explain the extra residual correlation,
# so fit indices degrade. Yet T1 does not reject, because the
# conditional mean structure E(X_i | Z) is still rank 1.

z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.4 * (z == 1) + 0.9 * (z == 2) + rnorm(n)
method <- rnorm(n, sd = 0.7)  # shared method factor for X1, X2
X <- make_X(eta)
X[, 1] <- X[, 1] + method
X[, 2] <- X[, 2] + method
run(X, z, "7. Correlated residuals  (CFA bad, structural holds)")

# ── 8. Two correlated factors masquerading as one ────────────────────
# Two distinct latent factors (e.g., anxiety and depression) are
# correlated at r = 0.85. A 1-factor CFA fits acceptably because
# the factors are hard to distinguish from covariance alone.
# But Z (e.g., a treatment) affects eta_1 and eta_2 with DIFFERENT
# patterns across Z-levels, producing a rank-2 mean structure.
# T1 catches what CFA misses.

z <- sample(0:2, n, replace = TRUE)
eta1 <- 0.5 * (z == 1) + 1.0 * (z == 2) + rnorm(n)
eta2 <- 0.95 * eta1 + sqrt(1 - 0.95^2) * rnorm(n)  # cor(eta1, eta2) ~ 0.95
eta2 <- eta2 + 0.6 * (z == 1) + 0.1 * (z == 2)      # different Z-pattern
X <- cbind(
  mu[1] + 0.9 * eta1 + rnorm(n, sd = sd_eps),
  mu[2] + 0.8 * eta1 + rnorm(n, sd = sd_eps),
  mu[3] + 0.7 * eta1 + rnorm(n, sd = sd_eps),
  mu[4] + 0.8 * eta2 + rnorm(n, sd = sd_eps),
  mu[5] + 0.7 * eta2 + rnorm(n, sd = sd_eps)
)
colnames(X) <- paste0("X", 1:5)
run(X, z, "8. Two correlated factors (r=.95)  (CFA borderline, T1 rejects)")
