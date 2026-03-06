# Fit Indices and Structural Interpretation

## The problem

A one-factor CFA evaluates whether the **covariance structure** of the
indicators \\X_1, \ldots, X_d\\ is consistent with a single latent
factor. The \\T_1\\ test (VanderWeele & Vansteelandt, 2022) evaluates
whether the **conditional mean structure** \\E(X_i \mid Z = z_j)\\ has
rank 1—a necessary condition for the structural interpretation.

These are fundamentally different questions. This vignette walks through
eight simulated scenarios that show how the two can diverge.

### Why two channels matter

For \\T_1\\ to reject, the mean-difference matrix \\\Delta\_{ij} = E(X_i
\mid Z = z_j) - E(X_i \mid Z = z_1)\\ must have rank \\\> 1\\. If \\Z\\
affects \\X\\ through a **single channel**—whether through \\\eta\\ or
through direct effects—the result is always a single outer product (rank
1):

\\ \Delta\_{ij} = \delta_i \cdot f(z_j) \quad \Rightarrow \quad
\text{rank } 1 \\

For rank \\\> 1\\, \\Z\\ must operate through **two channels** with
different patterns across \\Z\\-levels:

\\ \Delta\_{ij} = \lambda_i \cdot g(z_j) + \delta_i \cdot h(z_j) \\

This is rank 2 when \\g\\ and \\h\\ are not proportional (different
\\Z\\-level ratios) and \\\lambda\\ and \\\delta\\ are not proportional
(different indicator patterns). All rejection scenarios below are
constructed this way.

## Setup

We use \\d = 5\\ indicators throughout. With \\d = 3\\, a one-factor CFA
is just-identified (0 df), so fit indices are trivially perfect. With
\\d = 5\\ the CFA has 5 degrees of freedom—enough for meaningful fit
evaluation.

``` r
library(structest)
library(lavaan)
#> This is lavaan 0.6-21
#> lavaan is FREE software! Please report any bugs.

set.seed(12345)
n <- 2000
d <- 5
mu     <- c(2, 3, 1, 4, 2.5)
lambda <- c(1.0, 0.9, 0.8, 0.7, 0.6)
sd_eps <- 0.5
```

We define two helper functions: `make_X()` generates indicators from a
one-factor model, and `cfa_fit()` fits a one-factor CFA and returns
standard fit indices.

``` r
make_X <- function(eta) {
  X <- sapply(seq_len(d), function(i) {
    mu[i] + lambda[i] * eta + rnorm(n, sd = sd_eps)
  })
  colnames(X) <- paste0("X", seq_len(d))
  X
}

cfa_fit <- function(X) {
  mod <- paste0("F =~ ", paste(colnames(X), collapse = " + "))
  fit <- cfa(mod, data = as.data.frame(X))
  fi <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
  data.frame(
    chisq = fi["chisq"], df = fi["df"], p = fi["pvalue"],
    CFI = fi["cfi"], TLI = fi["tli"], RMSEA = fi["rmsea"], SRMR = fi["srmr"],
    row.names = NULL
  )
}
```

## Scenario 1: Baseline — structural interpretation holds

\\Z \to \eta \to X\\. The standard structural model. \\Z\\ affects the
indicators only through the latent factor.

``` r
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.4 * (z == 1) + 0.9 * (z == 2) + rnorm(n)
X <- make_X(eta)

cfa_fit(X)
#>      chisq df       p       CFI       TLI      RMSEA        SRMR
#> 1 8.973646  5 0.11012 0.9995067 0.9990134 0.01993401 0.004916313
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 2.251, df = 4, p-value = 0.6898
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

CFA fits well. \\T_1\\ does not reject. Everything is consistent.

## Scenario 2: Each indicator affected by \\Z\\ with different effect sizes

\\Z\\ operates through \\\eta\\ (channel 1) **and** directly on each
\\X_i\\ with non-proportional loadings (channel 2). The two channels
have different \\Z\\-level ratios (\\0.3/0.7 \approx 0.43\\ for the
\\\eta\\ channel vs \\0.6/0.2 = 3.0\\ for the direct channel), producing
a rank-2 mean structure.

``` r
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
X <- make_X(eta)
delta <- c(0.3, 0.5, 0.1, 0.4, 0.2)  # not proportional to lambda
X <- X + outer(0.6 * (z == 1) + 0.2 * (z == 2), delta)

cfa_fit(X)
#>      chisq df         p CFI      TLI RMSEA        SRMR
#> 1 2.488976  5 0.7781543   1 1.000634     0 0.002414124
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 14.83, df = 4, p-value = 0.00506
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

CFA is excellent. \\T_1\\ rejects decisively.

## Scenario 3: Only \\X_1\\ directly affected by \\Z\\

\\Z\\ affects \\\eta\\ (all indicators shift proportionally), but
\\X_1\\ also gets an extra direct effect from \\Z\\ with a different
\\Z\\-level pattern.

``` r
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
X <- make_X(eta)
X[, 1] <- X[, 1] + 0.6 * (z == 1) + 0.2 * (z == 2)

cfa_fit(X)
#>      chisq df         p CFI      TLI RMSEA        SRMR
#> 1 1.709306  5 0.8877252   1 1.000867     0 0.002345261
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 153.5, df = 4, p-value = < 2.2e-16
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

## Scenario 4: \\Z\\ affects \\X_1\\, \\X_2\\ directly, bypassing \\\eta\\

\\Z\\ operates through \\\eta\\ and has additional direct effects on
\\X_1\\ and \\X_2\\. The direct effects have a different \\Z\\-level
pattern than the \\\eta\\ channel.

``` r
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
X <- make_X(eta)
X[, 1] <- X[, 1] + 0.5 * (z == 1) + 0.1 * (z == 2)
X[, 2] <- X[, 2] + 0.1 * (z == 1) + 0.5 * (z == 2)

cfa_fit(X)
#>      chisq df       p CFI      TLI RMSEA        SRMR
#> 1 3.122153  5 0.68116   1 1.000487     0 0.002663649
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 200.7, df = 4, p-value = < 2.2e-16
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

## Scenario 5: Confounded, but structural interpretation holds

A confounder \\C\\ causes both \\Z\\ and \\\eta\\. But the indicators
still depend on \\Z\\ only through \\\eta\\—there are no direct \\Z \to
X_i\\ paths. \\E(X_i \mid Z) = \gamma_i + \lambda_i E(\eta \mid Z)\\
remains rank 1.

``` r
C <- rnorm(n)
z_latent <- 0.8 * C + rnorm(n)
z <- as.integer(cut(z_latent,
  breaks = quantile(z_latent, c(0, 1/3, 2/3, 1)),
  include.lowest = TRUE
)) - 1L
eta <- 1 + 0.6 * C + rnorm(n)
X <- make_X(eta)

cfa_fit(X)
#>      chisq df         p       CFI       TLI      RMSEA        SRMR
#> 1 6.044587  5 0.3019045 0.9998885 0.9997769 0.01022051 0.003399797
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 2.538, df = 4, p-value = 0.6378
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

Confounding does not break the structural interpretation—\\Z\\ still
reaches \\X\\ only through \\\eta\\.

## Scenario 6: \\T_1\\ blind spot — single direct channel

\\Z\\ affects only \\X_1\\ directly. There is no \\Z \to \eta\\ path.
The structural interpretation is clearly wrong (there is no latent
factor mediating \\Z \to X\\), yet the mean structure is rank 1 because
there is only one channel:

\\ E(X_i \mid Z = z_j) - E(X_i \mid Z = 0) = \delta_i \cdot f(z_j),
\quad \delta = (c, 0, 0, 0, 0) \\

One channel = one outer product = rank 1. \\T_1\\ cannot detect this.

``` r
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + rnorm(n)  # independent of Z
X <- make_X(eta)
X[, 1] <- X[, 1] + 0.6 * (z == 1) + 1.2 * (z == 2)

cfa_fit(X)
#>      chisq df         p       CFI       TLI      RMSEA        SRMR
#> 1 8.598378  5 0.1261963 0.9994278 0.9988556 0.01896939 0.005932542
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 4.038, df = 4, p-value = 0.4009
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

This is a fundamental limitation: \\T_1\\ is a **necessary condition**
test. Rejection proves the structural interpretation is wrong, but
failure to reject does not prove it is right. Using multiple auxiliary
variables \\Z\\ from different domains can help mitigate this blind
spot.

## Scenario 7: Correlated residuals — CFA fails, structural holds

\\Z \to \eta \to X\\ is the true model (structural interpretation
valid). But \\X_1\\ and \\X_2\\ share method variance—e.g., similar item
wording or the same response format. The 1-factor CFA cannot explain the
extra residual correlation, so fit indices degrade. Yet the conditional
mean structure \\E(X_i \mid Z)\\ is still rank 1 because the method
factor is independent of \\Z\\.

``` r
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.4 * (z == 1) + 0.9 * (z == 2) + rnorm(n)
method <- rnorm(n, sd = 0.7)
X <- make_X(eta)
X[, 1] <- X[, 1] + method
X[, 2] <- X[, 2] + method

cfa_fit(X)
#>     chisq df p       CFI       TLI     RMSEA       SRMR
#> 1 607.184  5 0 0.9198378 0.8396756 0.2453944 0.05210794
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 3.335, df = 4, p-value = 0.5033
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

CFA fit is poor. A researcher might reject this measurement model. But
\\T_1\\ does not reject—the structural interpretation is fine. The poor
CFA fit reflects a covariance misspecification (correlated residuals),
not a failure of the structural interpretation.

## Scenario 8: Two correlated factors masquerading as one

Two distinct latent factors (\\\eta_1\\ and \\\eta_2\\, e.g., anxiety
and depression) are correlated at \\r = 0.95\\. A 1-factor CFA fits
borderline acceptably because the factors are hard to distinguish from
covariance alone. But a treatment \\Z\\ affects \\\eta_1\\ and
\\\eta_2\\ with different patterns across \\Z\\-levels, producing a
rank-2 mean structure.

``` r
z <- sample(0:2, n, replace = TRUE)
eta1 <- 0.5 * (z == 1) + 1.0 * (z == 2) + rnorm(n)
eta2 <- 0.95 * eta1 + sqrt(1 - 0.95^2) * rnorm(n)
eta2 <- eta2 + 0.6 * (z == 1) + 0.1 * (z == 2)
X <- cbind(
  mu[1] + 0.9 * eta1 + rnorm(n, sd = sd_eps),
  mu[2] + 0.8 * eta1 + rnorm(n, sd = sd_eps),
  mu[3] + 0.7 * eta1 + rnorm(n, sd = sd_eps),
  mu[4] + 0.8 * eta2 + rnorm(n, sd = sd_eps),
  mu[5] + 0.7 * eta2 + rnorm(n, sd = sd_eps)
)
colnames(X) <- paste0("X", 1:5)

cfa_fit(X)
#>   chisq df p       CFI       TLI     RMSEA       SRMR
#> 1 123.9  5 0 0.9838097 0.9676195 0.1090413 0.01997892
test_t1(X, z)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 224.1, df = 4, p-value = < 2.2e-16
#> n = 2000, d = 5 indicators, p = 3 Z-levels
```

The CFA chi-square rejects (large \\n\\ makes it sensitive), but CFI and
TLI are in the borderline acceptable range. A researcher might keep the
1-factor model. \\T_1\\ rejects clearly—the indicators do not move
together across \\Z\\-levels in a way consistent with a single
structural factor.

To see why, the ratio of mean shifts across \\Z\\-levels differs by
indicator group:

|                    | \\Z = 1\\  | \\Z = 2\\ | Ratio    |
|--------------------|------------|-----------|----------|
| \\\eta_1\\ channel | \\+0.5\\   | \\+1.0\\  | \\0.50\\ |
| \\\eta_2\\ channel | \\+1.075\\ | \\+1.05\\ | \\1.02\\ |

\\X_1\\–\\X_3\\ (loading on \\\eta_1\\) double from \\Z = 1\\ to \\Z =
2\\. \\X_4\\–\\X_5\\ (loading on \\\eta_2\\) stay roughly flat. The
indicators respond to \\Z\\ in inconsistent ratios—not a single factor.

## References

VanderWeele, T. J., & Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
*84*(5), 2032–2054. <https://doi.org/10.1111/rssb.12555>
