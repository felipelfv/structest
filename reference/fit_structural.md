# Fit a Structural Model with Standard Errors

Fits the structural model from VanderWeele & Vansteelandt (2022) and
returns a coefficient table with point estimates, GMM sandwich standard
errors, z-values, p-values, and confidence intervals. The
overidentification test statistic (T0 or T1) is included as a model
diagnostic.

## Usage

``` r
fit_structural(
  X,
  z,
  type = c("t1", "t0"),
  conf_level = 0.95,
  na.rm = TRUE,
  max_iter = 1000L,
  verbose = FALSE
)
```

## Arguments

- X:

  numeric matrix (n x d) of indicator variables. For `type = "t1"`, d
  \>= 2; for `type = "t0"`, d \>= 3.

- z:

  numeric vector or factor of length n encoding a discrete auxiliary
  variable. For `type = "t1"`, p \>= 3 levels; for `type = "t0"`, p \>=
  2 levels.

- type:

  character; `"t1"` for the reliability-independent model or `"t0"` for
  the reliability-dependent model.

- conf_level:

  numeric; confidence level for confidence intervals (default 0.95).

- na.rm:

  logical; if `TRUE`, rows with any `NA` are removed.

- max_iter:

  integer; maximum iterations for `nlm`.

- verbose:

  logical; if `TRUE`, print progress information.

## Value

An object of class `c("structest_fit_t1", "structest_fit")` or
`c("structest_fit_t0", "structest_fit")` containing:

- coefficients:

  data.frame with columns `Estimate`, `Std.Err`, `z_value`, `p_value`,
  `CI_lower`, `CI_upper`.

- vcov:

  variance-covariance matrix of the parameter estimates.

- estimates:

  list with named vectors of all parameters (including constrained
  ones): `gamma`, `alpha`/`lambda`, `beta`.

- test:

  list with `statistic`, `df`, `p.value` for the overidentification
  test.

- type:

  character; `"t1"` or `"t0"`.

- n_obs:

  number of observations.

- d:

  number of indicators.

- p:

  number of Z-levels.

- data.name:

  character string describing the data.

- conf_level:

  confidence level used.

- convergence:

  convergence code from `nlm`.

## Details

Standard errors are computed using the GMM sandwich variance estimator:
\$\$\mathrm{Var}(\hat\theta) = \frac{1}{n} (G' \Omega^{-1} G)^{-1}\$\$
where \\G = E\[\partial g / \partial \theta'\]\\ is the analytically
derived Jacobian of moment conditions and \\\Omega = \mathrm{Var}(g_k)\\
is the variance of per-observation moment contributions evaluated at the
final estimates.

For `type = "t0"`, the variance \\\Omega\\ is additionally adjusted for
uncertainty in the reliability estimates, as described in the paper (p.
2041).

## References

VanderWeele, T. J., & Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
*84*(5), 2032–2054.
[doi:10.1111/rssb.12555](https://doi.org/10.1111/rssb.12555)

## See also

[`test_t1`](https://felipelfv.github.io/structest/reference/test_t1.md)
and
[`test_t0`](https://felipelfv.github.io/structest/reference/test_t0.md)
for the overidentification tests alone;
[`estimate_reliability`](https://felipelfv.github.io/structest/reference/estimate_reliability.md)
for reliability estimation.

## Examples

``` r
set.seed(12345)
n <- 1000
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
lambda <- c(1.0, 0.8, 0.6)
X <- cbind(
  2 + lambda[1] * eta + rnorm(n, sd = 0.5),
  3 + lambda[2] * eta + rnorm(n, sd = 0.5),
  1 + lambda[3] * eta + rnorm(n, sd = 0.5)
)

# Reliability-independent fit
fit1 <- fit_structural(X, z, type = "t1")
fit1
#> 
#> Structural Model Fit (T1: Reliability-independent)
#> 
#> data:   X and z 
#> n = 1000, d = 3 indicators, p = 3 Z-levels
#> 
#> Coefficients:
#>           Estimate Std.Err z value  Pr(>|z|)
#> gamma[X1]   2.9461  0.0604 48.7392 < 2.2e-16
#> gamma[X2]   3.7879  0.0502 75.4901 < 2.2e-16
#> gamma[X3]   1.6110  0.0403 39.9312 < 2.2e-16
#> alpha[X2]   0.7226  0.0624 11.5740 < 2.2e-16
#> alpha[X3]   0.5236  0.0593  8.8268 < 2.2e-16
#> beta[z=1]   0.2918  0.0845  3.4526 0.0005552
#> beta[z=2]   0.7396  0.0865  8.5472 < 2.2e-16
#> 
#> Overidentification test: J = 0.4784, df = 2, p = 0.7872
#> 
summary(fit1)
#> 
#> Structural Model Fit (T1: Reliability-independent)
#> 
#> data:   X and z 
#> n = 1000, d = 3 indicators, p = 3 Z-levels
#> 
#> Coefficients:
#>           Estimate Std.Err z value  Pr(>|z|)           95% CI
#> gamma[X1]   2.9461  0.0604 48.7392 < 2.2e-16 [2.8276, 3.0646]
#> gamma[X2]   3.7879  0.0502 75.4901 < 2.2e-16 [3.6896, 3.8863]
#> gamma[X3]   1.6110  0.0403 39.9312 < 2.2e-16 [1.5320, 1.6901]
#> alpha[X2]   0.7226  0.0624 11.5740 < 2.2e-16 [0.6003, 0.8450]
#> alpha[X3]   0.5236  0.0593  8.8268 < 2.2e-16 [0.4073, 0.6398]
#> beta[z=1]   0.2918  0.0845  3.4526 0.0005552 [0.1262, 0.4575]
#> beta[z=2]   0.7396  0.0865  8.5472 < 2.2e-16 [0.5700, 0.9092]
#> 
#> Overidentification test: J = 0.4784, df = 2, p = 0.7872
#> (Fail to reject structural interpretation)
#> 
#> Convergence code: 2
#> 
coef(fit1)
#> gamma[X1] gamma[X2] gamma[X3] alpha[X2] alpha[X3] beta[z=1] beta[z=2] 
#> 2.9461134 3.7879127 1.6110381 0.7226236 0.5235738 0.2918066 0.7395983 
confint(fit1)
#>               2.5 %    97.5 %
#> gamma[X1] 2.8276406 3.0645863
#> gamma[X2] 3.6895664 3.8862590
#> gamma[X3] 1.5319627 1.6901134
#> alpha[X2] 0.6002534 0.8449939
#> alpha[X3] 0.4073157 0.6398319
#> beta[z=1] 0.1261550 0.4574583
#> beta[z=2] 0.5699997 0.9091969

# Reliability-dependent fit
fit0 <- fit_structural(X, z, type = "t0")
summary(fit0)
#> 
#> Structural Model Fit (T0: Reliability-dependent)
#> 
#> data:   X and z 
#> n = 1000, d = 3 indicators, p = 3 Z-levels
#> 
#> Coefficients:
#>           Estimate Std.Err z value  Pr(>|z|)           95% CI
#> gamma[X1]   2.9566  0.0596 49.6135 < 2.2e-16 [2.8398, 3.0734]
#> gamma[X2]   3.7767  0.0482 78.2886 < 2.2e-16 [3.6821, 3.8712]
#> gamma[X3]   1.5998  0.0370 43.2250 < 2.2e-16 [1.5273, 1.6724]
#> beta[z=1]   0.2806  0.0823  3.4098 0.0006501 [0.1193, 0.4419]
#> beta[z=2]   0.7168  0.0839  8.5387 < 2.2e-16 [0.5523, 0.8813]
#> 
#> Overidentification test: J = 1.661, df = 4, p = 0.7978
#> (Fail to reject structural interpretation)
#> 
#> Estimated reliabilities (lambda):
#> [1] 1.0496 0.8207 0.6023
#> 
#> Convergence code: 2
#> 
```
