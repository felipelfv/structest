# Estimate Reliability Coefficients (Section 3.1)

Estimates the reliability coefficients \\\lambda_i\\ for each indicator
in a univariate latent factor model, using pairwise covariances and a
quasi-Poisson GLM with log link.

## Usage

``` r
estimate_reliability(X, na.rm = TRUE)
```

## Arguments

- X:

  numeric matrix (n x d) of indicator variables, d \>= 3.

- na.rm:

  logical; if `TRUE`, rows with any `NA` are removed.

## Value

A list with components:

- lambda:

  numeric vector of length d with estimated reliabilities.

- glm_fit:

  the fitted `glm` object.

- V_k:

  matrix (n x d) of per-subject estimating function contributions for
  the reliability parameters (used internally by `test_t0`).

- U_pairs:

  matrix (n x choose(d,2)) of pairwise estimating function contributions
  (used internally by `test_t0`).

- pairs:

  matrix (choose(d,2) x 2) of indicator index pairs.

- Xc:

  matrix (n x d) of mean-centred indicators.

## Details

Under the structural model \\X_i = \mu_i + \lambda_i \eta +
\varepsilon_i\\, we have \\\mathrm{Cov}(X_i, X_j) = \lambda_i
\lambda_j\\. The log of each pairwise covariance is modelled as a linear
combination of \\\log(\lambda_i)\\, yielding a quasi-Poisson GLM.

## Examples

``` r
# Simulate data from a one-factor model
set.seed(42)
n <- 1000
eta <- rnorm(n)
lambda_true <- c(1.0, 0.8, 0.6)
X <- cbind(
  lambda_true[1] * eta + rnorm(n, sd = 0.5),
  lambda_true[2] * eta + rnorm(n, sd = 0.5),
  lambda_true[3] * eta + rnorm(n, sd = 0.5)
)

rel <- estimate_reliability(X)
rel$lambda
#> [1] 1.0221319 0.7752837 0.6143667

# Ratios lambda_i / lambda_1 should be close to truth
rel$lambda / rel$lambda[1]
#> [1] 1.0000000 0.7584967 0.6010640

set.seed(12345)
n <- 1000
eta <- rnorm(n)
lambda_true <- c(1.0, 0.8, 0.6)
X <- cbind(
  lambda_true[1] * eta + rnorm(n, sd = 0.5),
  lambda_true[2] * eta + rnorm(n, sd = 0.5),
  lambda_true[3] * eta + rnorm(n, sd = 0.5)
)

rel <- estimate_reliability(X)
rel$lambda
#> [1] 1.0258103 0.7859982 0.5831385

# Ratios lambda_i / lambda_1 should be close to truth
rel$lambda / rel$lambda[1]
#> [1] 1.0000000 0.7662218 0.5684662
```
