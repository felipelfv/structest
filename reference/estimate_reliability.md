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

## References

VanderWeele, T. J., & Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
*84*(5), 2032–2054.
[doi:10.1111/rssb.12555](https://doi.org/10.1111/rssb.12555)

## See also

[`test_t0`](https://felipelfv.github.io/structest/reference/test_t0.md)
which uses these estimates internally;
[`test_t1`](https://felipelfv.github.io/structest/reference/test_t1.md)
for a reliability-independent alternative.

## Examples

``` r
# Simulate data from a one-factor model
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
