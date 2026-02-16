# Estimate Reliability Coefficients (Section 3.1)

Estimates the reliability coefficients (loadings) \\\lambda_i\\ for each
indicator in a univariate latent factor model, using pairwise
covariances and a quasi-Poisson GLM with log link.

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

## Details

Under the structural model \\X_i = \alpha_i + \lambda_i Y +
\varepsilon_i\\, we have \\\mathrm{Cov}(X_i, X_j) = \lambda_i
\lambda_j\\. The log of each pairwise covariance is modelled as a linear
combination of \\\log(\lambda_i)\\, yielding a quasi-Poisson GLM.
