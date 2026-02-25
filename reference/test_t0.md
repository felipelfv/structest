# T0 Test: Reliability-Dependent Test of Structural Interpretation (Section 3.2)

Tests whether the structural interpretation of a univariate latent
factor model can be rejected, using estimated reliability coefficients.
The variance of the estimating equations is adjusted for uncertainty in
the reliability estimates (see paper, p. 2041).

## Usage

``` r
test_t0(X, z, na.rm = TRUE, max_iter = 1000L, verbose = FALSE)
```

## Arguments

- X:

  numeric matrix (n x d) of indicator variables, d \>= 3.

- z:

  numeric vector or factor of length n encoding a discrete auxiliary
  variable. Must have at least 2 distinct levels. If `z` is continuous,
  discretise it (e.g., into quantile groups) before use.

- na.rm:

  logical; if `TRUE`, rows with any `NA` are removed.

- max_iter:

  integer; maximum iterations for `nlm`.

- verbose:

  logical; if `TRUE`, print progress information.

## Value

An object of class `c("structest_t0", "structest", "htest")` containing:

- statistic:

  the T0 test statistic.

- parameter:

  degrees of freedom, \\(d-1)(p-1)\\.

- p.value:

  p-value from chi-squared distribution.

- method:

  description of the test.

- data.name:

  name of the data objects.

- estimates:

  list with `gamma` (intercepts), `beta` (Z effects), and `lambda`
  (reliabilities).

- n_obs:

  number of observations used.

- d:

  number of indicators.

- p:

  number of Z-levels.

- convergence:

  convergence code from `nlm`.

- optim_details:

  full output from `nlm`.

## Details

Under the structural model, \\E\[X_i \| Z = z_j\] = \gamma_i +
(\lambda_i / \lambda_1) \beta_j\\. This gives \\d \times p\\ moment
conditions with \\d + p - 1\\ free parameters (\\\gamma_1, \ldots,
\gamma_d\\ and \\\beta_2, \ldots, \beta_p\\ with \\\beta_1 = 0\\),
yielding \\(d-1)(p-1)\\ degrees of freedom.

The variance of the estimating equations is adjusted for the uncertainty
in the reliability estimates \\\lambda_i\\ (see paper, p. 2041).
Consistent generalised methods of moments estimators (Newey & McFadden,
1994) are obtained by minimising the resulting distance metric
statistic.

## References

VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
84, 2032–2054.
