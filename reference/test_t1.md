# T1 Test: Reliability-Independent Test of Structural Interpretation (Section 3.3)

Tests whether the structural interpretation of a univariate latent
factor model can be rejected, without requiring estimation of
reliability coefficients.

## Usage

``` r
test_t1(X, z, na.rm = TRUE, max_iter = 1000L, tol = 1e-25, verbose = FALSE)
```

## Arguments

- X:

  numeric matrix (n x d) of indicator variables, d \>= 2.

- z:

  numeric vector of length n encoding the auxiliary variable. Must have
  at least 3 distinct levels.

- na.rm:

  logical; if `TRUE`, rows with any `NA` are removed.

- max_iter:

  integer; maximum iterations for `nlm`.

- tol:

  numeric; convergence tolerance for the alternating LS initializer.

- verbose:

  logical; if `TRUE`, print progress information.

## Value

An object of class `c("structest_t1", "structest", "htest")` containing:

- statistic:

  the T1 test statistic.

- parameter:

  degrees of freedom, \\(d-1)(p-2)\\.

- p.value:

  p-value from chi-squared distribution.

- method:

  description of the test.

- data.name:

  name of the data objects.

- estimates:

  list with `gamma` (\\\gamma_i\\), `alpha` (\\\alpha_i\\, with
  \\\alpha_1 = 1\\), and `beta` (\\\beta_j\\, with \\\beta_1 = 0\\).

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

Under the structural model \\X_i = \mu_i + \lambda_i \eta +
\varepsilon_i\\, the conditional expectations satisfy \$\$E(X_i \mid Z =
z_j) = \gamma_i + \alpha_i \beta_j\$\$ where \\\gamma_i\\ are intercepts
absorbing the reference-level means, \\\alpha_i\\ are parameters (with
\\\alpha_1 = 1\\ for identification), and \\\beta_j\\ are parameters
(with \\\beta_1 = 0\\ for the reference level).

This gives \\d \times p\\ moment conditions: \$\$E\[I(Z = z_j)(X_i -
\gamma_i - \alpha_i \beta_j)\] = 0\$\$ with \\2d + p - 2\\ free
parameters (\\d\\ intercepts, \\d - 1\\ alphas, \\p - 1\\ betas),
yielding \\dp - (2d + p - 2) = (d-1)(p-2)\\ degrees of freedom.

The test checks whether the mean-difference matrix \\\Delta\_{ij} =
E(X_i \mid Z = z_j) - E(X_i \mid Z = z_1)\\ has rank \\\le 1\\, which is
the testable implication of the structural model (Theorem 2 in the
paper).

Consistent generalised methods of moments estimators (Newey & McFadden,
1994) are obtained by minimising a distance metric statistic. The
procedure uses two-step estimation to obtain stable initial estimates,
then minimises the criterion with the weight matrix recomputed at each
parameter value. The test statistic is asymptotically
\\\chi^2\_{(d-1)(p-2)}\\ under the null.

## References

VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
84, 2032–2054.
