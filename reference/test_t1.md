# T1 Test: Reliability-Independent Test of Structural Interpretation (Section 3.3)

Tests whether the structural interpretation of a univariate latent
factor model can be rejected, without requiring estimation of
reliability coefficients. Uses centered indicators and a two-step
efficient GMM procedure.

## Usage

``` r
test_t1(X, z, na.rm = TRUE, max_iter = 1000L, tol = 1e-25, verbose = FALSE)
```

## Arguments

- X:

  numeric matrix (n x d) of indicator variables, d \>= 3.

- z:

  numeric vector of length n encoding the auxiliary variable. Must have
  at least 3 distinct levels.

- na.rm:

  logical; if `TRUE`, rows with any `NA` are removed.

- max_iter:

  integer; maximum iterations for the alternating LS procedure.

- tol:

  numeric; convergence tolerance for the alternating LS procedure.

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

  list with `alpha` and `beta`.

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

The indicators are centered by their row mean: \\X^c\_{ik} = X\_{ik} -
\bar{X}\_{i\cdot}\\, which eliminates the latent factor \\Y\\. Under the
structural interpretation, \\E\[X^c_i \| Z = z_j\] = \alpha_i \beta_j\\
where \\\alpha_i\\ and \\\beta_j\\ are identifiable nuisance parameters.
The test has \\(d-1)(p-2)\\ degrees of freedom. The full model (Equation
3 in the paper) has \\d \times p\\ moment conditions and \\2d + p - 2\\
free parameters (\\d\\ intercepts, \\d\\ alphas, \\p-1\\ betas, minus 1
for the alpha-beta scale indeterminacy), giving \\dp - (2d+p-2) =
(d-1)(p-2)\\ degrees of freedom.

The procedure uses a two-step efficient GMM: first with a fixed weight
matrix from initial estimates, then with the efficient weight matrix.

## References

VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society Series B*, 84, 2032–2054.
