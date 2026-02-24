# Testing the Structural Interpretation of a Latent Factor Model

## Introduction

Factor analysis is one of the most widely used methods in the social and
behavioural sciences. When a single latent factor fits the covariance
among a set of observed indicators, researchers often treat the latent
variable as the causally relevant quantity in subsequent analyses. This
practice implicitly assumes a **structural** interpretation of the
latent factor model: that the latent variable itself, rather than the
individual indicators, is what is causally efficacious.

But a measurement model that fits well does *not* guarantee that this
structural interpretation is correct. The indicators could each have
their own direct causal effects on an outcome, even though their
covariance structure is perfectly consistent with a single latent
factor. If the structural interpretation is wrong, causal conclusions
treating the latent variable as the cause may be misleading.

VanderWeele and Vansteelandt (2022) showed that the structural
interpretation imposes testable empirical constraints. The **structest**
package implements two GMM-based tests—\\T_0\\ and \\T_1\\—that can
reject these constraints, providing researchers with a principled way to
evaluate whether the structural interpretation is tenable.

``` r
library(structest)
```

## The Model

Consider a reflective latent factor model with \\d\\ indicators:

\\ X_i = \mu_i + \lambda_i\\\eta + \varepsilon_i, \qquad i = 1, \ldots,
d \\

where \\\eta\\ is the latent variable, \\\lambda_i\\ are factor loadings
(reliabilities), \\\mu_i\\ are intercepts, and the \\\varepsilon_i\\ are
independent error terms.

Under the **structural** interpretation, the latent variable \\\eta\\ is
the causally relevant quantity. If we have an external variable \\Z\\
that affects \\\eta\\ (and the indicators only *through* \\\eta\\), then
by Theorem 1 of VanderWeele and Vansteelandt (2022), for any indicators
\\i\\, \\j\\ and any values \\z\\, \\z^\*\\:

\\ \lambda_i \bigl\\E(X_j \mid Z=z) - E(X_j \mid Z=z^\*)\bigr\\ =
\lambda_j \bigl\\E(X_i \mid Z=z) - E(X_i \mid Z=z^\*)\bigr\\ \\

This identity is testable with observed data. The \\T_0\\ and \\T_1\\
tests formalise it as overidentifying restrictions in a GMM framework.

## The \\T_0\\ Test (Reliability-Dependent)

The \\T_0\\ test (Section 3.2 of the paper) models the conditional
expectations as:

\\ E(X_i \mid Z = z_j) = \gamma_i + \frac{\lambda_i}{\lambda_1}\\\beta_j
\\

with \\\beta_1 = 0\\ (reference level). This gives \\d \times p\\ moment
conditions with \\d + p - 1\\ free parameters, yielding \\(d-1)(p-1)\\
degrees of freedom. The test requires estimated reliability coefficients
\\\lambda_i\\, which are obtained from pairwise covariances via a
quasi-Poisson GLM.

**When to use:** \\d \geq 3\\ indicators and \\p \geq 2\\ levels of
\\Z\\.

### Example: data consistent with the structural model

``` r
set.seed(12345)
n <- 5000
d <- 5

# Simulate under the structural model
z <- rbinom(n, 1, 0.5)                        # binary auxiliary variable
eta <- 1 + 0.5 * z + rnorm(n)                 # latent factor affected by Z
lambda <- c(1.0, 0.8, 0.6, 0.9, 0.7)         # factor loadings
X <- matrix(NA, n, d)
for (i in 1:d) {
  X[, i] <- 2 + lambda[i] * eta + rnorm(n, sd = 0.5)
}

# T0 test (should NOT reject)
result_t0 <- test_t0(X, z)
print(result_t0)
#> 
#>   T0: Reliability-dependent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 6.075, df = 4, p-value = 0.1936
#> n = 5000, d = 5 indicators, p = 2 Z-levels
```

Because the data were generated under the structural model, the test
statistic is small and the \\p\\-value is large—we fail to reject
\\H_0\\.

### Example: data violating the structural model

Now consider data where \\Z\\ affects each indicator disproportionately,
which violates the structural interpretation:

``` r
set.seed(12345)
eta <- rnorm(n)
X_bad <- cbind(
  2 + 1.0 * eta + 0.5 * z + rnorm(n, sd = 0.3),
  3 + 0.8 * eta + 2.0 * z + rnorm(n, sd = 0.3),   # disproportionate Z effect
  1 + 0.6 * eta - 0.3 * z + rnorm(n, sd = 0.3),
  2 + 0.9 * eta + 0.4 * z + rnorm(n, sd = 0.3),
  1 + 0.7 * eta + 0.8 * z + rnorm(n, sd = 0.3)
)

# T0 test (should reject)
result_t0_bad <- test_t0(X_bad, z)
print(result_t0_bad)
#> 
#>   T0: Reliability-dependent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X_bad and z 
#> statistic = 2226, df = 4, p-value = < 2.2e-16
#> n = 5000, d = 5 indicators, p = 2 Z-levels
```

The test statistic is large and the \\p\\-value is very small, correctly
rejecting the structural interpretation.

We can inspect the full output with
[`summary()`](https://rdrr.io/r/base/summary.html):

``` r
summary(result_t0_bad)
#> 
#>   T0: Reliability-dependent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X_bad and z 
#> 
#> Test statistic:  2226 
#> Degrees of freedom:  4 
#> P-value:  < 2.2e-16 
#> 
#> Sample size:  5000 
#> Indicators (d):  5 
#> Z-levels (p):  2 
#> Convergence code:  1 
#> 
#> Parameter estimates:
#>   gamma (intercepts):
#>     X1     X2     X3     X4     X5 
#> 1.7006 2.7836 0.8274 1.7447 0.7968 
#>   beta (Z effects):
#>   z=0   z=1 
#> 0.000 1.728 
#>   lambda (reliabilities):
#> [1] 1.0099 1.0213 0.4555 0.9002 0.8222
```

## The \\T_1\\ Test (Reliability-Independent)

The \\T_1\\ test (Section 3.3 of the paper) avoids estimating
reliabilities entirely. Under the structural model, the conditional
expectations satisfy (Equation 3 in the paper):

\\ E(X_i \mid Z = z_j) = \gamma_i + \alpha_i\\\beta_j \\

with \\\alpha_1 = 1\\ for identification and \\\beta_1 = 0\\ for the
reference level. This means the matrix of mean differences

\\ \Delta\_{ij} = E(X_i \mid Z = z_j) - E(X_i \mid Z = z_1) \\

has rank \\\leq 1\\ under the structural model. The \\T_1\\ test checks
this rank constraint using a two-step efficient GMM procedure, with \\d
\times p\\ moment conditions, \\2d + p - 2\\ free parameters, and
\\(d-1)(p-2)\\ degrees of freedom.

**When to use:** \\d \geq 3\\ indicators and \\p \geq 3\\ levels of
\\Z\\. Preferred over \\T_0\\ because it does not rely on
error-structure assumptions needed for reliability estimation.

### Example: data consistent with the structural model

``` r
set.seed(12345)
n <- 5000
d <- 4

# Simulate under the structural model with a 4-level Z
z <- sample(0:3, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + 1.2 * (z == 3) + rnorm(n)
lambda <- c(1.0, 0.8, 0.6, 0.9)
X <- matrix(NA, n, d)
for (i in 1:d) {
  X[, i] <- 2 + lambda[i] * eta + rnorm(n, sd = 0.5)
}

# T1 test (should NOT reject)
result_t1 <- test_t1(X, z)
print(result_t1)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X and z 
#> statistic = 3.017, df = 6, p-value = 0.8068
#> n = 5000, d = 4 indicators, p = 4 Z-levels
```

Again, the test does not reject—consistent with data generated under the
structural model.

### Example: data violating the structural model

``` r
set.seed(12345)
eta <- rnorm(n)
X_bad <- cbind(
  2 + 1.0 * eta + 0.5 * (z == 1) + 0.3 * (z == 2) + 0.4 * (z == 3)
    + rnorm(n, sd = 0.3),
  3 + 0.8 * eta + 2.0 * (z == 1) + 0.1 * (z == 2) + 3.0 * (z == 3)
    + rnorm(n, sd = 0.3),
  1 + 0.6 * eta - 0.3 * (z == 1) + 0.8 * (z == 2) - 0.5 * (z == 3)
    + rnorm(n, sd = 0.3),
  2 + 0.9 * eta + 0.4 * (z == 1) + 0.6 * (z == 2) + 0.3 * (z == 3)
    + rnorm(n, sd = 0.3)
)

# T1 test (should reject)
result_t1_bad <- test_t1(X_bad, z)
print(result_t1_bad)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X_bad and z 
#> statistic = 1400, df = 6, p-value = < 2.2e-16
#> n = 5000, d = 4 indicators, p = 4 Z-levels
```

The structural interpretation is rejected, as expected.

``` r
summary(result_t1_bad)
#> 
#>   T1: Reliability-independent test of structural interpretation (VanderWeele & Vansteelandt, 2022) 
#> 
#> data:   X_bad and z 
#> 
#> Test statistic:  1400 
#> Degrees of freedom:  6 
#> P-value:  < 2.2e-16 
#> 
#> Sample size:  5000 
#> Indicators (d):  4 
#> Z-levels (p):  4 
#> Convergence code:  1 
#> 
#> Parameter estimates:
#>   gamma (intercepts):
#>     X1     X2     X3     X4 
#> 2.0874 3.0879 1.0657 2.1165 
#>   alpha:
#>      X1      X2      X3      X4 
#>  1.0000  7.7241 -1.5415  0.6062 
#>   beta (Z effects):
#>     z=0     z=1     z=2     z=3 
#>  0.0000  0.2375 -0.0676  0.3791
```

## Estimating Reliabilities

The
[`estimate_reliability()`](https://felipelfv.github.io/structest/reference/estimate_reliability.md)
function estimates the factor loadings \\\lambda_i\\ from pairwise
covariances. Under the model \\X_i = \mu_i + \lambda_i\\\eta +
\varepsilon_i\\, we have \\\text{Cov}(X_i, X_j) =
\lambda_i\\\lambda_j\\. Taking logs gives \\\log \text{Cov}(X_i, X_j) =
\log\lambda_i + \log\lambda_j\\, which is a linear model in the
log-reliabilities. The function fits this as a quasi-Poisson GLM with
log link.

``` r
set.seed(12345)
n <- 5000
d <- 5
eta <- rnorm(n)
lambda_true <- c(1.0, 0.8, 0.6, 0.9, 0.7)
X <- matrix(NA, n, d)
for (i in 1:d) {
  X[, i] <- lambda_true[i] * eta + rnorm(n, sd = 0.5)
}

rel <- estimate_reliability(X)

# Compare estimated vs true
data.frame(
  true = lambda_true,
  estimated = round(rel$lambda, 3)
)
#>   true estimated
#> 1  1.0     0.974
#> 2  0.8     0.808
#> 3  0.6     0.594
#> 4  0.9     0.898
#> 5  0.7     0.695
```

This function is called internally by
[`test_t0()`](https://felipelfv.github.io/structest/reference/test_t0.md),
but you can use it directly if you need the reliability estimates for
other purposes.

## Choosing Between \\T_0\\ and \\T_1\\

|                                                      | \\T_0\\        | \\T_1\\        |
|------------------------------------------------------|----------------|----------------|
| **Relies on reliability estimates**                  | Yes            | No             |
| **Sensitive to error distribution misspecification** | Yes            | No             |
| **Minimum indicators (\\d\\)**                       | 3              | 3              |
| **Minimum \\Z\\-levels (\\p\\)**                     | 2              | 3              |
| **Degrees of freedom**                               | \\(d-1)(p-1)\\ | \\(d-1)(p-2)\\ |

**Guidance:**

- **\\T_1\\ is generally preferred** because it does not depend on the
  error-structure assumptions needed to estimate reliabilities (e.g.,
  independent, homoscedastic errors).
- **\\T_0\\ is necessary** when only 2 levels of \\Z\\ are available,
  since \\T_1\\ requires \\p \geq 3\\.
- When both tests are applicable (\\d \geq 3\\, \\p \geq 3\\), running
  both provides a useful **sensitivity check**. Agreement strengthens
  the conclusion; disagreement may signal misspecification of the error
  distribution (which affects \\T_0\\ but not \\T_1\\).

## Reference

VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
84, 2032–2054.
[doi:10.1111/rssb.12555](https://doi.org/10.1111/rssb.12555)
