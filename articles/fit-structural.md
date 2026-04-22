# Fitting the Structural Model with Standard Errors

## Introduction

The
[`test_t0()`](https://felipelfv.github.io/structest/reference/test_t0.md)
and
[`test_t1()`](https://felipelfv.github.io/structest/reference/test_t1.md)
functions estimate the structural model parameters via GMM and return
both a test statistic and the point estimates (accessible via
`$estimates`). However, they do not provide standard errors for those
estimates—they tell you whether the structural interpretation can be
rejected, but not how precisely each individual parameter is estimated
or whether individual effects are statistically significant.

The
[`fit_structural()`](https://felipelfv.github.io/structest/reference/fit_structural.md)
function fills this gap. It runs the same GMM estimation as the test
functions, then additionally computes sandwich standard errors for every
free parameter. The result is a coefficient table—analogous to
[`summary.lm()`](https://rdrr.io/r/stats/summary.lm.html)—with point
estimates, standard errors, \\z\\-values, \\p\\-values, and confidence
intervals. The overidentification test statistic (equivalent to \\T_0\\
or \\T_1\\) is reported as a model diagnostic.

This vignette assumes familiarity with the model and notation introduced
in the *Getting Started* vignette
([`vignette("structest")`](https://felipelfv.github.io/structest/articles/structest.md)).
We briefly recall the key definitions below before developing the
standard-error machinery.

``` r
library(structest)
```

## Recap: The Structural Model and Moment Conditions

Recall from the *Getting Started* vignette that under the basic latent
factor model (Equation 1 there), each indicator satisfies \\X_i =
\lambda_i\\\eta + \varepsilon_i\\. The structural interpretation
requires that any external discrete variable \\Z\\ affects the
indicators only through \\\eta\\, which imposes testable restrictions on
the conditional means \\E(X_i \mid Z = z_j)\\.

### \\T_1\\ parameterisation (reliability-independent)

Under Equation (3) of the *Getting Started* vignette, the restricted
model is

\\ E(X_i \mid Z = z_j) = \gamma_i + \alpha_i\\\beta_j \\

with identification constraints \\\alpha_1 = 1\\ and \\\beta_1 = 0\\.
The free parameter vector is

\\ \theta = (\gamma_1, \ldots, \gamma_d,\\ \alpha_2, \ldots, \alpha_d,\\
\beta_2, \ldots, \beta_p) \\

with \\2d + p - 2\\ free parameters in total. Here:

- \\\gamma_i = E(X_i \mid Z = z_1)\\: intercept for indicator \\i\\ at
  the reference level of \\Z\\.
- \\\alpha_i = \lambda_i / \lambda_1\\: the relative loading of
  indicator \\i\\ compared to the first indicator (\\\alpha_1 = 1\\ is
  fixed for identification).
- \\\beta_j\\: the shift in the latent factor at \\Z\\-level \\j\\
  relative to the reference (\\\beta_1 = 0\\).

### \\T_0\\ parameterisation (reliability-dependent)

Under Equation (2) of the *Getting Started* vignette, the restricted
model is

\\ E(X_i \mid Z = z_j) = \gamma_i + \frac{\lambda_i}{\lambda_1}\\\beta_j
\\

with \\\beta_1 = 0\\ and the ratios \\\lambda_i / \lambda_1\\
pre-estimated from pairwise covariances. Here the free parameter vector
is

\\ \theta = (\gamma_1, \ldots, \gamma_d,\\ \beta_2, \ldots, \beta_p) \\

with \\d + p - 1\\ free parameters. The loading ratios are treated as
known constants plugged in from the reliability estimation step.

## Moment Conditions

Both parameterisations are estimated by the generalised method of
moments (GMM). Using the same notation as the *Getting Started*
vignette, the moment conditions are defined by the residuals: for each
indicator \\i\\ and \\Z\\-level \\j\\, observation \\k\\ contributes

\\ U\_{k,(i,j)} = I(Z_k = z_j)\\\left( X\_{ik} - E\_\theta(X_i \mid Z =
z_j)\right) \tag{4} \\

where \\E\_\theta(\cdot)\\ is the predicted conditional mean under the
model. For the \\T_1\\ parameterisation, this is \\E\_\theta(X_i \mid Z
= z_j) = \gamma_i + \alpha_i\\\beta_j\\; for \\T_0\\, it is
\\E\_\theta(X_i \mid Z = z_j) = \gamma_i +
(\lambda_i/\lambda_1)\\\beta_j\\.

Stacking all \\d \times p\\ moment conditions gives a vector
\\U_k(\theta) \in \mathbb{R}^{d \times p}\\. Under correct
specification, \\E\[U_k(\theta_0)\] = 0\\ — the residuals average to
zero at the true parameter values.

### From moments to estimation

The GMM estimator finds the \\\hat\theta\\ that makes the sample average
of these moment conditions as close to zero as possible. Specifically,
it minimises the distance metric statistic

\\ Q_N(\theta) = N\\\bar U(\theta)^\top\\\Sigma^{-1}\\\bar U(\theta)
\tag{5} \\

where \\\bar U(\theta) = N^{-1}\sum\_{k=1}^{N} U_k(\theta)\\ is the
sample average and \\\Sigma\\ is the empirical covariance matrix of the
\\U_k\\ (with the variance adjustment for \\T_0\\ described in the
*Getting Started* vignette). The weight matrix \\\Sigma^{-1}\\
standardises each moment condition by its variance and accounts for
correlations between them.

### From estimation to testing

At the minimiser \\\hat\theta\\, the minimised value \\Q_N(\hat\theta)\\
is exactly the test statistic \\T_0\\ or \\T_1\\. The system has \\d
\times p\\ moment conditions but only \\m\\ free parameters (\\m = 2d +
p - 2\\ for \\T_1\\; \\m = d + p - 1\\ for \\T_0\\). The excess \\d
\times p - m\\ conditions provide the overidentification test: if the
model is correct, \\Q_N(\hat\theta)\\ is asymptotically \\\chi^2\\ with
\\d \times p - m\\ degrees of freedom.

### From estimation to standard errors

The same GMM framework also provides standard errors for \\\hat\theta\\.
These come from asking: how sensitive is \\\hat\theta\\ to perturbations
in the data? The answer involves two quantities:

- **The Jacobian** \\G\\: how the moment conditions \\E\[U_k\]\\ respond
  to changes in \\\theta\\. This captures the “curvature” of the
  objective function.
- **The variance** \\\Sigma\\: how noisy the individual moment
  contributions \\U_k\\ are.

The next section derives the sandwich variance formula from these two
ingredients.

## Sandwich Standard Errors

### The general sandwich formula

For a GMM estimator with weight matrix \\W\\, the asymptotic variance is
(Newey & McFadden, 1994, Theorem 3.4):

\\ \mathrm{Var}(\hat\theta) = \frac{1}{N} \left(G^\top W G\right)^{-1}
G^\top W\\\Sigma\\W G \left(G^\top W G\right)^{-1} \\

where \\G = E\\\left\[\partial U_k(\theta) / \partial
\theta^\top\right\]\\ is the \\q \times m\\ Jacobian (\\q = d \cdot p\\
moments, \\m\\ free parameters) and \\\Sigma = \mathrm{Var}(U_k)\\ is
the \\q \times q\\ variance of per-observation moment contributions.

### Simplification under efficient GMM

Our estimator uses the efficient weight matrix \\W = \Sigma^{-1}\\
(Hansen, 1982). Substituting into the general formula, the middle term
collapses:

\\ (G^\top \Sigma^{-1} G)^{-1}\\ G^\top
\underbrace{\Sigma^{-1}\\\Sigma\\\Sigma^{-1}}\_{\Sigma^{-1}}\\G\\
(G^\top \Sigma^{-1} G)^{-1} = (G^\top \Sigma^{-1} G)^{-1} \\

giving the efficient sandwich variance:

\\ \mathrm{Var}(\hat\theta) = \frac{1}{N}
\left(G^\top\\\Sigma^{-1}\\G\right)^{-1} \tag{6} \\

This is what
[`fit_structural()`](https://felipelfv.github.io/structest/reference/fit_structural.md)
computes. The same sandwich formula is used in the broader SEM
literature—for instance, Bollen, Kolenikov, & Bauldry (2014) apply it to
general latent variable models via model-implied instrumental variables.
Our setting is simpler: the moment conditions come directly from the
restricted conditional means, so \\G\\ has a closed-form expression and
no instrumental variable selection is required.

The key computational task is deriving \\G\\ analytically and estimating
\\\Sigma\\ from the data.

### Jacobian for the \\T_1\\ model

Under the \\T_1\\ parameterisation, the moment condition for indicator
\\i\\ and \\Z\\-level \\j\\ at observation \\k\\ is

\\ U\_{k,(i,j)} = I(Z_k = z_j)\\\left( X\_{ik} - \gamma_i -
\alpha_i\\\beta_j\right). \\

The free parameter vector is \\\theta = (\gamma_1, \ldots, \gamma_d,\\
\alpha_2, \ldots, \alpha_d,\\ \beta_2, \ldots, \beta_p)\\, with
\\\alpha_1 = 1\\ and \\\beta_1 = 0\\ fixed. The Jacobian \\G\\ has
entries \\G\_{(i,j),\\\ell} = E\[\partial U\_{k,(i,j)} / \partial
\theta\_\ell\]\\. Let \\\pi_j = P(Z = z_j)\\.

**Derivative with respect to \\\gamma_s\\** (\\s = 1, \ldots, d\\):

\\ E\\\left\[\frac{\partial U\_{k,(i,j)}}{\partial\\\gamma_s}\right\] =
-\pi_j\\\mathbf{1}(i = s). \\

Only the moment conditions for indicator \\s\\ are affected.

**Derivative with respect to \\\alpha_s\\** (\\s = 2, \ldots, d\\):

\\ E\\\left\[\frac{\partial U\_{k,(i,j)}}{\partial\\\alpha_s}\right\] =
-\pi_j\\\beta_j\\\mathbf{1}(i = s). \\

This involves the product with \\\beta_j\\: when the \\Z\\-effect
\\\beta_j\\ is near zero, the moment conditions carry little information
about \\\alpha_s\\, and its standard error will be large.

**Derivative with respect to \\\beta_t\\** (\\t = 2, \ldots, p\\):

\\ E\\\left\[\frac{\partial U\_{k,(i,j)}}{\partial\\\beta_t}\right\] =
-\pi_j\\\alpha_i\\\mathbf{1}(j = t). \\

Indicators with larger loadings (\\\alpha_i\\) contribute more
information about the \\Z\\-effect.

**Assembling \\G\\.** The Jacobian is a \\(d \cdot p) \times (2d + p -
2)\\ matrix. Row \\(i, j)\\ corresponds to the moment condition for
indicator \\i\\ at \\Z\\-level \\j\\. In practice, \\\pi_j\\ is
estimated by the sample proportion \\\hat\pi_j = n_j / N\\.

### Jacobian for the \\T_0\\ model

Under the \\T_0\\ parameterisation, the moment condition is

\\ U\_{k,(i,j)} = I(Z_k = z_j)\\\left( X\_{ik} - \gamma_i -
r_i\\\beta_j\right) \\

where \\r_i = \lambda_i / \lambda_1\\ is the pre-estimated loading
ratio, treated as a known constant. The free parameter vector is
\\\theta = (\gamma_1, \ldots, \gamma_d,\\ \beta_2, \ldots, \beta_p)\\.

**Derivative with respect to \\\gamma_s\\:** Identical to the \\T_1\\
case: \\-\pi_j\\\mathbf{1}(i = s)\\.

**Derivative with respect to \\\beta_t\\** (\\t = 2, \ldots, p\\):

\\ E\\\left\[\frac{\partial U\_{k,(i,j)}}{\partial\\\beta_t}\right\] =
-\pi_j\\r_i\\\mathbf{1}(j = t). \\

The Jacobian is a \\(d \cdot p) \times (d + p - 1)\\ matrix—smaller than
the \\T_1\\ Jacobian because there are no \\\alpha\\ parameters. Since
\\r_i\\ is treated as fixed, the Jacobian does not capture uncertainty
in the reliability estimates; that uncertainty is instead absorbed into
\\\Sigma\\ through the variance adjustment described in the *Getting
Started* vignette (see *Variance adjustment for estimated
reliabilities*).

## Example: Fitting the \\T_1\\ Model

``` r
set.seed(12345)
n <- 1000

# Simulate under the structural model
z <- sample(0:2, n, replace = TRUE)
eta <- 1 + 0.3 * (z == 1) + 0.7 * (z == 2) + rnorm(n)
lambda <- c(1.0, 0.8, 0.6)
X <- cbind(
  2 + lambda[1] * eta + rnorm(n, sd = 0.5),
  3 + lambda[2] * eta + rnorm(n, sd = 0.5),
  1 + lambda[3] * eta + rnorm(n, sd = 0.5)
)

fit1 <- fit_structural(X, z, type = "t1")
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
```

The coefficient table reports:

- **gamma**: the intercept for each indicator, equal to \\E(X_i \mid Z =
  z_1)\\ at the reference level. These combine the indicator mean
  \\\mu_i\\ and the baseline latent factor level \\\lambda_i E(\eta \mid
  Z = z_1)\\.
- **alpha**: the relative loading of each indicator, with \\\alpha_1 =
  1\\ fixed for identification. Under the structural model, \\\alpha_i =
  \lambda_i / \lambda_1\\.
- **beta**: the shift in the latent factor at each \\Z\\-level relative
  to the reference, with \\\beta_1 = 0\\.

The overidentification test (J statistic) is equivalent to the \\T_1\\
test from
[`test_t1()`](https://felipelfv.github.io/structest/reference/test_t1.md).
A large \\p\\-value indicates that the structural model fits the data
adequately.

## Example: Fitting the \\T_0\\ Model

``` r
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
```

The \\T_0\\ model estimates fewer parameters (no \\\alpha_i\\), since
the loading ratios are pre-estimated from pairwise covariances. The
estimated reliabilities \\\hat\lambda\\ are reported in the summary
output. The standard errors account for the uncertainty in these
reliability estimates.

## Extracting Results

The `structest_fit` object supports standard R model-object methods:

``` r
# Point estimates (named vector of free parameters)
coef(fit1)
#> gamma[X1] gamma[X2] gamma[X3] alpha[X2] alpha[X3] beta[z=1] beta[z=2] 
#> 2.9461134 3.7879127 1.6110381 0.7226236 0.5235738 0.2918066 0.7395983

# Confidence intervals (default: 95%)
confint(fit1)
#>               2.5 %    97.5 %
#> gamma[X1] 2.8276406 3.0645863
#> gamma[X2] 3.6895664 3.8862590
#> gamma[X3] 1.5319627 1.6901134
#> alpha[X2] 0.6002534 0.8449939
#> alpha[X3] 0.4073157 0.6398319
#> beta[z=1] 0.1261550 0.4574583
#> beta[z=2] 0.5699997 0.9091969

# 90% confidence intervals
confint(fit1, level = 0.90)
#>                 5 %      95 %
#> gamma[X1] 2.8466879 3.0455390
#> gamma[X2] 3.7053779 3.8704476
#> gamma[X3] 1.5446760 1.6774002
#> alpha[X2] 0.6199273 0.8253200
#> alpha[X3] 0.4260070 0.6211407
#> beta[z=1] 0.1527874 0.4308259
#> beta[z=2] 0.5972667 0.8819300

# Variance-covariance matrix
vcov(fit1)
#>               gamma[X1]     gamma[X2]     gamma[X3]     alpha[X2]     alpha[X3]
#> gamma[X1]  0.0036537722  0.0024985253  0.0017794394  0.0002331679  0.0002173940
#> gamma[X2]  0.0024985253  0.0025177922  0.0014479055 -0.0011669238 -0.0002347877
#> gamma[X3]  0.0017794394  0.0014479055  0.0016277438 -0.0002711279 -0.0011252933
#> alpha[X2]  0.0002331679 -0.0011669238 -0.0002711279  0.0038981227  0.0010610621
#> alpha[X3]  0.0002173940 -0.0002347877 -0.0011252933  0.0010610621  0.0035184414
#> beta[z=1] -0.0036202036 -0.0026325650 -0.0018532008 -0.0001561185 -0.0001632195
#> beta[z=2] -0.0036666453 -0.0024502008 -0.0017480050 -0.0007357683 -0.0005342247
#>               beta[z=1]     beta[z=2]
#> gamma[X1] -0.0036202036 -0.0036666453
#> gamma[X2] -0.0026325650 -0.0024502008
#> gamma[X3] -0.0018532008 -0.0017480050
#> alpha[X2] -0.0001561185 -0.0007357683
#> alpha[X3] -0.0001632195 -0.0005342247
#> beta[z=1]  0.0071432404  0.0036698175
#> beta[z=2]  0.0036698175  0.0074877012
```

The [`confint()`](https://rdrr.io/r/stats/confint.html) method accepts a
`parm` argument for selecting specific parameters:

``` r
# CI for beta parameters only
confint(fit1, parm = c("beta[z=1]", "beta[z=2]"))
#>               2.5 %    97.5 %
#> beta[z=1] 0.1261550 0.4574583
#> beta[z=2] 0.5699997 0.9091969
```

## Relationship to the Tests

The
[`fit_structural()`](https://felipelfv.github.io/structest/reference/fit_structural.md)
function and the
[`test_t0()`](https://felipelfv.github.io/structest/reference/test_t0.md)/[`test_t1()`](https://felipelfv.github.io/structest/reference/test_t1.md)
functions share the same underlying GMM estimation. The point estimates
are identical:

``` r
t1_test <- test_t1(X, z)
fit1 <- fit_structural(X, z, type = "t1")

# Same point estimates
all.equal(t1_test$estimates$gamma, fit1$estimates$gamma)
#> [1] TRUE
all.equal(t1_test$estimates$alpha, fit1$estimates$alpha)
#> [1] TRUE
all.equal(t1_test$estimates$beta, fit1$estimates$beta)
#> [1] TRUE

# Same test statistic
all.equal(unname(t1_test$statistic), fit1$test$statistic)
#> [1] TRUE
```

The point estimates and test statistic are already available from the
test functions. The difference is that
[`fit_structural()`](https://felipelfv.github.io/structest/reference/fit_structural.md)
additionally computes the Jacobian \\G\\ and the sandwich variance
(Equation 6), providing standard errors, confidence intervals, and
per-parameter significance tests that the test functions do not.

## References

Bollen, K. A., Kolenikov, S., & Bauldry, S. (2014). Model-implied
instrumental variable-generalized method of moments (MIIV-GMM)
estimators for latent variable models. *Psychometrika*, *79*(1), 20–50.
<https://doi.org/10.1007/s11336-013-9335-3>

Hansen, L. P. (1982). Large sample properties of generalized method of
moments estimators. *Econometrica*, *50*(4), 1029–1054.
<https://doi.org/10.2307/1912775>

Newey, W. K., & McFadden, D. (1994). Large sample estimation and
hypothesis testing. In R. F. Engle & D. L. McFadden (Eds.), *Handbook of
Econometrics* (Vol. 4, pp. 2111–2245). Elsevier.
<https://doi.org/10.1016/S1573-4412(05)80005-4>

VanderWeele, T. J., & Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
*84*(5), 2032–2054. <https://doi.org/10.1111/rssb.12555>
