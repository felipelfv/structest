# Equivalence Test for the Structural Interpretation

Inverts the burden of proof of
[`test_t0`](https://felipelfv.github.io/structest/reference/test_t0.md)
/
[`test_t1`](https://felipelfv.github.io/structest/reference/test_t1.md).
Those tests can only *reject* the structural interpretation; a large
p-value never confirms it (absence of evidence is not evidence of
absence). This equivalence test asks the complementary question: are the
deviations from the structural model small enough to be considered
negligible, relative to a pre-specified margin?

## Usage

``` r
equiv_test(object, margin, alpha = 0.05)
```

## Arguments

- object:

  a fitted `"structest"` object, as returned by
  [`test_t0`](https://felipelfv.github.io/structest/reference/test_t0.md)
  or
  [`test_t1`](https://felipelfv.github.io/structest/reference/test_t1.md).

- margin:

  numeric scalar, the equivalence margin on the standardised misfit
  \\\varepsilon\\ (an RMSEA-style index, see Details). There is no
  default: the margin is a substantive judgement about how much
  departure from the structural interpretation is tolerable, and must be
  supplied explicitly.

- alpha:

  numeric significance level for the equivalence test (default 0.05).

## Value

An object of class `c("structest_equiv", "structest", "htest")`
containing:

- statistic:

  the underlying T0/T1 distance-metric statistic \\J\\.

- parameter:

  degrees of freedom \\df\\.

- p.value:

  equivalence p-value (lower-tail probability of a non-central
  \\\chi^2\\; small values reject poor fit).

- margin:

  the supplied equivalence margin \\\varepsilon_0\\.

- epsilon:

  point estimate of the standardised misfit \\\hat\varepsilon\\.

- epsilon_upper:

  one-sided \\(1-\alpha)\\ upper confidence limit on \\\varepsilon\\.
  Equivalence holds iff this is below `margin`.

- equivalence:

  logical; `TRUE` if poor fit is rejected at level `alpha`.

- method, data.name, n_obs, d, p:

  carried over from `object`.

## Details

The T0/T1 distance-metric statistic \\J\\ is asymptotically
\\\chi^2\_{df}\\ under the null that the structural interpretation
holds, and non-central \\\chi^2\_{df}(\lambda)\\ under local departures
from it (Newey, 1985). The non-centrality \\\lambda\\ indexes population
misfit. We summarise it on an RMSEA-style scale, \$\$\varepsilon =
\sqrt{\lambda / (n\\df)},\$\$ which is roughly invariant to sample size.
The margin \\\varepsilon_0\\ corresponds to a boundary non-centrality
\\\lambda_0 = n\\df\\\varepsilon_0^2\\.

The test is a one-sided ("close-fit") equivalence test: \$\$H_0:\\
\varepsilon \ge \varepsilon_0 \quad\text{vs}\quad H_1:\\ \varepsilon \<
\varepsilon_0.\$\$ Rejecting \\H_0\\ licenses the positive claim that
departures from the structural interpretation are smaller than the
margin (so the
[`fit_structural`](https://felipelfv.github.io/structest/reference/fit_structural.md)
estimates may be read structurally). The equivalence p-value is
\\P(\chi^2\_{df}(\lambda_0) \le J\_{\text{obs}})\\, and the decision is
equivalent to checking whether `epsilon_upper` falls below `margin`.

Note: unlike SEM, there is no canonical \\\varepsilon\\ threshold for
"approximate structural interpretation". The margin is yours to justify,
and the non-central \\\chi^2\\ approximation is asymptotic; calibration
by simulation under known local misfit is advisable.

## References

VanderWeele, T. J., & Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
*84*(5), 2032–2054.
[doi:10.1111/rssb.12555](https://doi.org/10.1111/rssb.12555)

Newey, W. K. (1985). Generalized method of moments specification
testing. *Journal of Econometrics*, *29*(3), 229–256.

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods*, *1*(2), 130–149.

## See also

[`test_t0`](https://felipelfv.github.io/structest/reference/test_t0.md),
[`test_t1`](https://felipelfv.github.io/structest/reference/test_t1.md)
for the rejection tests this inverts;
[`fit_structural`](https://felipelfv.github.io/structest/reference/fit_structural.md)
for the structural estimates.

## Examples

``` r
set.seed(12345)
n <- 1000
z <- rbinom(n, 1, 0.5)
eta <- 1 + 0.5 * z + rnorm(n)
lambda <- c(1.0, 0.8, 0.6)
X <- cbind(
  2 + lambda[1] * eta + rnorm(n, sd = 0.5),
  3 + lambda[2] * eta + rnorm(n, sd = 0.5),
  1 + lambda[3] * eta + rnorm(n, sd = 0.5)
)

fit <- test_t0(X, z)
equiv_test(fit, margin = 0.05)
#> 
#>   Equivalence test for the structural interpretation (close-fit, based on T0) 
#> 
#> data:   X and z 
#> H0: misfit epsilon >= 0.05  vs  H1: epsilon < 0.05
#> T0 = 6.452, df = 2, equivalence p-value = 0.5349
#> epsilon (estimate) = 0.04718, 95% upper limit = 0.09005
#> 
#> Equivalence NOT supported at alpha = 0.05: cannot rule out misfit at or beyond the margin.
#> 
```
