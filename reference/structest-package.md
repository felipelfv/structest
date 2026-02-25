# structest: Testing the Structural Interpretation of a Latent Factor Model

Implements two statistical tests from VanderWeele and Vansteelandt
(2022) to reject the structural interpretation of a univariate latent
factor model.

## Main functions

- [`test_t0`](https://felipelfv.github.io/structest/reference/test_t0.md):

  Reliability-dependent test (Section 3.2). Requires d \>= 3 indicators
  and p \>= 2 Z-levels.

- [`test_t1`](https://felipelfv.github.io/structest/reference/test_t1.md):

  Reliability-independent test (Section 3.3). Requires d \>= 2
  indicators and p \>= 3 Z-levels.

- [`estimate_reliability`](https://felipelfv.github.io/structest/reference/estimate_reliability.md):

  Estimate reliability coefficients via quasi-Poisson GLM (Section 3.1).

## References

VanderWeele, T. J. and Vansteelandt, S. (2022). A statistical test to
reject the structural interpretation of a latent factor model. *Journal
of the Royal Statistical Society: Series B (Statistical Methodology)*,
84, 2032–2054.

## See also

Useful links:

- <https://felipelfv.github.io/structest/>

- Report bugs at <https://github.com/felipelfv/structest/issues>

## Author

**Maintainer**: Felipe Fontana Vieira <felipelfv98@gmail.com>
