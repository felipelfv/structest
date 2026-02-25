# Contributing to structest

Thank you for considering contributing to structest! This package
implements the statistical tests from VanderWeele and Vansteelandt
(2022) for rejecting the structural interpretation of a latent factor
model. Contributions that improve the package are welcome.

## Reporting Issues

If you encounter a bug or have a suggestion, please [open an
issue](https://github.com/felipelfv/structest/issues) on GitHub. When
reporting a bug:

- Provide a clear and descriptive title.
- Include a minimal reproducible example demonstrating the problem.
- Describe the expected behaviour and what actually happened.
- Include your R version and
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) output if
  relevant.

## Making Changes

### Getting started

1.  Fork the repository on GitHub.

2.  Clone your fork to your local machine.

    ``` bash
    git clone https://github.com/your-username/structest.git
    cd structest
    ```

3.  Install development dependencies.

    ``` r
    install.packages("devtools")
    devtools::install_deps()
    ```

4.  Make sure the package passes checks before you start.

    ``` r
    devtools::check()
    devtools::test()
    ```

### Creating a branch

Always create a new branch for your work.

``` bash
git checkout -b your-branch-name
```

### Making your changes

- Make your changes in the codebase.
- Write tests for your changes if applicable. We use
  [testthat](https://cran.r-project.org/package=testthat) for unit
  tests.
- We use [roxygen2](https://cran.r-project.org/package=roxygen2) with
  Markdown syntax for documentation. Edit the `.R` source files, not
  `.Rd` files directly.
- Run `devtools::document()` after changing any roxygen comments.
- Run `devtools::test()` to ensure your changes do not break existing
  functionality.

### Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the
documentation directly using the GitHub web interface, as long as the
changes are made in the source `.R` file (not the generated `.Rd` file).

## Submitting a Pull Request

1.  Push your changes to your fork.

    ``` bash
    git push origin your-branch-name
    ```

2.  Open a pull request on GitHub against the `main` branch.

3.  Provide a clear and descriptive title for your pull request.

4.  Describe the changes you made and why they are necessary.

5.  Reference any related issues (e.g., `Fixes #issue-number`).

6.  Ensure all tests pass and there are no merge conflicts.

## LLM Usage

We are not against the use of LLMs (e.g., ChatGPT, Claude, Copilot) to
assist with contributions. However, if you use an LLM to generate code,
please:

- Provide a clear explanation of the choices made in your code. You
  should be able to justify why the code does what it does.
- Be transparent about it — if the LLM wrote the entirety of your
  contribution, say so in the pull request.
- Review and understand the generated code before submitting. You are
  responsible for the correctness of your contribution.

## Code of Conduct

Please note that the structest project is released with a [Contributor
Code of
Conduct](https://felipelfv.github.io/structest/CODE_OF_CONDUCT.md). By
contributing to this project you agree to abide by its terms.
