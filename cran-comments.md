# cran-comments

## Test environments

* local macOS (Apple Silicon), R release
* GitHub Actions (via R-CMD-check workflow): ubuntu-latest (release, devel, oldrel), windows-latest (release), macOS-latest (release)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Comments

Implements the statistical tests of VanderWeele and Vansteelandt (2022)
<doi:10.1111/rssb.12555>. The implementation has been verified numerically
against the authors' published replication code.
